import numpy as np
import pandas as pd
import pytest
from numpy.testing import assert_allclose

from saturna import EmpiricalNullError, ModelFitError, dtu, fit_dtu, test_dtu
from saturna._api import _bh
from saturna._empirical import estimate_empirical_null
from saturna._moderation import moderate_variances


@pytest.fixture
def small():
    return (
        np.array([[10, 15, 12, 20, 35, 31], [30, 25, 31, 20, 15, 20],
                  [22, 26, 28, 17, 19, 15], [20, 19, 22, 28, 24, 30]]),
        ["a", "a", "b", "b"],
        np.column_stack([np.ones(6), [0, 0, 0, 1, 1, 1]]),
    )


def test_contrast_sign_and_convenience(small):
    fit = fit_dtu(*small)
    plus = test_dtu(fit, [0, 1], empirical=False)
    minus = test_dtu(fit, [0, -1], empirical=False)
    assert_allclose(plus.estimate, -minus.estimate)
    assert_allclose(plus.pvalue, minus.pvalue)
    pd.testing.assert_frame_equal(plus, dtu(*small, [0, 1], empirical=False))
    assert plus.estimate.iloc[0] > 0
    assert np.all((plus.pvalue >= 0) & (plus.pvalue <= 1))


def test_small_feature_empirical_policy(small):
    fit = fit_dtu(*small)
    with pytest.warns(UserWarning, match="Fewer than 500"):
        result = test_dtu(fit, [0, 1])
    assert result.empirical_pvalue.isna().all()
    with pytest.raises(EmpiricalNullError, match="500"):
        test_dtu(fit, [0, 1], empirical=True)


@pytest.mark.parametrize("bad", [-1, np.nan, np.inf, True, "2", 2**53 + 1])
def test_invalid_counts(small, bad):
    counts, genes, design = small
    mutated = counts.astype(object)
    mutated[0, 0] = bad
    with pytest.raises(ValueError):
        fit_dtu(mutated, genes, design)


def test_mixed_list_precision_rejected(small):
    counts, genes, design = small
    mutated = counts.astype(object).tolist()
    mutated[0][0], mutated[0][1] = 2**53 + 1, 1.5
    with pytest.raises(ValueError, match="precision"):
        fit_dtu(mutated, genes, design)


def test_zero_gene_and_one_isoform_rejected(small):
    counts, genes, design = small
    with pytest.raises(ValueError, match="only one transcript"):
        fit_dtu(counts, ["a", "b", "c", "d"], design)
    counts[:2] = 0
    with pytest.raises(ModelFitError, match="insufficient samples"):
        fit_dtu(counts, genes, design)


def test_gene_zero_samples_can_destroy_design_rank(small):
    counts, genes, design = small
    counts[:2, 3:] = 0
    with pytest.raises(ModelFitError, match="rank deficient"):
        fit_dtu(counts, genes, design)


def test_alignment_and_labels(small):
    counts, genes, design = small
    c = pd.DataFrame(counts, index=list("abcd"), columns=list("uvwxyz"))
    g = pd.Series(genes, index=c.index)
    x = pd.DataFrame(design, index=c.columns, columns=["intercept", "condition"])
    fit = fit_dtu(c, g, x)
    result = test_dtu(fit, pd.Series([0, 1], index=x.columns), empirical=False)
    assert list(result.index) == list("abcd")
    with pytest.raises(ValueError, match="sample labels"):
        fit_dtu(c, g, x.iloc[::-1])
    with pytest.raises(ValueError, match="Gene mapping index"):
        fit_dtu(c, g.iloc[::-1], x)
    with pytest.raises(ValueError, match="Contrast index"):
        test_dtu(fit, pd.Series([1, 0], index=x.columns[::-1]))


@pytest.mark.parametrize("contrast", [[0, 0], [1], [0, np.nan], [0, np.inf]])
def test_invalid_contrast(small, contrast):
    with pytest.raises(ValueError):
        test_dtu(fit_dtu(*small), contrast, empirical=False)


def test_nonconvergence_and_design_validation(small):
    counts, genes, design = small
    with pytest.raises(ModelFitError, match="did not converge"):
        fit_dtu(*small, max_iter=1)
    with pytest.raises(ValueError, match="full column rank"):
        fit_dtu(counts, genes, np.ones((6, 2)))
    with pytest.raises(ValueError, match="finite"):
        fit_dtu(counts, genes, design * np.nan)


def test_all_constant_counts_are_rejected_as_degenerate():
    counts = np.tile([[10], [20], [30], [40]], (1, 8))
    design = np.column_stack([np.ones(8), np.repeat([0, 1], 4)])
    with pytest.raises(ModelFitError, match="effectively zero"):
        fit_dtu(counts, ["a", "a", "b", "b"], design)


def test_two_transcripts_remain_unmoderated(small):
    counts, genes, design = small
    fit = fit_dtu(counts[:2], genes[:2], design)
    assert fit.moderation.method == "unmoderated"
    assert_allclose(fit.variance_posterior, fit.dispersion)
    assert_allclose(fit.df_posterior, fit.residual_df)


def test_bh_ties_and_missing_values():
    assert_allclose(_bh(np.array([0.01, 0.04, 0.04, np.nan, 0.8])),
                    [0.04, 0.053333333333, 0.053333333333, np.nan, 0.8], equal_nan=True)


@pytest.mark.parametrize("values,degrees", [([1, -1], 3), ([1, np.nan], 3), ([1, 2], 0),
                                          ([0, 0, 0], [3, 4, 5])])
def test_invalid_moderation(values, degrees):
    with pytest.raises(ValueError):
        moderate_variances(values, degrees)


def test_degenerate_empirical_null_rejected():
    with pytest.raises(EmpiricalNullError, match="spread"):
        estimate_empirical_null(np.ones(600), np.zeros(600))


def test_invalid_empirical_setting(small):
    with pytest.raises(ValueError, match="empirical must"):
        test_dtu(fit_dtu(*small), [0, 1], empirical=1)


@pytest.mark.parametrize("scale", [1e-200, 1e200, -1e-200, -1e200])
def test_contrast_scale_preserves_inference(small, scale):
    fit = fit_dtu(*small)
    expected = test_dtu(fit, [0, 1], empirical=False)
    actual = test_dtu(fit, [0, scale], empirical=False)
    assert_allclose(actual.pvalue, expected.pvalue, atol=0, rtol=1e-14)
    assert_allclose(actual.statistic, np.sign(scale) * expected.statistic, rtol=1e-14)
    assert_allclose(actual.estimate / scale, expected.estimate, rtol=1e-14)


def test_correlated_design_preserves_complementary_transcript_tests():
    rng = np.random.default_rng(1)
    counts = rng.poisson(100, (4, 20))
    group = np.repeat([0, 1], 10)
    design = np.column_stack([np.ones(20), group, group + 1e-7 * rng.normal(size=20)])
    fit = fit_dtu(counts, ["a", "a", "b", "b"], design)
    result = test_dtu(fit, [0, 1, 1], empirical=False)
    assert_allclose(result.pvalue.iloc[::2], result.pvalue.iloc[1::2], rtol=2e-7, atol=1e-9)
    assert_allclose(result.standard_error.iloc[::2], result.standard_error.iloc[1::2], rtol=2e-7)
