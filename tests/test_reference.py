"""Compare native calculations with committed outputs of genuine R execution."""

import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
from numpy.testing import assert_allclose

from saturna import fit_dtu, test_dtu
from saturna._empirical import empirical_pvalues, estimate_empirical_null
from saturna._moderation import moderate_variances

ROOT = Path(__file__).resolve().parents[1]


def read_reference(name):
    return json.loads((ROOT / "reference/outputs" / f"{name}.json").read_text())


def input_data(name):
    counts = pd.read_csv(ROOT / "datasets" / name / "counts.csv", index_col=0)
    genes = counts.pop("gene_id")
    design = pd.read_csv(ROOT / "datasets" / name / "design.csv", index_col=0)
    return counts, genes, design


@pytest.mark.parametrize("scenario", ["balanced", "unequal_df", "empirical", "covariate"])
def test_full_r_pipeline(scenario):
    reference = read_reference(scenario)
    fitted = fit_dtu(*input_data(scenario))
    for key in ("coefficients", "covariance_unscaled", "dispersion", "residual_df",
                "variance_posterior", "df_posterior"):
        assert_allclose(getattr(fitted, key), reference[key], rtol=2e-8, atol=2e-10)
    result = test_dtu(fitted, reference["contrast"], empirical=scenario == "empirical")
    assert list(result.index) == reference["transcript_ids"]
    columns = {"estimate": "estimates", "standard_error": "se", "df": "df",
               "statistic": "t", "pvalue": "pval", "fdr": "regular_FDR"}
    for native, upstream in columns.items():
        assert_allclose(result[native], reference["result"][upstream], rtol=2e-8, atol=2e-10)
    if scenario == "empirical":
        for native, upstream in (("empirical_pvalue", "empirical_pval"),
                                 ("empirical_fdr", "empirical_FDR")):
            assert_allclose(result[native], reference["result"][upstream], rtol=3e-6, atol=3e-8)
    if scenario == "unequal_df":
        assert np.unique(fitted.residual_df).size == 3
        assert fitted.moderation.method == "unequal_df"


def test_independent_shifted_empirical_null():
    reference = read_reference("empirical_null")
    fitted = estimate_empirical_null(reference["pvalues"], reference["statistics"])
    corrected = empirical_pvalues(reference["pvalues"], reference["statistics"], fitted)
    assert_allclose(corrected, reference["result"]["pval"], rtol=3e-6, atol=3e-8)
    assert 0.24 < fitted.location < 0.26
    assert 1.29 < fitted.scale < 1.32


@pytest.mark.parametrize("name", ["equal", "infinite_prior", "unequal", "zeros", "two"])
def test_moderation_r_cases(name):
    case = read_reference("moderation")[name]
    fit = moderate_variances(case["input"]["var"], case["input"]["df"])
    reference = case["output"]
    assert_allclose(fit.posterior, reference["var.post"], rtol=1e-8, atol=1e-9)
    assert_allclose(fit.prior_variance, reference["var.prior"], rtol=1e-8, atol=1e-9)
    assert_allclose(fit.prior_df, float(reference["df.prior"]), rtol=1e-8, atol=1e-9)


def test_reference_integrity_and_runtime():
    manifest = json.loads((ROOT / "reference/manifest.json").read_text())
    for path, expected in manifest["sha256"].items():
        assert hashlib.sha256((ROOT / path).read_bytes()).hexdigest() == expected, path
    assert manifest["runtime"]["packages"]["satuRn"] == "1.20.0"
    assert manifest["runtime"]["packages"]["limma"] == "3.68.4"
    assert manifest["generator"].startswith("Actual Rscript")
