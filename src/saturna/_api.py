# SPDX-License-Identifier: GPL-3.0-or-later
"""Public bulk DTU APIs with explicit design and sample-alignment checks."""

import warnings
from dataclasses import dataclass
from typing import Literal

import numpy as np
import pandas as pd
from numpy.typing import ArrayLike, NDArray
from scipy.linalg import solve_triangular
from scipy.stats import t as student_t

from ._empirical import empirical_pvalues, estimate_empirical_null
from ._glm import ModelFitError, fit_quasibinomial
from ._moderation import ModerationResult, moderate_variances


@dataclass(frozen=True)
class DTUFit:
    """Fitted transcript models; each row corresponds to transcript_ids in order."""

    coefficients: NDArray[np.float64]
    covariance_unscaled: NDArray[np.float64]
    dispersion: NDArray[np.float64]
    residual_df: NDArray[np.float64]
    variance_posterior: NDArray[np.float64]
    df_posterior: NDArray[np.float64]
    iterations: NDArray[np.int64]
    transcript_ids: tuple[str, ...]
    gene_ids: tuple[str, ...]
    coefficient_names: tuple[str, ...]
    moderation: ModerationResult
    information_factor: NDArray[np.float64]


def _counts_array(counts: ArrayLike) -> NDArray[np.float64]:
    original = np.asarray(counts, dtype=object)
    for value in original.flat:
        if isinstance(value, (bool, np.bool_)) or not isinstance(
            value, (int, float, np.integer, np.floating)
        ):
            raise ValueError("counts must contain real numbers, without booleans or strings")
        if isinstance(value, (int, np.integer)):
            try:
                converted = float(value)
            except OverflowError as error:
                raise ValueError("Integer count exceeds float64 precision") from error
            if not np.isfinite(converted) or int(converted) != int(value):
                raise ValueError("Integer count exceeds float64 precision")
    values = original.astype(float)
    if values.ndim != 2 or 0 in values.shape:
        raise ValueError("counts must be a nonempty transcript-by-sample matrix")
    if not np.all(np.isfinite(values)) or np.any(values < 0):
        raise ValueError("counts must be finite and nonnegative; missing values are unsupported")
    return values


def fit_dtu(
    counts: ArrayLike | pd.DataFrame,
    gene_ids: ArrayLike | pd.Series,
    design: ArrayLike | pd.DataFrame,
    *, transcript_ids: ArrayLike | None = None,
    max_iter: int = 100,
    tolerance: float = 1e-8,
) -> DTUFit:
    """Fit quasi-binomial usage models and nonrobust empirical-Bayes moderation.

    Counts are transcripts x independent bulk samples; design is samples x
    coefficients and must include any desired intercept explicitly. Genes must
    contain at least two transcripts. Missing, nonestimable and ill-conditioned
    models are rejected rather than silently receiving inferential results.
    """
    values = _counts_array(counts)
    m, n = values.shape
    if not isinstance(max_iter, int) or isinstance(max_iter, bool) or max_iter < 1:
        raise ValueError("max_iter must be a positive integer")
    if not np.isfinite(tolerance) or not 0 < tolerance < 1:
        raise ValueError("tolerance must be finite and strictly between zero and one")
    if isinstance(counts, pd.DataFrame):
        if counts.index.has_duplicates or counts.columns.has_duplicates:
            raise ValueError("Count transcript and sample labels must be unique")
        if isinstance(design, pd.DataFrame) and not counts.columns.equals(design.index):
            raise ValueError("Design sample labels must equal count columns in the same order")
        if isinstance(gene_ids, pd.Series) and not counts.index.equals(gene_ids.index):
            raise ValueError("Gene mapping index must equal count transcript labels in order")
        if transcript_ids is not None:
            raise ValueError("Use the DataFrame index for transcript IDs")
        transcript_ids = counts.index
    raw_ids = np.asarray(transcript_ids if transcript_ids is not None else range(m), dtype=object)
    if raw_ids.shape != (m,) or pd.isna(raw_ids).any():
        raise ValueError("transcript_ids must contain one nonmissing identifier per row")
    ids = tuple(str(x) for x in raw_ids)
    if len(ids) != m or len(set(ids)) != m:
        raise ValueError("transcript_ids must be unique and match count rows")
    mapping = np.asarray(gene_ids, dtype=object)
    if mapping.shape != (m,) or pd.isna(mapping).any():
        raise ValueError("gene_ids must supply one nonmissing gene identifier per transcript")
    genes = tuple(str(x) for x in mapping)
    if len(pd.unique(mapping)) != len(set(genes)):
        raise ValueError("Distinct gene identifiers must remain distinct when converted to strings")
    matrix = np.asarray(design, dtype=float)
    if matrix.ndim != 2 or matrix.shape[0] != n or matrix.shape[1] < 1:
        raise ValueError("design must be a sample-by-coefficient matrix")
    if not np.all(np.isfinite(matrix)):
        raise ValueError("design must contain finite numeric values")
    p = matrix.shape[1]
    if n <= p or np.linalg.matrix_rank(matrix) != p:
        raise ValueError(
            "design must have full column rank and positive residual degrees of freedom"
        )
    names = tuple(str(x) for x in (
        design.columns if isinstance(design, pd.DataFrame) else range(p)
    ))
    if len(set(names)) != len(names):
        raise ValueError("Design coefficient names must be unique")
    coefficients = np.empty((m, p))
    covariance = np.empty((m, p, p))
    information_factor = np.empty((m, p, p))
    dispersions, degrees = np.empty(m), np.empty(m)
    iterations = np.empty(m, dtype=np.int64)
    for gene in dict.fromkeys(genes):
        rows = np.flatnonzero(np.array(genes) == gene)
        if rows.size < 2:
            raise ValueError(f"Gene {gene!r} has only one transcript; filter single-isoform genes")
        totals = values[rows].sum(axis=0)
        if np.any(totals > 2**52):
            raise ValueError("Gene totals exceed the safe float64 count range (2**52)")
        keep = totals > 0
        if keep.sum() <= p:
            raise ModelFitError(f"Gene {gene!r} has insufficient samples with nonzero gene totals")
        for row in rows:
            successes = values[row, keep] + 1
            failures = totals[keep] - values[row, keep] + 1
            try:
                model = fit_quasibinomial(
                    successes, failures, matrix[keep], max_iter=max_iter, tolerance=tolerance
                )
            except ModelFitError as error:
                raise ModelFitError(f"Transcript {ids[row]!r}: {error}") from error
            if model.dispersion <= 1e-12:
                raise ModelFitError(
                    f"Transcript {ids[row]!r}: residual dispersion is effectively zero; "
                    "a perfectly explained model cannot support this inferential workflow"
                )
            coefficients[row] = model.coefficients
            covariance[row] = model.covariance
            information_factor[row] = model.information_factor
            dispersions[row] = model.dispersion
            degrees[row] = model.residual_df
            iterations[row] = model.iterations
    moderation = moderate_variances(dispersions, degrees)
    return DTUFit(
        coefficients, covariance, dispersions, degrees, moderation.posterior,
        degrees + moderation.prior_df, iterations, ids, genes, names, moderation,
        information_factor,
    )


def _bh(pvalues: NDArray) -> NDArray:
    adjusted = np.full(pvalues.shape, np.nan)
    valid = np.flatnonzero(np.isfinite(pvalues))
    order = valid[np.argsort(pvalues[valid], kind="stable")]
    if order.size:
        ranked = pvalues[order] * order.size / np.arange(1, order.size + 1)
        adjusted[order] = np.minimum(1, np.minimum.accumulate(ranked[::-1])[::-1])
    return adjusted


def test_dtu(
    fitted: DTUFit, contrast: ArrayLike, *, empirical: bool | Literal["auto"] = "auto"
) -> pd.DataFrame:
    """Test one coefficient contrast and return transcript-level statistics.

    `auto` runs empirical calibration for >=500 valid tests, as in satuRn1.20.0.
    `False` explicitly requests only the regular moderated test. `True` requires
    >=500 valid tests. All BH corrections are across this contrast's transcripts.
    """
    if not isinstance(fitted, DTUFit):
        raise TypeError("fitted must be the result of fit_dtu")
    if isinstance(contrast, pd.Series):
        if tuple(str(x) for x in contrast.index) != fitted.coefficient_names:
            raise ValueError("Contrast index must match fitted coefficient names in order")
    vector = np.asarray(contrast, dtype=float)
    if vector.shape != (fitted.coefficients.shape[1],) or not np.all(np.isfinite(vector)):
        raise ValueError("contrast must be a finite vector matching the fitted coefficients")
    if not np.any(vector != 0):
        raise ValueError("An all-zero contrast is not a statistical test")
    if not (type(empirical) is bool or isinstance(empirical, str) and empirical == "auto"):
        raise ValueError("empirical must be True, False, or 'auto'")
    # Normalize before the quadratic form: squaring an otherwise finite contrast
    # can overflow/underflow and change the test. Scale only reported quantities.
    contrast_scale = np.max(np.abs(vector))
    normalized = vector / contrast_scale
    normalized_estimates = fitted.coefficients @ normalized
    # Direct triangular solves avoid cancellation in c' (X'WX)^-1 c when
    # design columns are nearly collinear. The full covariance is for inspection.
    transformed = np.array([
        solve_triangular(factor.T, normalized, lower=True)
        for factor in fitted.information_factor
    ])
    unscaled = np.sum(transformed**2, axis=1)
    normalized_error = np.sqrt(unscaled * fitted.variance_posterior)
    with np.errstate(divide="ignore", invalid="ignore"):
        statistic = normalized_estimates / normalized_error
    with np.errstate(over="ignore", under="ignore"):
        estimates = normalized_estimates * contrast_scale
        standard_error = normalized_error * contrast_scale
    if (not np.all(np.isfinite(estimates)) or not np.all(np.isfinite(standard_error))
            or np.any(standard_error == 0)):
        raise ValueError(
            "Contrast scale cannot be reported safely in float64; rescale the contrast"
        )
    pvalues = 2 * student_t.sf(np.abs(statistic), fitted.df_posterior)
    result = pd.DataFrame({
        "gene_id": fitted.gene_ids,
        "estimate": estimates,
        "standard_error": standard_error,
        "df": fitted.df_posterior,
        "statistic": statistic,
        "pvalue": pvalues,
        "fdr": _bh(pvalues),
        "empirical_pvalue": np.nan,
        "empirical_fdr": np.nan,
    }, index=pd.Index(fitted.transcript_ids, name="transcript_id"))
    eligible = int(np.count_nonzero(np.isfinite(pvalues))) >= 500
    if empirical is True or (empirical == "auto" and eligible):
        null = estimate_empirical_null(pvalues, statistic)
        corrected = empirical_pvalues(pvalues, statistic, null)
        result["empirical_pvalue"] = corrected
        result["empirical_fdr"] = _bh(corrected)
        result.attrs["empirical_null"] = {
            "location": null.location, "scale": null.scale,
            "observations": null.observations, "method": null.method,
        }
    elif empirical == "auto":
        warnings.warn("Fewer than 500 valid tests: empirical correction was not applied",
                      UserWarning, stacklevel=2)
        result.attrs["empirical_null"] = {"applied": False, "reason": "fewer than 500 valid tests"}
    else:
        result.attrs["empirical_null"] = {"applied": False, "reason": "explicitly disabled"}
    result.attrs["effect_scale"] = "natural-log odds ratio for the specified contrast"
    result.attrs["moderation"] = fitted.moderation.method
    return result


test_dtu.__test__ = False


def dtu(
    counts: ArrayLike | pd.DataFrame, gene_ids: ArrayLike | pd.Series,
    design: ArrayLike | pd.DataFrame, contrast: ArrayLike,
    *, empirical: bool | Literal["auto"] = "auto", **fit_options,
) -> pd.DataFrame:
    """Fit and test one DTU contrast; see fit_dtu and test_dtu for supported inputs."""
    return test_dtu(fit_dtu(counts, gene_ids, design, **fit_options), contrast, empirical=empirical)
