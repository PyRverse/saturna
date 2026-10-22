# SPDX-License-Identifier: GPL-3.0-or-later
# Adapted from limma 3.68.5: fitFDist.R, fitFDistUnequalDF1.R, squeezeVar.R.
# Original algorithms: Gordon Smyth, Belinda Phipson, Lizhong Chen.
# Upstream license: GPL (>= 2); distributed here under GPL-3.0-or-later.
# Python adaptation and modifications: PyRverse contributors, 2026.
# See NOTICE.md and docs/provenance.md for source revisions and scope.
"""Nonrobust empirical-Bayes dispersion moderation without a covariate trend."""

from dataclasses import dataclass

import numpy as np
from numpy.typing import ArrayLike, NDArray
from scipy.optimize import brentq, minimize_scalar
from scipy.special import digamma, gammaln, polygamma


@dataclass(frozen=True)
class ModerationResult:
    posterior: NDArray[np.float64]
    prior_variance: float | NDArray[np.float64]
    prior_df: float
    method: str


def _log_minus_digamma(value: ArrayLike) -> NDArray[np.float64]:
    x = np.asarray(value, dtype=float)
    # Avoid cancellation at very large shape parameters.
    with np.errstate(invalid="ignore", divide="ignore"):
        ordinary = np.log(x) - digamma(x)
        asymptotic = 1 / (2 * x) + 1 / (12 * x**2) - 1 / (120 * x**4)
    return np.where(x > 1e6, asymptotic, ordinary)


def _inverse_trigamma(value: float) -> float:
    if value > 1e7:
        return float(1 / np.sqrt(value))
    if value < 1e-6:
        return 1 / value
    return float(brentq(lambda x: float(polygamma(1, x)) - value, 1e-12, 1 / value + 1))


def _equal_df_prior(variances: NDArray, df: NDArray) -> tuple[float, float]:
    median = float(np.median(variances))
    if median == 0:
        median = 1.0
    positive = np.maximum(variances, median * 1e-5)
    adjusted_log = np.log(positive) + _log_minus_digamma(df / 2)
    excess = float(np.var(adjusted_log, ddof=1) - np.mean(polygamma(1, df / 2)))
    if excess <= 0:
        return float(np.mean(positive)), np.inf
    prior_df = 2 * _inverse_trigamma(excess)
    scale = np.exp(np.mean(adjusted_log) - _log_minus_digamma(prior_df / 2))
    return float(scale), float(prior_df)


def _unequal_df_prior(variances: NDArray, df: NDArray) -> tuple[float, float]:
    positive_values = variances[variances > 0]
    if positive_values.size < 2:
        raise ValueError("Unequal-df moderation needs at least two positive dispersions")
    positive = np.maximum(variances, 1e-12 * np.median(positive_values))
    d1 = df / 2
    weights = 1 / polygamma(1, d1)
    adjusted = np.log(positive) + _log_minus_digamma(d1)
    mean = float(np.sum(weights * adjusted) / weights.sum())
    d1x = d1 * positive

    def objective(parameter: float) -> float:
        d2 = parameter / (1 - parameter)
        scale = np.exp(mean - _log_minus_digamma(d2))
        product = d2 * scale
        terms = -(d1 + d2) * np.log1p(d1x / product) - d1 * np.log(product)
        terms += gammaln(d1 + d2) - gammaln(d2)
        return float(-2 * terms.sum())

    # R optimize() default tolerance; the parameterization and bounds are limma's.
    fit = minimize_scalar(
        objective, bounds=(0.5, 0.9998), method="bounded",
        options={"xatol": np.finfo(float).eps**0.25},
    )
    if not fit.success or not np.isfinite(fit.fun):
        raise RuntimeError("Unequal-df dispersion moderation did not converge")
    d2 = fit.x / (1 - fit.x)
    return float(np.exp(mean - _log_minus_digamma(d2))), float(2 * d2)


def moderate_variances(variances: ArrayLike, df: ArrayLike) -> ModerationResult:
    """Moderate finite nonnegative variances using limma 3.68.5's nonrobust policy.

    Equal degrees of freedom use moment estimation; unequal positive degrees of
    freedom use the weighted likelihood estimator. Covariate trends, robust
    moderation, missing variances, and nonpositive degrees of freedom are not
    supported. For fewer than three features, return unmoderated variances.
    """
    values = np.asarray(variances, dtype=float)
    degrees = np.broadcast_to(np.asarray(df, dtype=float), values.shape).copy()
    if values.ndim != 1 or not values.size:
        raise ValueError("variances must be a nonempty one-dimensional array")
    if not np.all(np.isfinite(values)) or np.any(values < 0):
        raise ValueError("variances must be finite and nonnegative")
    if not np.all(np.isfinite(degrees)) or np.any(degrees <= 0):
        raise ValueError("degrees of freedom must be finite and positive")
    if values.size < 3:
        return ModerationResult(values.copy(), values.copy(), 0.0, "unmoderated")
    equal = np.all(degrees == degrees[0])
    scale, prior_df = (
        _equal_df_prior(values, degrees) if equal else _unequal_df_prior(values, degrees)
    )
    posterior = (
        np.full(values.shape, scale) if np.isinf(prior_df)
        else (degrees * values + prior_df * scale) / (degrees + prior_df)
    )
    return ModerationResult(posterior, scale, prior_df, "equal_df" if equal else "unequal_df")
