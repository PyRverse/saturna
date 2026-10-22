# SPDX-License-Identifier: GPL-3.0-or-later
"""Grouped quasi-binomial fitting by iteratively reweighted least squares."""

from dataclasses import dataclass

import numpy as np
from numpy.typing import NDArray
from scipy.linalg import lstsq, qr, solve_triangular
from scipy.special import expit, xlogy


class ModelFitError(ValueError):
    """A transcript model is outside the supported, estimable input contract."""


@dataclass(frozen=True)
class BinomialFit:
    coefficients: NDArray[np.float64]
    covariance: NDArray[np.float64]
    dispersion: float
    residual_df: int
    iterations: int
    information_factor: NDArray[np.float64]


def _deviance(successes: NDArray, total: NDArray, mean: NDArray) -> float:
    failures = total - successes
    return float(2 * np.sum(xlogy(successes, successes / (total * mean))
                           + xlogy(failures, failures / (total * (1 - mean)))))


def fit_quasibinomial(
    successes: NDArray, failures: NDArray, design: NDArray,
    *, max_iter: int = 100, tolerance: float = 1e-8,
) -> BinomialFit:
    """Fit paired successes/failures with logit link and Pearson dispersion.

    The input already includes satuRn's one-count additions to each category.
    The starting mean and relative-deviance criterion match R's GLM conventions.
    """
    total = successes + failures
    y = successes / total
    n, p = design.shape
    if n <= p:
        raise ModelFitError("A gene has insufficient nonzero samples for positive residual df")
    if np.linalg.matrix_rank(design) != p:
        raise ModelFitError("A gene's retained-sample design is rank deficient")
    if np.linalg.cond(design) > 1e8:
        raise ModelFitError("Design is too ill-conditioned for the supported float64 workflow")
    mean = (successes + 0.5) / (total + 1)
    eta = np.log(mean / (1 - mean))
    previous_deviance = _deviance(successes, total, mean)
    for iteration in range(1, max_iter + 1):
        variance = mean * (1 - mean)
        weights = total * variance
        working_response = eta + (y - mean) / variance
        weighted_design = design * np.sqrt(weights)[:, None]
        if np.linalg.cond(weighted_design) > 1e8:
            raise ModelFitError("Weighted design is too ill-conditioned for this workflow")
        beta, _, rank, _ = lstsq(
            weighted_design, working_response * np.sqrt(weights), lapack_driver="gelsy"
        )
        if rank != p:
            raise ModelFitError("Weighted design became rank deficient")
        eta = design @ beta
        mean = expit(eta)
        if np.any((mean <= 0) | (mean >= 1)) or not np.all(np.isfinite(beta)):
            raise ModelFitError("Fitted probabilities reached a numerical boundary")
        deviance = _deviance(successes, total, mean)
        if abs(deviance - previous_deviance) / (0.1 + abs(deviance)) < tolerance:
            # R glm's final stored IRLS weights come from the preceding mean;
            # its working residuals use the final fitted mean.
            residuals = (y - mean) / (mean * (1 - mean))
            dispersion = float(np.sum(weights * residuals**2) / (n - p))
            _, factor = qr(weighted_design, mode="economic")
            inverse_factor = solve_triangular(factor, np.eye(p))
            covariance = inverse_factor @ inverse_factor.T
            return BinomialFit(beta, covariance, dispersion, n - p, iteration, factor)
        previous_deviance = deviance
    raise ModelFitError(f"Quasi-binomial fitting did not converge within {max_iter} iterations")
