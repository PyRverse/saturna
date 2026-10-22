# SPDX-License-Identifier: GPL-3.0-or-later
# Original PyRverse implementation of truncated-normal likelihood estimation.
# No locfdr source is copied or translated. See NOTICE.md for the GPL-2 boundary.
"""Empirical-null calibration using an independent likelihood optimizer."""

from dataclasses import dataclass

import numpy as np
from numpy.typing import ArrayLike, NDArray
from scipy.optimize import minimize
from scipy.special import ndtr
from scipy.stats import norm


class EmpiricalNullError(ValueError):
    """The central null distribution cannot be estimated reliably."""


@dataclass(frozen=True)
class EmpiricalNull:
    location: float
    scale: float
    observations: int
    method: str = "two-stage truncated-normal maximum likelihood (SciPy)"


def _truncated_normal_mle(
    z: NDArray, lower: float, upper: float, start: tuple[float, float]
) -> tuple[float, float]:
    central = z[(z >= lower) & (z <= upper)]
    if central.size < 5 or np.ptp(central) < 1e-8:
        raise EmpiricalNullError("Too few distinct central statistics to fit the empirical null")
    average = float(np.mean(central))
    mean_square = float(np.mean(central**2))

    def objective(parameters: NDArray) -> tuple[float, NDArray]:
        location, log_scale = parameters
        scale = np.exp(log_scale)
        a, b = (lower - location) / scale, (upper - location) / scale
        mass = ndtr(b) - ndtr(a)
        if mass <= 1e-15 or not np.isfinite(mass):
            return 1e100, np.zeros(2)
        pa, pb = norm.pdf(a), norm.pdf(b)
        second = mean_square - 2 * location * average + location**2
        value = log_scale + np.log(mass) + second / (2 * scale**2)
        gradient = np.array([
            (location - average) / scale**2 + (pa - pb) / (scale * mass),
            1 - second / scale**2 + (a * pa - b * pb) / mass,
        ])
        return float(value), gradient

    fit = minimize(
        objective, np.array([start[0], np.log(start[1])]), jac=True, method="L-BFGS-B",
        bounds=[(-20, 20), (np.log(0.02), np.log(20))],
        options={"maxiter": 1000, "ftol": 1e-14, "gtol": 1e-10, "maxls": 50},
    )
    location, log_scale = fit.x
    if (not fit.success or not np.isfinite(fit.fun) or fit.fun > 1e50
            or abs(location) > 19.9 or not np.log(0.021) < log_scale < np.log(19.9)):
        raise EmpiricalNullError("Empirical-null likelihood is ill-conditioned or did not converge")
    return float(location), float(np.exp(log_scale))


def estimate_empirical_null(pvalues: ArrayLike, statistics: ArrayLike) -> EmpiricalNull:
    """Estimate the same central-null likelihood target as satuRn.

    This is not a translation of locfdr's Newton iterations. It uses SciPy's
    optimizer; fit diagnostics and R fixture tolerances are reported separately.
    """
    p = np.asarray(pvalues, dtype=float)
    t = np.asarray(statistics, dtype=float)
    if p.ndim != 1 or p.shape != t.shape:
        raise ValueError("pvalues and statistics must be equally sized vectors")
    valid = np.isfinite(p) & ~np.isnan(t)
    if np.any((p[valid] < 0) | (p[valid] > 1)):
        raise ValueError("pvalues must lie between zero and one")
    if valid.sum() < 500:
        raise EmpiricalNullError("satuRn's empirical correction requires at least 500 valid tests")
    z = norm.ppf(p[valid] / 2) * np.sign(t[valid])
    z = z[np.isfinite(z) & (np.abs(z) < 10)]
    if z.size < 20:
        raise EmpiricalNullError("Too few finite central statistics after extreme-tail filtering")
    width_factor = 4.3 * np.exp(-0.26 * np.log10(z.size))
    median = float(np.median(z))
    scale_start = float(np.diff(np.quantile(z, [0.25, 0.75]))[0] / (2 * norm.ppf(0.75)))
    if not np.isfinite(scale_start) or scale_start <= 1e-8:
        raise EmpiricalNullError("Empirical-null statistics have no usable central spread")
    location, scale = _truncated_normal_mle(
        z, median - width_factor * scale_start, median + width_factor * scale_start,
        (median, scale_start),
    )
    location, scale = _truncated_normal_mle(
        z, location - width_factor * scale, location + width_factor * scale, (location, scale),
    )
    return EmpiricalNull(location, scale, int(z.size))


def empirical_pvalues(
    pvalues: ArrayLike, statistics: ArrayLike, fitted_null: EmpiricalNull
) -> NDArray[np.float64]:
    p, t = np.asarray(pvalues, dtype=float), np.asarray(statistics, dtype=float)
    with np.errstate(invalid="ignore"):
        z = norm.ppf(p / 2) * np.sign(t)
        return 2 * norm.sf(np.abs((z - fitted_null.location) / fitted_null.scale))
