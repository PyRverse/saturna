# SPDX-License-Identifier: GPL-3.0-or-later
"""Native bulk differential transcript usage, with reference-based validation."""

from ._api import DTUFit, dtu, fit_dtu, test_dtu
from ._empirical import EmpiricalNullError
from ._glm import ModelFitError

__version__ = "0.1.0a1"
__all__ = ["DTUFit", "EmpiricalNullError", "ModelFitError", "dtu", "fit_dtu", "test_dtu"]
