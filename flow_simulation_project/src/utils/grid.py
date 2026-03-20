"""Grid generation helpers."""

from __future__ import annotations

import numpy as np


def create_grid(
    x_min: float = -2.0,
    x_max: float = 2.0,
    y_min: float = -2.0,
    y_max: float = 2.0,
    n: int = 200,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Create Cartesian meshgrid and complex coordinates."""
    x = np.linspace(x_min, x_max, n)
    y = np.linspace(y_min, y_max, n)
    X, Y = np.meshgrid(x, y)
    Z = X + 1j * Y
    return X, Y, Z
