"""Grid generation helpers."""

from __future__ import annotations

import numpy as np


def create_grid(
    x_min: float = -3.0,
    x_max: float = 3.0,
    y_min: float = -3.0,
    y_max: float = 3.0,
    n: int = 201,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """功能：创建二维笛卡尔网格及对应复坐标。

    输入：
        x_min, x_max: x 方向范围。
        y_min, y_max: y 方向范围。
        n: 单方向采样点数。

    输出：
        X, Y, Z，其中 Z = X + iY。

    物理参数说明：
        网格范围决定可视化窗口，n 决定数值分辨率。
    """
    if n <= 1:
        raise ValueError("n must be greater than 1.")

    x = np.linspace(x_min, x_max, n)
    y = np.linspace(y_min, y_max, n)
    X, Y = np.meshgrid(x, y)
    Z = X + 1j * Y
    return X, Y, Z
