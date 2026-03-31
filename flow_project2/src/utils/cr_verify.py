"""Cauchy-Riemann verification helpers."""

from __future__ import annotations

import numpy as np


def verify_cr(
    phi: np.ndarray,
    psi: np.ndarray,
    dx: float,
    dy: float,
    mask: np.ndarray | None = None,
) -> dict[str, float]:
    """功能：验证数值场是否满足 C-R 方程。

    输入：
        phi: 速度势数组。
        psi: 流函数数组。
        dx: x 方向网格间距。
        dy: y 方向网格间距。
        mask: 有效区域掩码（可选）。

    输出：
        C-R 残差统计字典。

    物理参数说明：
        势流模型在解析区域应满足 C-R 条件，残差应接近 0。
    """
    if dx <= 0 or dy <= 0:
        raise ValueError("dx and dy must be positive.")

    dphi_dy, dphi_dx = np.gradient(phi, dy, dx)
    dpsi_dy, dpsi_dx = np.gradient(psi, dy, dx)

    r1 = dphi_dx - dpsi_dy
    r2 = dphi_dy + dpsi_dx

    if mask is not None:
        if mask.shape != phi.shape:
            raise ValueError("mask shape must match phi/psi shape.")
        r1 = np.where(mask, r1, np.nan)
        r2 = np.where(mask, r2, np.nan)

    return {
        "max_abs_r1": float(np.nanmax(np.abs(r1))),
        "max_abs_r2": float(np.nanmax(np.abs(r2))),
        "mean_abs_r1": float(np.nanmean(np.abs(r1))),
        "mean_abs_r2": float(np.nanmean(np.abs(r2))),
    }
