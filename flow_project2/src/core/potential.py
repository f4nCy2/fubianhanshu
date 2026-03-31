"""Complex potential definitions for Project 2."""

from __future__ import annotations

import numpy as np


def cylinder_potential(z: np.ndarray, U: float = 1.0, a: float = 1.0) -> np.ndarray:
    """功能：计算无环量圆柱绕流复势。

    输入：
        z: 复平面坐标数组，z = x + iy。
        U: 来流速度，单位 m/s。
        a: 圆柱半径，单位 m。

    输出：
        复势 W0(z) = U*(z + a^2/z)。

    物理参数说明：
        该项对应项目 1 的基准模型，压力分布上下对称。
    """
    if a <= 0:
        raise ValueError("a must be positive.")
    z_safe = np.where(np.abs(z) < 1e-12, 1e-12 + 0j, z)
    return U * (z_safe + (a**2) / z_safe)


def vortex_potential(z: np.ndarray, gamma: float = 0.0) -> np.ndarray:
    """功能：计算点涡项复势。

    输入：
        z: 复平面坐标数组。
        gamma: 环量，单位 m^2/s。

    输出：
        复势 Wv(z) = i*Gamma/(2*pi)*ln(z)。

    物理参数说明：
        环量项用于打破对称性并产生升力。
    """
    z_safe = np.where(np.abs(z) < 1e-12, 1e-12 + 0j, z)
    return 1j * gamma / (2.0 * np.pi) * np.log(z_safe)


def total_potential(z: np.ndarray, U: float = 1.0, a: float = 1.0, gamma: float = 0.0) -> np.ndarray:
    """功能：叠加得到项目 2 总复势。

    输入：
        z: 复平面坐标数组。
        U: 来流速度，单位 m/s。
        a: 圆柱半径，单位 m。
        gamma: 环量，单位 m^2/s。

    输出：
        复势 W(z) = U*(z+a^2/z)+i*Gamma/(2*pi)*ln(z)。

    物理参数说明：
        Gamma 决定驻点迁移和升力大小。
    """
    return cylinder_potential(z=z, U=U, a=a) + vortex_potential(z=z, gamma=gamma)
