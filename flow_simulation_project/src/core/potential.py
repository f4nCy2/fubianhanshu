"""Complex potential definitions for ideal flow."""

from __future__ import annotations

import numpy as np


def doublet_strength(U: float, a: float) -> float:
    """功能：根据无穿透边界条件计算偶极子强度。

    输入：
        U: 远场来流速度，单位 m/s。
        a: 圆柱半径，单位 m。

    输出：
        kappa: 偶极子强度，满足 kappa = 2*pi*U*a^2。

    物理参数说明：
        U 表征均匀来流强度；a 决定圆柱几何尺度。
    """
    if a <= 0:
        raise ValueError("a must be positive.")
    return 2.0 * np.pi * U * (a**2)


def uniform_flow(z: np.ndarray, U: float = 1.0, alpha: float = 0.0) -> np.ndarray:
    """功能：计算均匀来流的复势 W(z)。

    输入：
        z: 复平面坐标数组，z = x + i y。
        U: 远场来流速度，单位 m/s。
        alpha: 来流与 x 轴夹角，单位 rad。

    输出：
        复势数组 W(z) = U*exp(-i*alpha)*z。

    物理参数说明：
        alpha 控制来流方向，U 控制速度尺度。
    """
    return U * np.exp(-1j * alpha) * z


def doublet_potential(z: np.ndarray, kappa: float) -> np.ndarray:
    """功能：计算原点偶极子的复势。

    输入：
        z: 复平面坐标数组，z = x + i y。
        kappa: 偶极子强度。

    输出：
        偶极子复势数组 W_d = kappa/(2*pi*z)。

    物理参数说明：
        kappa 由圆柱无穿透条件确定，决定扰动势强弱。
    """
    z_safe = np.where(np.abs(z) < 1e-12, 1e-12 + 0j, z)
    return kappa / (2.0 * np.pi * z_safe)


def cylinder_flow(z: np.ndarray, U: float = 1.0, a: float = 1.0, alpha: float = 0.0) -> np.ndarray:
    """功能：计算圆柱绕流（均匀流 + 偶极子）的复势。

    输入：
        z: 复平面坐标数组。
        U: 远场来流速度，单位 m/s。
        a: 圆柱半径，单位 m。
        alpha: 来流角，单位 rad。

    输出：
        圆柱绕流复势 W = U*exp(-i*alpha)*(z + a^2/z)。

    物理参数说明：
        a 体现障碍物几何尺寸；alpha 影响流场旋转方向。
    """
    kappa = doublet_strength(U=U, a=a)
    return uniform_flow(z, U=U, alpha=alpha) + np.exp(-1j * alpha) * doublet_potential(z, kappa=kappa)


def vortex(z: np.ndarray, gamma: float = 1.0) -> np.ndarray:
    """功能：计算原点点涡的复势。

    输入：
        z: 复平面坐标数组。
        gamma: 环量，单位 m^2/s。

    输出：
        点涡复势数组 W_v = -i*gamma/(2*pi)*log(z)。

    物理参数说明：
        gamma 的符号决定旋转方向，绝对值决定涡强。
    """
    z_safe = np.where(np.abs(z) < 1e-12, 1e-12 + 0j, z)
    return -1j * gamma / (2 * np.pi) * np.log(z_safe)


def total_potential(
    z: np.ndarray,
    U: float = 1.0,
    a: float = 1.0,
    alpha: float = 0.0,
    gamma: float = 0.0,
) -> np.ndarray:
    """功能：叠加圆柱绕流与环量项，得到总复势。

    输入：
        z: 复平面坐标数组。
        U: 远场来流速度，单位 m/s。
        a: 圆柱半径，单位 m。
        alpha: 来流角，单位 rad。
        gamma: 环量，单位 m^2/s。

    输出：
        总复势数组 W_total(z)。

    物理参数说明：
        gamma != 0 时可模拟升力相关的非对称流场效应。
    """
    return cylinder_flow(z, U=U, a=a, alpha=alpha) + vortex(z, gamma=gamma)
