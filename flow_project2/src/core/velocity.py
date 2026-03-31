"""Velocity and stagnation-point tools for Project 2."""

from __future__ import annotations

import numpy as np

try:
    from numba import njit
except Exception:  # pragma: no cover
    def njit(*args, **kwargs):
        def _decorator(func):
            return func

        return _decorator


@njit(cache=True)
def _velocity_components_kernel(
    X: np.ndarray,
    Y: np.ndarray,
    U: float,
    a: float,
    gamma: float,
    eps: float,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Numba kernel for velocity components outside cylinder."""
    ny, nx = X.shape
    u = np.empty_like(X)
    v = np.empty_like(Y)
    mask = np.empty_like(X, dtype=np.bool_)

    for i in range(ny):
        for j in range(nx):
            x = X[i, j]
            y = Y[i, j]
            r2 = x * x + y * y

            if r2 <= a * a:
                mask[i, j] = False
                u[i, j] = np.nan
                v[i, j] = np.nan
                continue

            z = complex(x, y)
            z_safe = z if abs(z) > eps else complex(eps, 0.0)
            dwdz = U * (1.0 - (a * a) / (z_safe * z_safe)) + 1j * gamma / (2.0 * np.pi * z_safe)

            mask[i, j] = True
            u[i, j] = dwdz.real
            v[i, j] = -dwdz.imag

    return u, v, mask


def complex_velocity(z: np.ndarray, U: float = 1.0, a: float = 1.0, gamma: float = 0.0) -> np.ndarray:
    """功能：计算复速度 dW/dz。

    输入：
        z: 复平面坐标数组。
        U: 来流速度，单位 m/s。
        a: 圆柱半径，单位 m。
        gamma: 环量，单位 m^2/s。

    输出：
        复速度 dW/dz = u - iv。

    物理参数说明：
        复速度对应速度势和流函数导数，反映局部速度分布。
    """
    if a <= 0:
        raise ValueError("a must be positive.")

    z_safe = np.where(np.abs(z) < 1e-12, 1e-12 + 0j, z)
    return U * (1.0 - (a**2) / (z_safe**2)) + 1j * gamma / (2.0 * np.pi * z_safe)


def velocity_field_numba(
    X: np.ndarray,
    Y: np.ndarray,
    U: float = 1.0,
    a: float = 1.0,
    gamma: float = 0.0,
    eps: float = 1e-12,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """功能：在网格上计算速度场并屏蔽圆柱内部。

    输入：
        X, Y: 网格坐标。
        U: 来流速度，单位 m/s。
        a: 圆柱半径，单位 m。
        gamma: 环量，单位 m^2/s。
        eps: 数值奇异保护阈值。

    输出：
        u, v, v_mag, valid_mask。

    物理参数说明：
        环量项改变切向速度分布，导致压力不对称。
    """
    u, v, valid_mask = _velocity_components_kernel(X, Y, U, a, gamma, eps)
    v_mag = np.sqrt(u**2 + v**2)
    return u, v, v_mag, valid_mask


def surface_tangential_velocity(theta: np.ndarray, U: float, a: float, gamma: float) -> np.ndarray:
    """功能：计算圆柱表面切向速度。

    输入：
        theta: 圆周角数组，单位 rad。
        U: 来流速度，单位 m/s。
        a: 圆柱半径，单位 m。
        gamma: 环量，单位 m^2/s。

    输出：
        切向速度 v_theta(theta) = -2U sin(theta) - Gamma/(2*pi*a)。

    物理参数说明：
        该表达式用于驻点条件与压力分布分析。
    """
    if a <= 0:
        raise ValueError("a must be positive.")
    return -2.0 * U * np.sin(theta) - gamma / (2.0 * np.pi * a)


def stagnation_angles(U: float, a: float, gamma: float) -> np.ndarray:
    """功能：求解圆柱表面驻点角度。

    输入：
        U: 来流速度，单位 m/s。
        a: 圆柱半径，单位 m。
        gamma: 环量，单位 m^2/s。

    输出：
        驻点角度数组（rad）。若不存在表面驻点则返回空数组。

    物理参数说明：
        驻点条件 sin(theta_s) = -Gamma/(4*pi*U*a)，
        当 |Gamma|>4*pi*U*a 时驻点脱离圆柱表面。
    """
    if U == 0:
        raise ValueError("U must be non-zero.")
    if a <= 0:
        raise ValueError("a must be positive.")

    rhs = -gamma / (4.0 * np.pi * U * a)
    if abs(rhs) > 1.0:
        return np.array([], dtype=float)

    theta_1 = np.arcsin(rhs)
    theta_2 = np.pi - theta_1
    return np.array([theta_1, theta_2], dtype=float)


def stagnation_points(U: float, a: float, gamma: float) -> np.ndarray:
    """功能：返回表面驻点复坐标。

    输入：
        U: 来流速度，单位 m/s。
        a: 圆柱半径，单位 m。
        gamma: 环量，单位 m^2/s。

    输出：
        驻点复坐标数组；若驻点离开表面则返回空数组。

    物理参数说明：
        用于可视化驻点迁移与临界环量现象。
    """
    angles = stagnation_angles(U=U, a=a, gamma=gamma)
    if angles.size == 0:
        return np.array([], dtype=complex)
    return a * np.exp(1j * angles)
