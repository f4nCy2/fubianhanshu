"""Velocity computation utilities."""

from __future__ import annotations

import numpy as np

try:
    from numba import njit
except Exception:  # pragma: no cover - optional fallback when numba is unavailable
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
    alpha: float,
    gamma: float,
    eps: float,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Numba kernel for masked velocity field calculation on meshgrid."""
    ny, nx = X.shape
    u = np.empty_like(X)
    v = np.empty_like(Y)
    mask = np.empty_like(X, dtype=np.bool_)

    cos_a = np.cos(alpha)
    sin_a = np.sin(alpha)

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

            # U*exp(-i*alpha)*(1-a^2/z^2)
            base = U * complex(cos_a, -sin_a) * (1.0 - (a * a) / (z_safe * z_safe))
            circulation = -1j * gamma / (2.0 * np.pi * z_safe)
            dwdz = base + circulation

            mask[i, j] = True
            u[i, j] = dwdz.real
            v[i, j] = -dwdz.imag

    return u, v, mask


def complex_velocity(
    z: np.ndarray,
    U: float = 1.0,
    a: float = 1.0,
    alpha: float = 0.0,
    gamma: float = 0.0,
) -> np.ndarray:
    """功能：计算复速度 dW/dz（圆柱绕流，可叠加环量）。

    输入：
        z: 复平面坐标数组。
        U: 远场来流速度，单位 m/s。
        a: 圆柱半径，单位 m。
        alpha: 来流角，单位 rad。
        gamma: 环量，单位 m^2/s。

    输出：
        复速度数组 dW/dz = u - i v。

    物理参数说明：
        gamma 控制环量引起的切向速度附加项。
    """
    z_safe = np.where(np.abs(z) < 1e-12, 1e-12 + 0j, z)
    base = U * np.exp(-1j * alpha) * (1 - (a**2) / (z_safe**2))
    circulation = -1j * gamma / (2.0 * np.pi * z_safe)
    return base + circulation


def velocity_field_numba(
    X: np.ndarray,
    Y: np.ndarray,
    U: float = 1.0,
    a: float = 1.0,
    alpha: float = 0.0,
    gamma: float = 0.0,
    eps: float = 1e-12,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """功能：在网格上计算速度场并屏蔽圆柱内部区域。

    输入：
        X, Y: 网格坐标（meshgrid 结果）。
        U: 远场来流速度，单位 m/s。
        a: 圆柱半径，单位 m。
        alpha: 来流角，单位 rad。
        gamma: 环量，单位 m^2/s。
        eps: 奇异点保护阈值。

    输出：
        u, v, vmag, valid_mask:
        其中 valid_mask=False 表示圆柱内部无效点。

    物理参数说明：
        a 定义障碍物边界；eps 用于避免 z=0 导致除零。
    """
    u, v, valid_mask = _velocity_components_kernel(X, Y, U, a, alpha, gamma, eps)
    vmag = np.sqrt(u**2 + v**2)
    return u, v, vmag, valid_mask


def uv_from_complex_velocity(dwdz: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """功能：将复速度转换为实部速度分量。

    输入：
        dwdz: 复速度数组，满足 dW/dz = u - i v。

    输出：
        (u, v): x 与 y 方向速度分量。

    物理参数说明：
        该变换对应势流复变量表示与笛卡尔速度分量映射。
    """
    u = np.real(dwdz)
    v = -np.imag(dwdz)
    return u, v


def speed(u: np.ndarray, v: np.ndarray) -> np.ndarray:
    """功能：计算速度模长。 

    输入：
        u, v: 速度分量数组。

    输出：
        |V| = sqrt(u^2 + v^2) 数组。

    物理参数说明：
        |V| 常用于伯努利关系和压力系数计算。
    """
    return np.sqrt(u**2 + v**2)


def stagnation_points(U: float = 1.0, a: float = 1.0, alpha: float = 0.0, gamma: float = 0.0) -> np.ndarray:
    """功能：求解驻点位置（dW/dz = 0）。

    输入：
        U: 远场来流速度，单位 m/s。
        a: 圆柱半径，单位 m。
        alpha: 来流角，单位 rad。
        gamma: 环量，单位 m^2/s。

    输出：
        复平面驻点数组（通常为两个根）。

    物理参数说明：
        gamma=0 且 alpha=0 时，驻点理论上位于 z = ±a。
    """
    if U == 0:
        raise ValueError("U must be non-zero for stagnation-point calculation.")
    if a <= 0:
        raise ValueError("a must be positive.")

    c2 = U * np.exp(-1j * alpha)
    c1 = -1j * gamma / (2.0 * np.pi)
    c0 = -U * np.exp(-1j * alpha) * (a**2)
    roots = np.roots(np.array([c2, c1, c0], dtype=complex))
    return roots
