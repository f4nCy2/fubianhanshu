"""Cauchy-Riemann equation verification tools."""

from __future__ import annotations

import numpy as np
import sympy as sp


def verify_cr(
    phi: np.ndarray,
    psi: np.ndarray,
    dx: float,
    dy: float,
    mask: np.ndarray | None = None,
) -> dict[str, float]:
    """Return C-R residual statistics on a grid, optionally restricted by a mask."""
    dphi_dy, dphi_dx = np.gradient(phi, dy, dx)
    dpsi_dy, dpsi_dx = np.gradient(psi, dy, dx)

    r1 = dphi_dx - dpsi_dy
    r2 = dphi_dy + dpsi_dx

    if mask is not None:
        if mask.shape != r1.shape:
            raise ValueError("mask shape must match phi/psi shape.")
        r1 = np.where(mask, r1, np.nan)
        r2 = np.where(mask, r2, np.nan)

    return {
        "max_abs_r1": float(np.nanmax(np.abs(r1))),
        "max_abs_r2": float(np.nanmax(np.abs(r2))),
        "mean_abs_r1": float(np.nanmean(np.abs(r1))),
        "mean_abs_r2": float(np.nanmean(np.abs(r2))),
    }


def cylinder_phi_psi_polar(U: float | sp.Symbol, a: float | sp.Symbol) -> tuple[sp.Expr, sp.Expr, sp.Symbol, sp.Symbol]:
    """Return symbolic phi(r,theta), psi(r,theta) for W=U(z+a^2/z)."""
    r, theta = sp.symbols("r theta", positive=True, real=True)
    phi = U * (r + (a**2) / r) * sp.cos(theta)
    psi = U * (r - (a**2) / r) * sp.sin(theta)
    return phi, psi, r, theta


def symbolic_cr_residuals(U: float = 1.0, a: float = 1.0) -> dict[str, sp.Expr]:
    """Return simplified polar C-R residuals, which should be identically zero for r>a."""
    phi, psi, r, theta = cylinder_phi_psi_polar(U=sp.Float(U), a=sp.Float(a))

    dphi_dr = sp.diff(phi, r)
    dphi_dtheta = sp.diff(phi, theta)
    dpsi_dr = sp.diff(psi, r)
    dpsi_dtheta = sp.diff(psi, theta)

    # Polar C-R:
    # dphi/dr = (1/r) dpsi/dtheta
    # dpsi/dr = -(1/r) dphi/dtheta
    res1 = sp.simplify(dphi_dr - (1 / r) * dpsi_dtheta)
    res2 = sp.simplify(dpsi_dr + (1 / r) * dphi_dtheta)

    return {
        "residual_1": res1,
        "residual_2": res2,
        "dphi_dr": dphi_dr,
        "dphi_dtheta": dphi_dtheta,
        "dpsi_dr": dpsi_dr,
        "dpsi_dtheta": dpsi_dtheta,
    }


def random_point_cr_check(
    U: float = 1.0,
    a: float = 1.0,
    n_samples: int = 200,
    r_min_factor: float = 1.05,
    r_max_factor: float = 4.0,
    seed: int = 42,
) -> dict[str, float]:
    """Evaluate polar C-R residuals on random points in r>a using SymPy lambdify."""
    if a <= 0:
        raise ValueError("a must be positive.")
    if n_samples <= 0:
        raise ValueError("n_samples must be positive.")
    if r_min_factor <= 1.0:
        raise ValueError("r_min_factor must be greater than 1.0 to stay outside the cylinder.")
    if r_max_factor <= r_min_factor:
        raise ValueError("r_max_factor must be greater than r_min_factor.")

    r, theta = sp.symbols("r theta", positive=True, real=True)
    residuals = symbolic_cr_residuals(U=U, a=a)
    f_r1 = sp.lambdify((r, theta), residuals["residual_1"], "numpy")
    f_r2 = sp.lambdify((r, theta), residuals["residual_2"], "numpy")

    rng = np.random.default_rng(seed)
    rs = rng.uniform(r_min_factor * a, r_max_factor * a, n_samples)
    thetas = rng.uniform(0.0, 2.0 * np.pi, n_samples)

    r1_vals = np.asarray(f_r1(rs, thetas), dtype=float)
    r2_vals = np.asarray(f_r2(rs, thetas), dtype=float)

    return {
        "n_samples": float(n_samples),
        "max_abs_r1": float(np.max(np.abs(r1_vals))),
        "max_abs_r2": float(np.max(np.abs(r2_vals))),
        "mean_abs_r1": float(np.mean(np.abs(r1_vals))),
        "mean_abs_r2": float(np.mean(np.abs(r2_vals))),
    }
