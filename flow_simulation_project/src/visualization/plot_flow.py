"""Flow field plotting utilities."""

from __future__ import annotations

import matplotlib.pyplot as plt
import numpy as np


def plot_stream_and_potential(
    X: np.ndarray,
    Y: np.ndarray,
    phi: np.ndarray,
    psi: np.ndarray,
    a: float = 1.0,
    stagnation: np.ndarray | None = None,
    cp: np.ndarray | None = None,
    u: np.ndarray | None = None,
    v: np.ndarray | None = None,
    show_potential: bool = True,
    show_stagnation: bool = True,
    show_pressure: bool = False,
    show_quiver: bool = False,
    levels: int = 40,
) -> tuple[plt.Figure, plt.Axes]:
    """Plot streamlines, equipotential lines and optional overlays."""
    fig, ax = plt.subplots(figsize=(8, 6))
    ax.set_facecolor("#0b2340")

    if show_pressure and cp is not None:
        cp_plot = np.where(np.isfinite(cp), cp, np.nan)
        ax.contourf(X, Y, cp_plot, levels=50, cmap="turbo", alpha=0.75)

    psi_plot = np.where(np.isfinite(psi), psi, np.nan)
    ax.contour(X, Y, psi_plot, levels=levels, colors="white", linewidths=1.0, alpha=0.8)

    if show_potential:
        phi_plot = np.where(np.isfinite(phi), phi, np.nan)
        ax.contour(
            X,
            Y,
            phi_plot,
            levels=levels,
            colors="yellow",
            linewidths=0.9,
            linestyles="--",
            alpha=0.6,
        )

    if show_quiver and u is not None and v is not None:
        step = max(1, X.shape[0] // 24)
        ax.quiver(
            X[::step, ::step],
            Y[::step, ::step],
            u[::step, ::step],
            v[::step, ::step],
            color="cyan",
            alpha=0.65,
            scale=45,
        )

    cyl = plt.Circle((0.0, 0.0), a, facecolor="red", edgecolor="black", alpha=0.4, linewidth=1.2)
    ax.add_patch(cyl)

    if show_stagnation and stagnation is not None and stagnation.size > 0:
        ax.plot(np.real(stagnation), np.imag(stagnation), "b*", markersize=12)

    ax.set_aspect("equal", adjustable="box")
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_title("Cylinder Flow: Streamlines and Equipotential Lines")
    return fig, ax
