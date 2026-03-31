"""Flow-field visualization for Project 2."""

from __future__ import annotations

import matplotlib.pyplot as plt
import numpy as np


def plot_stream_and_potential(
    X: np.ndarray,
    Y: np.ndarray,
    phi: np.ndarray,
    psi: np.ndarray,
    a: float,
    cp: np.ndarray | None = None,
    stagnation: np.ndarray | None = None,
    show_pressure: bool = False,
    show_potential: bool = True,
    levels: int = 45,
) -> tuple[plt.Figure, plt.Axes]:
    """Plot streamlines/equipotential lines with optional pressure background."""
    fig, ax = plt.subplots(figsize=(9, 7))
    ax.set_facecolor("#0c2333")

    if show_pressure and cp is not None:
        cp_plot = np.where(np.isfinite(cp), cp, np.nan)
        ax.contourf(X, Y, cp_plot, levels=55, cmap="turbo", alpha=0.7)

    ax.contour(X, Y, psi, levels=levels, colors="white", linewidths=0.9, alpha=0.85)

    if show_potential:
        ax.contour(X, Y, phi, levels=levels, colors="#ffd166", linewidths=0.85, linestyles="--", alpha=0.7)

    cyl = plt.Circle((0.0, 0.0), a, facecolor="#f94144", edgecolor="black", linewidth=1.1, alpha=0.45)
    ax.add_patch(cyl)

    if stagnation is not None and stagnation.size > 0:
        ax.plot(np.real(stagnation), np.imag(stagnation), "c*", markersize=11)

    ax.set_aspect("equal", adjustable="box")
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_title("Project2: Circulation-Controlled Cylinder Flow")
    return fig, ax
