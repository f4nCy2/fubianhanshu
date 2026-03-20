"""Pressure contour visualization."""

from __future__ import annotations

import matplotlib.pyplot as plt
import numpy as np


def plot_pressure_contour(
    X: np.ndarray,
    Y: np.ndarray,
    cp: np.ndarray,
    levels: int = 50,
) -> tuple[plt.Figure, plt.Axes]:
    """Plot pressure coefficient contour map."""
    fig, ax = plt.subplots(figsize=(8, 6))
    c = ax.contourf(X, Y, cp, levels=levels, cmap="coolwarm")
    fig.colorbar(c, ax=ax, label="Cp")
    ax.set_aspect("equal", adjustable="box")
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_title("Pressure Coefficient Contour")
    return fig, ax


def plot_cp_polar(theta: np.ndarray, cp_num: np.ndarray, cp_theo: np.ndarray) -> tuple[plt.Figure, plt.Axes]:
    """Plot numerical and theoretical Cp on a polar chart."""
    fig, ax = plt.subplots(subplot_kw={"projection": "polar"}, figsize=(7, 7))
    ax.plot(theta, cp_num, color="tab:blue", linewidth=2.0, label="Cp numerical")
    ax.plot(theta, cp_theo, color="tab:orange", linewidth=1.8, linestyle="--", label="Cp theory")

    key_theta = np.deg2rad([0, 90, 180, 270])
    key_cp = np.array([1.0, -3.0, 1.0, -3.0])
    ax.scatter(key_theta, key_cp, c="red", s=28, zorder=5)

    ax.set_title("Cylinder Surface Pressure Coefficient")
    ax.grid(True, alpha=0.4)
    ax.legend(loc="upper right")
    return fig, ax
