"""Pressure-related visualization for Project 2."""

from __future__ import annotations

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


def plot_surface_cp(theta: np.ndarray, cp: np.ndarray) -> tuple[plt.Figure, plt.Axes]:
    """Plot cylinder surface pressure coefficient against angle."""
    fig, ax = plt.subplots(figsize=(9, 5))
    ax.plot(np.degrees(theta), cp, color="#0077b6", linewidth=2.0)
    ax.axhline(0.0, color="black", linewidth=0.9, linestyle="--", alpha=0.7)
    ax.set_xlabel("theta (deg)")
    ax.set_ylabel("Cp")
    ax.set_title("Surface Pressure Coefficient Distribution")
    ax.grid(alpha=0.35)
    return fig, ax


def plot_gamma_lift_curve(df: pd.DataFrame) -> tuple[plt.Figure, plt.Axes]:
    """Plot numerical/theoretical lift versus circulation."""
    fig, ax = plt.subplots(figsize=(9, 5.5))
    ax.plot(df["gamma"], df["lift_num"], color="#e76f51", linewidth=2.2, label="L_num")
    ax.plot(df["gamma"], df["lift_theory"], color="#264653", linewidth=1.8, linestyle="--", label="L_theory")
    ax.set_xlabel("Gamma (m^2/s)")
    ax.set_ylabel("Lift per span (N/m)")
    ax.set_title("Gamma-Lift Verification (Kutta-Joukowski)")
    ax.grid(alpha=0.35)
    ax.legend()
    return fig, ax
