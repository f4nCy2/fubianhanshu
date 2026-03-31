"""Unified helper module for Project 2 submission."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from src.core.potential import total_potential
from src.core.pressure import (
    cylinder_surface_pressure,
    lift_theory_kj,
    pressure_coefficient,
    safe_gamma_range,
)
from src.core.velocity import complex_velocity, velocity_field_numba
from src.ui.simulator import FlowSimulator
from src.utils.cr_verify import verify_cr


def complex_potential_function(
    z: np.ndarray,
    U: float = 1.0,
    a: float = 1.0,
    gamma: float = 0.0,
) -> np.ndarray:
    """Return complex potential W(z) for Project 2 model."""
    return total_potential(z=z, U=U, a=a, gamma=gamma)


def complex_velocity_function(
    z: np.ndarray,
    U: float = 1.0,
    a: float = 1.0,
    gamma: float = 0.0,
) -> np.ndarray:
    """Return complex velocity dW/dz = u-iv."""
    return complex_velocity(z=z, U=U, a=a, gamma=gamma)


def verify_cr_equations(
    phi: np.ndarray,
    psi: np.ndarray,
    dx: float,
    dy: float,
    mask: np.ndarray | None = None,
) -> dict[str, float]:
    """Validate C-R residuals for potential/stream-function fields."""
    return verify_cr(phi=phi, psi=psi, dx=dx, dy=dy, mask=mask)


def generate_flow_field(
    U: float = 1.0,
    a: float = 1.0,
    gamma: float = 0.0,
    n: int = 201,
    x_min: float = -3.0,
    x_max: float = 3.0,
    y_min: float = -3.0,
    y_max: float = 3.0,
) -> dict[str, np.ndarray]:
    """Generate full flow field outputs using FlowSimulator."""
    sim = FlowSimulator(
        U=U,
        a=a,
        gamma=gamma,
        n=n,
        x_min=x_min,
        x_max=x_max,
        y_min=y_min,
        y_max=y_max,
    )
    return sim.run()


def bernoulli_pressure_coefficient(speed: np.ndarray, U_inf: float = 1.0) -> np.ndarray:
    """Return Bernoulli pressure coefficient Cp."""
    return pressure_coefficient(speed, U_inf=U_inf)


@dataclass
class BernoulliCalculator:
    """Compute pressure metrics from Bernoulli relation."""

    rho: float = 1.225
    p_inf: float = 101325.0
    u_inf: float = 1.0

    def dynamic_pressure(self) -> float:
        return 0.5 * self.rho * self.u_inf**2

    def pressure_from_speed(self, speed: np.ndarray) -> np.ndarray:
        return self.p_inf + 0.5 * self.rho * (self.u_inf**2 - speed**2)

    def pressure_coefficient(self, speed: np.ndarray) -> np.ndarray:
        return pressure_coefficient(speed, U_inf=self.u_inf)


@dataclass
class FlowFieldGenerator:
    """Generate velocity fields for Project 2."""

    U: float = 1.0
    a: float = 1.0
    gamma: float = 0.0
    n: int = 201
    x_min: float = -3.0
    x_max: float = 3.0
    y_min: float = -3.0
    y_max: float = 3.0

    def generate(self) -> dict[str, np.ndarray]:
        return generate_flow_field(
            U=self.U,
            a=self.a,
            gamma=self.gamma,
            n=self.n,
            x_min=self.x_min,
            x_max=self.x_max,
            y_min=self.y_min,
            y_max=self.y_max,
        )

    def generate_velocity_only(self) -> dict[str, np.ndarray]:
        out = self.generate()
        return {
            "X": out["X"],
            "Y": out["Y"],
            "u": out["u"],
            "v": out["v"],
            "V": out["V"],
            "valid_mask": out["valid_mask"],
        }


def project2_key_metrics(U: float = 1.0, a: float = 1.0, gamma: float = 0.0, rho: float = 1.225) -> dict[str, float]:
    """Return key engineering metrics used in Project 2 report."""
    gamma_critical, gamma_safe_pos, gamma_safe_neg = safe_gamma_range(U=U, a=a, safety_factor=1.5)
    lift_theory = lift_theory_kj(rho=rho, U=U, gamma=gamma)
    surface = cylinder_surface_pressure(U=U, a=a, gamma=gamma, rho=rho)
    cp_min = float(np.nanmin(surface["cp"]))
    cp_max = float(np.nanmax(surface["cp"]))

    return {
        "gamma_critical": gamma_critical,
        "gamma_safe_pos": gamma_safe_pos,
        "gamma_safe_neg": gamma_safe_neg,
        "lift_theory": lift_theory,
        "cp_min": cp_min,
        "cp_max": cp_max,
    }
