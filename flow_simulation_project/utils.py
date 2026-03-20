"""Unified helper module for coursework submission.

This module provides:
1) C-R equation verification wrappers
2) Bernoulli equation calculator class
3) Flow-field generator class
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from src.core.pressure import pressure_coefficient
from src.core.potential import total_potential
from src.core.velocity import complex_velocity as core_complex_velocity
from src.core.velocity import velocity_field_numba
from src.ui.simulator import FlowSimulator
from src.utils.cr_verify import verify_cr
from src.utils.grid import create_grid


def complex_potential_function(
    z: np.ndarray,
    U: float = 1.0,
    a: float = 1.0,
    alpha: float = 0.0,
    gamma: float = 0.0,
) -> np.ndarray:
    """Return complex potential W(z) for cylinder flow with optional circulation."""
    return total_potential(z=z, U=U, a=a, alpha=alpha, gamma=gamma)


def complex_velocity_function(
    z: np.ndarray,
    U: float = 1.0,
    a: float = 1.0,
    alpha: float = 0.0,
    gamma: float = 0.0,
) -> np.ndarray:
    """Return complex velocity dW/dz = u - iv."""
    return core_complex_velocity(z=z, U=U, a=a, alpha=alpha, gamma=gamma)


def verify_cr_equations(
    phi: np.ndarray,
    psi: np.ndarray,
    dx: float,
    dy: float,
    mask: np.ndarray | None = None,
) -> dict[str, float]:
    """Validate C-R residuals on a grid.

    Args:
        phi: Velocity potential field.
        psi: Stream function field.
        dx: Grid spacing along x.
        dy: Grid spacing along y.
        mask: Optional boolean mask for valid region.

    Returns:
        A dictionary containing max/mean absolute residuals.
    """
    return verify_cr(phi=phi, psi=psi, dx=dx, dy=dy, mask=mask)


def generate_flow_field(
    U: float = 1.0,
    a: float = 1.0,
    gamma: float = 0.0,
    alpha: float = 0.0,
    n: int = 201,
    x_min: float = -3.0,
    x_max: float = 3.0,
    y_min: float = -3.0,
    y_max: float = 3.0,
) -> dict[str, np.ndarray]:
    """Generate full flow-field outputs with a functional interface."""
    sim = FlowSimulator(
        U=U,
        a=a,
        gamma=gamma,
        alpha=alpha,
        n=n,
        x_min=x_min,
        x_max=x_max,
        y_min=y_min,
        y_max=y_max,
    )
    return sim.run()


def bernoulli_pressure_coefficient(speed: np.ndarray, U_inf: float = 1.0) -> np.ndarray:
    """Return Bernoulli pressure coefficient Cp = 1 - (V/U_inf)^2."""
    return pressure_coefficient(speed, U_inf=U_inf)


@dataclass
class BernoulliCalculator:
    """Compute Bernoulli-derived pressure metrics for incompressible steady flow.

    Args:
        rho: Fluid density in kg/m^3.
        p_inf: Free-stream static pressure in Pa.
        u_inf: Free-stream velocity magnitude in m/s.
    """

    rho: float = 1.225
    p_inf: float = 101325.0
    u_inf: float = 1.0

    def dynamic_pressure(self) -> float:
        """Return free-stream dynamic pressure q_inf = 0.5*rho*u_inf^2."""
        return 0.5 * self.rho * self.u_inf**2

    def pressure_from_speed(self, speed: np.ndarray) -> np.ndarray:
        """Return static pressure field from local speed using Bernoulli relation."""
        return self.p_inf + 0.5 * self.rho * (self.u_inf**2 - speed**2)

    def pressure_coefficient(self, speed: np.ndarray) -> np.ndarray:
        """Return pressure coefficient Cp = 1-(V/U_inf)^2."""
        return pressure_coefficient(speed, U_inf=self.u_inf)


@dataclass
class FlowFieldGenerator:
    """Generate potential-flow field data using vectorized and accelerated kernels.

    Args:
        U: Free-stream speed in m/s.
        a: Cylinder radius in m.
        gamma: Circulation in m^2/s.
        alpha: Flow incidence angle in rad.
        n: Grid size (n x n).
        x_min: Left boundary of x-domain.
        x_max: Right boundary of x-domain.
        y_min: Bottom boundary of y-domain.
        y_max: Top boundary of y-domain.
    """

    U: float = 1.0
    a: float = 1.0
    gamma: float = 0.0
    alpha: float = 0.0
    n: int = 201
    x_min: float = -3.0
    x_max: float = 3.0
    y_min: float = -3.0
    y_max: float = 3.0

    def generate(self) -> dict[str, np.ndarray]:
        """Generate full flow-field outputs with defaults matching the project spec."""
        sim = FlowSimulator(
            U=self.U,
            a=self.a,
            gamma=self.gamma,
            alpha=self.alpha,
            n=self.n,
            x_min=self.x_min,
            x_max=self.x_max,
            y_min=self.y_min,
            y_max=self.y_max,
        )
        return sim.run()

    def generate_velocity_only(self) -> dict[str, np.ndarray]:
        """Generate only velocity-related outputs for lightweight workflows."""
        X, Y, _ = create_grid(self.x_min, self.x_max, self.y_min, self.y_max, self.n)
        u, v, vmag, mask = velocity_field_numba(
            X,
            Y,
            U=self.U,
            a=self.a,
            alpha=self.alpha,
            gamma=self.gamma,
        )
        return {
            "X": X,
            "Y": Y,
            "u": u,
            "v": v,
            "V": vmag,
            "valid_mask": mask,
        }
