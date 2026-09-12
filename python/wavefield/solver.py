"""
2D finite-difference time-domain (FDTD) solver for the scalar wave equation:

    ∂²u/∂t² = c² ∇²u

Central second differences in space and time give the explicit leapfrog update

    u^{n+1} = 2 u^n − u^{n−1} + r² (u_{i±1,j} + u_{i,j±1} − 4 u^n)

with Courant number r = c Δt / Δx.  In 2D the CFL stability bound is r ≤ 1/√2.
Optional first-order Mur absorbing boundary conditions reduce reflections at
the domain edges.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Iterator, Literal, Optional, Tuple

import numpy as np


Boundary = Literal["dirichlet", "mur"]


@dataclass(frozen=True)
class SimulationConfig:
    """Grid, material, and time-step parameters."""

    nx: int = 128
    ny: int = 128
    dx: float = 1.0
    dy: float = 1.0
    c: float = 1.0
    dt: Optional[float] = None  # chosen for CFL safety if None
    cfl: float = 0.5  # target Courant number when dt is auto
    boundary: Boundary = "mur"
    damping: float = 0.0  # optional interior sponge (0 = off)

    def resolved_dt(self) -> float:
        if self.dt is not None:
            return float(self.dt)
        # Use min spacing so anisotropic grids stay stable.
        h = min(self.dx, self.dy)
        return float(self.cfl * h / (self.c * np.sqrt(2.0)))

    @property
    def courrant(self) -> float:
        """Courant number r = c Δt / Δx (using dx)."""
        return self.c * self.resolved_dt() / self.dx

    def cfl_ok(self) -> bool:
        """True if the 2D CFL condition holds for isotropic dx=dy."""
        h = min(self.dx, self.dy)
        r = self.c * self.resolved_dt() / h
        return r <= 1.0 / np.sqrt(2.0) + 1e-12


class WaveSimulator:
    """In-place 2D FDTD wavefield with optional Mur ABC."""

    def __init__(self, config: Optional[SimulationConfig] = None) -> None:
        self.config = config or SimulationConfig()
        if self.config.nx < 4 or self.config.ny < 4:
            raise ValueError("grid must be at least 4×4")
        self.dt = self.config.resolved_dt()
        self.rx2 = (self.config.c * self.dt / self.config.dx) ** 2
        self.ry2 = (self.config.c * self.dt / self.config.dy) ** 2
        self.u = np.zeros((self.config.ny, self.config.nx), dtype=np.float64)
        self.u_prev = np.zeros_like(self.u)
        self.u_next = np.zeros_like(self.u)
        self.step_index = 0
        # Mur coefficients (first-order, face-normal).
        self._mur_x = (self.config.c * self.dt - self.config.dx) / (
            self.config.c * self.dt + self.config.dx
        )
        self._mur_y = (self.config.c * self.dt - self.config.dy) / (
            self.config.c * self.dt + self.config.dy
        )

    def reset(self) -> None:
        self.u.fill(0.0)
        self.u_prev.fill(0.0)
        self.u_next.fill(0.0)
        self.step_index = 0

    def add_gaussian_pulse(
        self,
        x0: float,
        y0: float,
        amplitude: float = 1.0,
        sigma: float = 3.0,
        into: str = "both",
    ) -> None:
        """Soft Gaussian initial condition centered at (x0, y0) in grid coords."""
        yy, xx = np.mgrid[0 : self.config.ny, 0 : self.config.nx]
        blob = amplitude * np.exp(-((xx - x0) ** 2 + (yy - y0) ** 2) / (2.0 * sigma**2))
        if into in ("u", "both"):
            self.u += blob
        if into in ("prev", "both"):
            # Stationary start: u_prev = u so first step has zero velocity.
            self.u_prev += blob

    def add_point_source(self, i: int, j: int, value: float) -> None:
        """Hard source at integer cell (i=x, j=y) applied to current field."""
        if not (0 <= i < self.config.nx and 0 <= j < self.config.ny):
            raise IndexError("source outside grid")
        self.u[j, i] += value

    def energy(self) -> float:
        """Discrete energy proxy: kinetic-ish + potential-ish (∑ (Δu)² + (∇u)²)."""
        # Approximate ∑ [(u−u_prev)/dt]² + c² |∇u|²  (scaled, relative use only)
        kinetic = np.sum((self.u - self.u_prev) ** 2)
        gx = np.diff(self.u, axis=1)
        gy = np.diff(self.u, axis=0)
        potential = np.sum(gx**2) + np.sum(gy**2)
        return float(kinetic + self.rx2 * potential)

    def max_abs(self) -> float:
        return float(np.max(np.abs(self.u)))

    def step(self) -> np.ndarray:
        """Advance one time step; returns view of the new field (ny, nx)."""
        cfg = self.config
        u, up, un = self.u, self.u_prev, self.u_next
        rx2, ry2 = self.rx2, self.ry2

        # Interior update (vectorized).
        lap = (
            rx2 * (u[1:-1, 2:] + u[1:-1, :-2])
            + ry2 * (u[2:, 1:-1] + u[:-2, 1:-1])
            - 2.0 * (rx2 + ry2) * u[1:-1, 1:-1]
        )
        un[1:-1, 1:-1] = 2.0 * u[1:-1, 1:-1] - up[1:-1, 1:-1] + lap

        if cfg.damping > 0.0:
            un[1:-1, 1:-1] *= 1.0 - cfg.damping

        if cfg.boundary == "dirichlet":
            un[0, :] = 0.0
            un[-1, :] = 0.0
            un[:, 0] = 0.0
            un[:, -1] = 0.0
        else:
            # First-order Mur ABC on four faces (corners left as averages).
            mx, my = self._mur_x, self._mur_y
            un[1:-1, 0] = u[1:-1, 1] + mx * (un[1:-1, 1] - u[1:-1, 0])
            un[1:-1, -1] = u[1:-1, -2] + mx * (un[1:-1, -2] - u[1:-1, -1])
            un[0, 1:-1] = u[1, 1:-1] + my * (un[1, 1:-1] - u[0, 1:-1])
            un[-1, 1:-1] = u[-2, 1:-1] + my * (un[-2, 1:-1] - u[-1, 1:-1])
            # Corners: average of adjacent ABC faces.
            un[0, 0] = 0.5 * (un[0, 1] + un[1, 0])
            un[0, -1] = 0.5 * (un[0, -2] + un[1, -1])
            un[-1, 0] = 0.5 * (un[-1, 1] + un[-2, 0])
            un[-1, -1] = 0.5 * (un[-1, -2] + un[-2, -1])

        # Rotate buffers.
        self.u_prev = u
        self.u = un
        self.u_next = up  # reuse storage
        self.u_next.fill(0.0)
        self.step_index += 1
        return self.u

    def run(self, steps: int) -> np.ndarray:
        for _ in range(steps):
            self.step()
        return self.u

    def frames(
        self, n_frames: int, stride: int = 1
    ) -> Iterator[Tuple[int, np.ndarray]]:
        """Yield (step_index, copy of field) every `stride` steps."""
        for _ in range(n_frames):
            for _ in range(stride):
                self.step()
            yield self.step_index, self.u.copy()


def cfl_limit_2d() -> float:
    """Analytic 2D CFL limit for the standard 5-point stencil: 1/√2."""
    return 1.0 / np.sqrt(2.0)
