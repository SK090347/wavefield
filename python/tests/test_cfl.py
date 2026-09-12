"""CFL stability tests for the 2D FDTD scheme."""

import numpy as np
import pytest

from wavefield.solver import SimulationConfig, WaveSimulator, cfl_limit_2d


def test_cfl_limit_constant():
    assert abs(cfl_limit_2d() - 1.0 / np.sqrt(2.0)) < 1e-15


def test_stable_cfl_stays_bounded():
    cfg = SimulationConfig(nx=64, ny=64, cfl=0.5, boundary="dirichlet")
    assert cfg.cfl_ok()
    sim = WaveSimulator(cfg)
    sim.add_gaussian_pulse(32, 32, amplitude=1.0, sigma=3.0)
    peak0 = sim.max_abs()
    sim.run(400)
    # Stable scheme: amplitude should not explode.
    assert sim.max_abs() < 50.0 * peak0
    assert np.isfinite(sim.u).all()


def test_unstable_cfl_blows_up():
    # Force dt well above CFL limit.
    h = 1.0
    c = 1.0
    bad_dt = 1.1 * h / (c * np.sqrt(2.0))  # ~10% over limit
    cfg = SimulationConfig(nx=48, ny=48, dx=h, dy=h, c=c, dt=bad_dt, boundary="dirichlet")
    assert not cfg.cfl_ok()
    sim = WaveSimulator(cfg)
    sim.add_gaussian_pulse(24, 24, amplitude=1.0, sigma=2.5)
    sim.run(120)
    # Expect numerical blow-up or NaNs.
    assert (not np.isfinite(sim.u).all()) or sim.max_abs() > 1e6


def test_auto_dt_respects_cfl():
    cfg = SimulationConfig(nx=32, ny=32, cfl=0.7)
    assert cfg.cfl_ok()
    assert cfg.courrant <= 1.0 / np.sqrt(2.0) + 1e-12
