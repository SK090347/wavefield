"""Unit tests for WaveSimulator API and symmetry."""

import numpy as np
import pytest

from wavefield.solver import SimulationConfig, WaveSimulator


def test_gaussian_is_symmetric():
    cfg = SimulationConfig(nx=65, ny=65, cfl=0.5, boundary="dirichlet")
    sim = WaveSimulator(cfg)
    sim.add_gaussian_pulse(32, 32, amplitude=1.0, sigma=4.0)
    u = sim.u
    assert np.allclose(u, u.T, atol=1e-12)
    assert np.allclose(u, np.flipud(u), atol=1e-12)
    assert np.allclose(u, np.fliplr(u), atol=1e-12)


def test_reset_clears_state():
    sim = WaveSimulator(SimulationConfig(nx=16, ny=16))
    sim.add_gaussian_pulse(8, 8)
    sim.run(10)
    sim.reset()
    assert sim.step_index == 0
    assert sim.max_abs() == 0.0


def test_point_source_index_error():
    sim = WaveSimulator(SimulationConfig(nx=10, ny=10))
    with pytest.raises(IndexError):
        sim.add_point_source(100, 0, 1.0)


def test_tiny_grid_rejected():
    with pytest.raises(ValueError):
        WaveSimulator(SimulationConfig(nx=2, ny=2))


def test_frames_iterator_length():
    sim = WaveSimulator(SimulationConfig(nx=24, ny=24, cfl=0.5))
    sim.add_gaussian_pulse(12, 12, sigma=2.0)
    frames = list(sim.frames(5, stride=2))
    assert len(frames) == 5
    assert frames[-1][0] == 10
