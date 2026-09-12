"""Energy-ish sanity checks (Dirichlet box conserves discrete energy roughly)."""

import numpy as np

from wavefield.solver import SimulationConfig, WaveSimulator


def test_energy_does_not_grow_unbounded_dirichlet():
    cfg = SimulationConfig(nx=64, ny=64, cfl=0.45, boundary="dirichlet", damping=0.0)
    sim = WaveSimulator(cfg)
    sim.add_gaussian_pulse(32, 32, amplitude=1.0, sigma=4.0)
    # Warm up a few steps so kinetic term is nonzero.
    sim.run(5)
    e0 = sim.energy()
    energies = [e0]
    for _ in range(200):
        sim.step()
        energies.append(sim.energy())
    # Relative growth should stay modest for a conservative interior scheme.
    e_max = max(energies)
    assert e_max < 5.0 * e0 + 1e-9
    assert np.isfinite(energies).all()


def test_mur_absorbs_outgoing_energy():
    """Pulse near center; after long time Mur ABC should leave little energy."""
    cfg = SimulationConfig(nx=80, ny=80, cfl=0.5, boundary="mur")
    sim = WaveSimulator(cfg)
    sim.add_gaussian_pulse(40, 40, amplitude=1.0, sigma=3.0)
    sim.run(5)
    e_early = sim.energy()
    sim.run(600)  # wave exits the domain
    e_late = sim.energy()
    assert e_late < 0.15 * e_early


def test_field_remains_finite():
    cfg = SimulationConfig(nx=40, ny=40, cfl=0.5, boundary="mur")
    sim = WaveSimulator(cfg)
    sim.add_gaussian_pulse(20, 20, amplitude=2.0, sigma=2.0)
    sim.run(300)
    assert np.isfinite(sim.u).all()
    assert sim.max_abs() < 100.0
