"""CLI: run a short simulation and optionally dump float32 frames."""

from __future__ import annotations

import argparse
import struct
from pathlib import Path

import numpy as np

from .solver import SimulationConfig, WaveSimulator


def main() -> None:
    p = argparse.ArgumentParser(description="2D FDTD wavefield (Python reference)")
    p.add_argument("--nx", type=int, default=128)
    p.add_argument("--ny", type=int, default=128)
    p.add_argument("--steps", type=int, default=200)
    p.add_argument("--cfl", type=float, default=0.5)
    p.add_argument("--boundary", choices=("mur", "dirichlet"), default="mur")
    p.add_argument("--dump-dir", type=Path, default=None)
    p.add_argument("--stride", type=int, default=4)
    args = p.parse_args()

    cfg = SimulationConfig(
        nx=args.nx, ny=args.ny, cfl=args.cfl, boundary=args.boundary
    )
    sim = WaveSimulator(cfg)
    sim.add_gaussian_pulse(args.nx * 0.5, args.ny * 0.5, amplitude=1.0, sigma=4.0)

    if args.dump_dir:
        args.dump_dir.mkdir(parents=True, exist_ok=True)
        # Header: nx, ny as uint32 little-endian
        meta = args.dump_dir / "meta.bin"
        with meta.open("wb") as f:
            f.write(struct.pack("<II", args.nx, args.ny))

    for k in range(args.steps):
        sim.step()
        if args.dump_dir and k % args.stride == 0:
            path = args.dump_dir / f"frame_{k:05d}.f32"
            sim.u.astype(np.float32).tofile(path)

    print(
        f"done steps={sim.step_index} max|u|={sim.max_abs():.6g} "
        f"energy={sim.energy():.6g} cfl={cfg.courrant:.4f} ok={cfg.cfl_ok()}"
    )


if __name__ == "__main__":
    main()
