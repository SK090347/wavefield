#!/usr/bin/env python3
"""Invoke the C++ wavefield_cpp binary via subprocess (no pybind required)."""

from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DEFAULT_BIN = ROOT / "cpp" / "build" / "wavefield_cpp"


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--bin", type=Path, default=DEFAULT_BIN)
    p.add_argument("args", nargs=argparse.REMAINDER)
    ns = p.parse_args()
    if not ns.bin.exists():
        print(
            f"missing binary: {ns.bin} — build with `make build-cpp`",
            file=sys.stderr,
        )
        return 1
    # REMAINDER may start with "--"
    args = ns.args
    if args and args[0] == "--":
        args = args[1:]
    cmd = [str(ns.bin), *args]
    print("+", " ".join(cmd))
    return subprocess.call(cmd)


if __name__ == "__main__":
    raise SystemExit(main())
