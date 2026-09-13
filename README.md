# wavefield

2D FDTD scalar wave simulator — NumPy reference, C++17 kernel, and a live Canvas heatmap. Same stencil and CFL rule in all three places; click the canvas to drop Gaussian pulses.

[![CI](https://github.com/SK090347/wavefield/actions/workflows/ci.yml/badge.svg)](https://github.com/SK090347/wavefield/actions/workflows/ci.yml)
[![License: MIT OR Apache-2.0](https://img.shields.io/badge/license-MIT%20OR%20Apache--2.0-blue.svg)](LICENSE)

Sumit Kumar Ta · Kolkata

## How it works

\[
\frac{\partial^2 u}{\partial t^2} = c^2 \nabla^2 u
\]

Leapfrog FDTD on a uniform grid with Courant numbers \(r_x = c\Delta t/\Delta x\), \(r_y = c\Delta t/\Delta y\). Isotropic CFL:

\[
\frac{c\,\Delta t}{\Delta x} \le \frac{1}{\sqrt{2}}
\]

Defaults aim for Courant ≈ 0.5. First-order Mur ABC drains outgoing waves (weak reflections remain). Discrete energy proxy + \(\max\|u\|\) catch blow-ups when you violate CFL on purpose.

Extended notes: [docs/MATH.md](docs/MATH.md).

## Quick start

**Python**

```bash
cd python
python -m pip install -r requirements.txt
python -m pytest -q
python -m wavefield --nx 128 --ny 128 --steps 300 --boundary mur
```

**C++**

```bash
cmake -S cpp -B cpp/build -DCMAKE_BUILD_TYPE=Release && cmake --build cpp/build -j
./cpp/build/wavefield_cpp --nx 128 --ny 128 --steps 300 --boundary mur
```

**Browser**

```bash
cd web && npm install && npm run dev   # http://localhost:5173
```

Toggle Mur vs Dirichlet, scrub CFL, switch colormaps.

## Layout

```
python/   NumPy solver + pytest (CFL / energy)
cpp/      C++17 kernel + CLI frame dump
web/      Vite + TypeScript heatmap
docs/     preview + MATH.md
```

## License

**MIT** OR **Apache-2.0**.
