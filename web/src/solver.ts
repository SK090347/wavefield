/** 2D FDTD scalar wave equation — same scheme as Python / C++. */

export type Boundary = "mur" | "dirichlet";

export interface SimConfig {
  nx: number;
  ny: number;
  dx: number;
  dy: number;
  c: number;
  cfl: number;
  boundary: Boundary;
  damping: number;
}

export function cflLimit2d(): number {
  return 1 / Math.SQRT2;
}

export function resolvedDt(cfg: SimConfig): number {
  const h = Math.min(cfg.dx, cfg.dy);
  return (cfg.cfl * h) / (cfg.c * Math.SQRT2);
}

export function courant(cfg: SimConfig): number {
  return (cfg.c * resolvedDt(cfg)) / cfg.dx;
}

export function cflOk(cfg: SimConfig): boolean {
  const h = Math.min(cfg.dx, cfg.dy);
  const r = (cfg.c * resolvedDt(cfg)) / h;
  return r <= cflLimit2d() + 1e-12;
}

export class WaveSimulator {
  readonly cfg: SimConfig;
  readonly dt: number;
  private rx2: number;
  private ry2: number;
  private murX: number;
  private murY: number;
  u: Float64Array;
  uPrev: Float64Array;
  uNext: Float64Array;
  stepIndex = 0;

  constructor(cfg: Partial<SimConfig> = {}) {
    this.cfg = {
      nx: 128,
      ny: 128,
      dx: 1,
      dy: 1,
      c: 1,
      cfl: 0.5,
      boundary: "mur",
      damping: 0,
      ...cfg,
    };
    if (this.cfg.nx < 4 || this.cfg.ny < 4) {
      throw new Error("grid must be at least 4×4");
    }
    this.dt = resolvedDt(this.cfg);
    const cx = (this.cfg.c * this.dt) / this.cfg.dx;
    const cy = (this.cfg.c * this.dt) / this.cfg.dy;
    this.rx2 = cx * cx;
    this.ry2 = cy * cy;
    this.murX =
      (this.cfg.c * this.dt - this.cfg.dx) /
      (this.cfg.c * this.dt + this.cfg.dx);
    this.murY =
      (this.cfg.c * this.dt - this.cfg.dy) /
      (this.cfg.c * this.dt + this.cfg.dy);
    const n = this.cfg.nx * this.cfg.ny;
    this.u = new Float64Array(n);
    this.uPrev = new Float64Array(n);
    this.uNext = new Float64Array(n);
  }

  private idx(i: number, j: number): number {
    return j * this.cfg.nx + i;
  }

  reset(): void {
    this.u.fill(0);
    this.uPrev.fill(0);
    this.uNext.fill(0);
    this.stepIndex = 0;
  }

  addGaussian(
    x0: number,
    y0: number,
    amplitude = 1,
    sigma = 3,
  ): void {
    const inv = 1 / (2 * sigma * sigma);
    const { nx, ny } = this.cfg;
    for (let j = 0; j < ny; j++) {
      for (let i = 0; i < nx; i++) {
        const dx = i - x0;
        const dy = j - y0;
        const v = amplitude * Math.exp(-(dx * dx + dy * dy) * inv);
        const k = this.idx(i, j);
        this.u[k] += v;
        this.uPrev[k] += v;
      }
    }
  }

  step(): void {
    const { nx, ny, boundary, damping } = this.cfg;
    const u = this.u;
    const up = this.uPrev;
    const un = this.uNext;
    const rx2 = this.rx2;
    const ry2 = this.ry2;

    for (let j = 1; j < ny - 1; j++) {
      for (let i = 1; i < nx - 1; i++) {
        const k = this.idx(i, j);
        const lap =
          rx2 * (u[k + 1] + u[k - 1]) +
          ry2 * (u[k + nx] + u[k - nx]) -
          2 * (rx2 + ry2) * u[k];
        let next = 2 * u[k] - up[k] + lap;
        if (damping > 0) next *= 1 - damping;
        un[k] = next;
      }
    }

    if (boundary === "dirichlet") {
      for (let i = 0; i < nx; i++) {
        un[this.idx(i, 0)] = 0;
        un[this.idx(i, ny - 1)] = 0;
      }
      for (let j = 0; j < ny; j++) {
        un[this.idx(0, j)] = 0;
        un[this.idx(nx - 1, j)] = 0;
      }
    } else {
      const mx = this.murX;
      const my = this.murY;
      for (let j = 1; j < ny - 1; j++) {
        const L = this.idx(0, j);
        const R = this.idx(nx - 1, j);
        un[L] = u[L + 1] + mx * (un[L + 1] - u[L]);
        un[R] = u[R - 1] + mx * (un[R - 1] - u[R]);
      }
      for (let i = 1; i < nx - 1; i++) {
        const B = this.idx(i, 0);
        const T = this.idx(i, ny - 1);
        un[B] = u[B + nx] + my * (un[B + nx] - u[B]);
        un[T] = u[T - nx] + my * (un[T - nx] - u[T]);
      }
      un[this.idx(0, 0)] = 0.5 * (un[this.idx(1, 0)] + un[this.idx(0, 1)]);
      un[this.idx(nx - 1, 0)] =
        0.5 * (un[this.idx(nx - 2, 0)] + un[this.idx(nx - 1, 1)]);
      un[this.idx(0, ny - 1)] =
        0.5 * (un[this.idx(1, ny - 1)] + un[this.idx(0, ny - 2)]);
      un[this.idx(nx - 1, ny - 1)] =
        0.5 * (un[this.idx(nx - 2, ny - 1)] + un[this.idx(nx - 1, ny - 2)]);
    }

    this.uPrev = u;
    this.u = un;
    this.uNext = up;
    this.uNext.fill(0);
    this.stepIndex += 1;
  }

  maxAbs(): number {
    let m = 0;
    for (let i = 0; i < this.u.length; i++) {
      const a = Math.abs(this.u[i]!);
      if (a > m) m = a;
    }
    return m;
  }

  energy(): number {
    let kinetic = 0;
    for (let i = 0; i < this.u.length; i++) {
      const d = this.u[i]! - this.uPrev[i]!;
      kinetic += d * d;
    }
    let potential = 0;
    const { nx, ny } = this.cfg;
    for (let j = 0; j < ny; j++) {
      for (let i = 0; i < nx - 1; i++) {
        const g = this.u[this.idx(i + 1, j)]! - this.u[this.idx(i, j)]!;
        potential += g * g;
      }
    }
    for (let j = 0; j < ny - 1; j++) {
      for (let i = 0; i < nx; i++) {
        const g = this.u[this.idx(i, j + 1)]! - this.u[this.idx(i, j)]!;
        potential += g * g;
      }
    }
    return kinetic + this.rx2 * potential;
  }
}
