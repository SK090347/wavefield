import { describe, expect, it } from "vitest";
import { WaveSimulator, cflLimit2d, cflOk, courant } from "./solver";

describe("WaveSimulator", () => {
  it("respects 2D CFL limit helper", () => {
    expect(cflLimit2d()).toBeCloseTo(1 / Math.SQRT2, 12);
  });

  it("stays finite under safe CFL", () => {
    const sim = new WaveSimulator({
      nx: 48,
      ny: 48,
      cfl: 0.5,
      boundary: "dirichlet",
    });
    expect(cflOk(sim.cfg)).toBe(true);
    sim.addGaussian(24, 24, 1, 3);
    for (let i = 0; i < 200; i++) sim.step();
    expect(Number.isFinite(sim.maxAbs())).toBe(true);
    expect(sim.maxAbs()).toBeLessThan(50);
  });

  it("auto dt tracks courant", () => {
    const sim = new WaveSimulator({ nx: 32, ny: 32, cfl: 0.4 });
    expect(courant(sim.cfg)).toBeCloseTo(0.4 / Math.SQRT2, 8);
  });
});
