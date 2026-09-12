/** Compact colormaps for bipolar / positive wave amplitude display. */

export type CmapName = "plasma" | "seismic" | "viridis";

function lerp(a: number, b: number, t: number): number {
  return a + (b - a) * t;
}

function sampleStops(
  stops: Array<[number, number, number]>,
  t: number,
): [number, number, number] {
  const x = Math.min(1, Math.max(0, t));
  const n = stops.length - 1;
  const f = x * n;
  const i = Math.min(n - 1, Math.floor(f));
  const u = f - i;
  const a = stops[i]!;
  const b = stops[i + 1]!;
  return [lerp(a[0], b[0], u), lerp(a[1], b[1], u), lerp(a[2], b[2], u)];
}

const PLASMA: Array<[number, number, number]> = [
  [13, 8, 135],
  [84, 2, 163],
  [139, 10, 165],
  [185, 50, 137],
  [219, 92, 104],
  [244, 136, 73],
  [254, 188, 43],
  [240, 249, 33],
];

const VIRIDIS: Array<[number, number, number]> = [
  [68, 1, 84],
  [59, 82, 139],
  [33, 145, 140],
  [94, 201, 98],
  [253, 231, 37],
];

/** Map signed amplitude in [-scale, scale] → RGB. */
export function colorize(
  value: number,
  scale: number,
  cmap: CmapName,
): [number, number, number] {
  const s = Math.max(scale, 1e-12);
  if (cmap === "seismic") {
    const t = Math.min(1, Math.max(-1, value / s));
    if (t < 0) {
      // blue → white
      const u = t + 1;
      return [
        Math.round(lerp(30, 245, u)),
        Math.round(lerp(70, 245, u)),
        Math.round(lerp(200, 245, u)),
      ];
    }
    // white → red
    return [
      Math.round(lerp(245, 200, t)),
      Math.round(lerp(245, 40, t)),
      Math.round(lerp(245, 40, t)),
    ];
  }
  // Map |value| for plasma/viridis so both crests & troughs light up.
  const t = Math.min(1, Math.abs(value) / s);
  const stops = cmap === "viridis" ? VIRIDIS : PLASMA;
  const [r, g, b] = sampleStops(stops, t);
  return [Math.round(r), Math.round(g), Math.round(b)];
}
