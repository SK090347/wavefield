import { WaveSimulator, type Boundary, cflOk, courant } from "./solver";
import { colorize, type CmapName } from "./colormap";

const canvas = document.getElementById("field") as HTMLCanvasElement;
const ctx = canvas.getContext("2d", { alpha: false })!;
const speedEl = document.getElementById("speed") as HTMLInputElement;
const cflEl = document.getElementById("cfl") as HTMLInputElement;
const cflVal = document.getElementById("cflVal")!;
const boundaryEl = document.getElementById("boundary") as HTMLSelectElement;
const cmapEl = document.getElementById("cmap") as HTMLSelectElement;
const pulseBtn = document.getElementById("pulse") as HTMLButtonElement;
const resetBtn = document.getElementById("reset") as HTMLButtonElement;
const pauseBtn = document.getElementById("pause") as HTMLButtonElement;
const statStep = document.getElementById("statStep")!;
const statMax = document.getElementById("statMax")!;
const statEnergy = document.getElementById("statEnergy")!;
const statFps = document.getElementById("statFps")!;

const NX = 160;
const NY = 160;
let sim = makeSim(0.5, "mur");
let paused = false;
let cmap: CmapName = "plasma";
let displayScale = 0.4;
let frames = 0;
let lastFpsTs = performance.now();
const image = ctx.createImageData(NX, NY);
const offscreen = document.createElement("canvas");
offscreen.width = NX;
offscreen.height = NY;
const offCtx = offscreen.getContext("2d")!;

function makeSim(cfl: number, boundary: Boundary): WaveSimulator {
  const s = new WaveSimulator({ nx: NX, ny: NY, cfl, boundary });
  s.addGaussian(NX * 0.5, NY * 0.42, 1.0, 4.5);
  s.addGaussian(NX * 0.35, NY * 0.6, 0.65, 3.2);
  return s;
}

function rebuild(): void {
  const cfl = Number(cflEl.value) / 100;
  const boundary = boundaryEl.value as Boundary;
  cflVal.textContent = cfl.toFixed(2);
  sim = makeSim(cfl, boundary);
  displayScale = 0.4;
  if (!cflOk(sim.cfg)) {
    cflVal.textContent = `${cfl.toFixed(2)} (unstable!)`;
  }
}

function render(): void {
  const u = sim.u;
  const data = image.data;
  // Slow adaptive scale so the heatmap stays lively.
  const m = Math.max(sim.maxAbs(), 1e-6);
  displayScale = displayScale * 0.97 + m * 0.03;
  for (let j = 0; j < NY; j++) {
    for (let i = 0; i < NX; i++) {
      const v = u[j * NX + i]!;
      const [r, g, b] = colorize(v, displayScale, cmap);
      const p = (j * NX + i) * 4;
      data[p] = r;
      data[p + 1] = g;
      data[p + 2] = b;
      data[p + 3] = 255;
    }
  }
  // Upscale with nearest-neighbor via an offscreen buffer.
  offCtx.putImageData(image, 0, 0);
  ctx.imageSmoothingEnabled = false;
  ctx.drawImage(offscreen, 0, 0, canvas.width, canvas.height);

  statStep.textContent = String(sim.stepIndex);
  statMax.textContent = sim.maxAbs().toFixed(3);
  statEnergy.textContent = sim.energy().toExponential(2);
}

function tick(now: number): void {
  if (!paused) {
    const steps = Number(speedEl.value);
    for (let s = 0; s < steps; s++) sim.step();
    render();
    frames += 1;
  }
  if (now - lastFpsTs >= 500) {
    const fps = (frames * 1000) / (now - lastFpsTs);
    statFps.textContent = fps.toFixed(0);
    frames = 0;
    lastFpsTs = now;
  }
  requestAnimationFrame(tick);
}

canvas.addEventListener("pointerdown", (ev) => {
  const rect = canvas.getBoundingClientRect();
  const x = ((ev.clientX - rect.left) / rect.width) * NX;
  const y = ((ev.clientY - rect.top) / rect.height) * NY;
  sim.addGaussian(x, y, 0.9, 3.5);
});

pulseBtn.addEventListener("click", () => {
  const x = NX * (0.25 + Math.random() * 0.5);
  const y = NY * (0.25 + Math.random() * 0.5);
  sim.addGaussian(x, y, 0.85 + Math.random() * 0.4, 3 + Math.random() * 2);
});

resetBtn.addEventListener("click", () => rebuild());
pauseBtn.addEventListener("click", () => {
  paused = !paused;
  pauseBtn.textContent = paused ? "Resume" : "Pause";
});
cflEl.addEventListener("input", () => rebuild());
boundaryEl.addEventListener("change", () => rebuild());
cmapEl.addEventListener("change", () => {
  cmap = cmapEl.value as CmapName;
});

console.info(
  `wavefield ready · courant=${courant(sim.cfg).toFixed(3)} cfl_ok=${cflOk(sim.cfg)}`,
);
render();
requestAnimationFrame(tick);
