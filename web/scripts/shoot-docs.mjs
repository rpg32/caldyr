// Docs screenshot tool — regenerates the tutorial images under docs/img/.
//
// Drives the running web app with headless Playwright (the Chrome/CDP path is
// known-flaky on the dev machine — see docs/UX_BACKLOG.md). Requires BOTH
// servers up: engine API on :8753 and the web app on :5273.
//
//   cd web && npm run dev                       # (in another terminal, + the API)
//   node scripts/shoot-docs.mjs
//
// Captures five UI shots. The sixth docs figure (distillation-reflux-tradeoff.png)
// is a data plot generated from the engine by scripts/reflux_chart.py, not a
// screenshot. Writes into ../docs/img. Does not touch git.
import { chromium } from "playwright";
import { fileURLToPath } from "node:url";
import { dirname, resolve } from "node:path";

const HERE = dirname(fileURLToPath(import.meta.url));
const IMG = resolve(HERE, "..", "..", "docs", "img");
const URL = "http://localhost:5273/";
const VIEWPORT = { width: 1680, height: 1000 };

// Flowsheets with no template in the gallery are injected via the app's autosave
// key (the store restores localStorage['caldyr.autosave.v1'] on load).
const RIGOROUS_COLUMN = {
  schema: "caldyr.flow/1",
  components: [{ id: "benzene" }, { id: "toluene" }],
  property_package: "thermo:PR",
  units: [{ id: "COL", type: "RigorousColumn", xy: [420, 220],
    params: { n_stages: 16, feed_stage: 8, reflux_ratio: 1.545, distillate_rate: 50.0, P: 101325.0 } }],
  streams: [
    { id: "FEED", from: null, to: "COL:in1",
      spec: { T: 365.0, P: 101325.0, molar_flow: 100.0, z: { benzene: 0.5, toluene: 0.5 } } },
    { id: "DIST", from: "COL:distillate", to: null },
    { id: "BOT", from: "COL:bottoms", to: null },
    { id: "QC", from: "COL:condenser_duty", to: null },
    { id: "QR", from: "COL:reboiler_duty", to: null },
  ],
};
const FLASH_RECYCLE = {
  schema: "caldyr.flow/1",
  components: [{ id: "n-pentane" }, { id: "n-octane" }],
  property_package: "thermo:PR",
  units: [
    { id: "MIX", type: "Mixer", params: { dP: 0.0 }, xy: [120, 200] },
    { id: "FL", type: "Flash", params: { T: 360.0, P: 101325.0 }, xy: [360, 200] },
    { id: "SP", type: "Splitter", params: { split: 0.6 }, xy: [600, 200] },
  ],
  streams: [
    { id: "FEED", from: null, to: "MIX:in1",
      spec: { T: 330.0, P: 101325.0, molar_flow: 10.0, z: { "n-pentane": 0.5, "n-octane": 0.5 } } },
    { id: "MIXOUT", from: "MIX:out", to: "FL:in1" },
    { id: "VAP", from: "FL:vapor", to: null },
    { id: "LIQ", from: "FL:liquid", to: "SP:in1" },
    { id: "RECY", from: "SP:out1", to: "MIX:in2" },
    { id: "BOT", from: "SP:out2", to: null },
    { id: "Q", from: "FL:duty", to: null },
  ],
};

const statusMatches = (re) =>
  `() => { const e = document.querySelector('[role=status]'); return e && ${re}.test(e.textContent); }`;

async function newPage(browser, autosave) {
  const ctx = await browser.newContext({ viewport: VIEWPORT, deviceScaleFactor: 2 });
  await ctx.addInitScript((seed) => {
    try {
      localStorage.setItem("caldyr.tour_seen", "1");
      if (seed) localStorage.setItem("caldyr.autosave.v1", seed);
    } catch { /* private mode */ }
  }, autosave ? JSON.stringify(autosave) : null);
  const page = await ctx.newPage();
  page.on("dialog", (d) => d.accept().catch(() => {}));   // "New" native confirm
  await page.goto(URL, { waitUntil: "domcontentloaded", timeout: 30000 });
  await page.waitForFunction(statusMatches("/Ready|Loaded/"), null, { timeout: 60000 });
  await page.waitForTimeout(400);
  // an empty canvas auto-opens the template gallery; dismiss any startup dialog
  if (await page.getByRole("dialog").isVisible().catch(() => false)) {
    await page.keyboard.press("Escape");
    await page.getByRole("dialog").waitFor({ state: "detached", timeout: 5000 }).catch(() => {});
  }
  return page;
}

// Clip the right inspector panel from the Params-tab left edge down to the bottom
// of its rendered content (leaf text / svg / canvas), so there's no empty tail.
async function inspectorClip(page, marginBottom = 26) {
  const pbox = await page.getByRole("tab", { name: "Params" }).boundingBox();
  const x = Math.max(0, pbox.x - 10);
  const y = pbox.y - 6;
  const bottom = await page.evaluate((panelX) => {
    let max = 0;
    for (const el of document.querySelectorAll("*")) {
      const r = el.getBoundingClientRect();
      if (r.width <= 0 || r.height <= 0 || r.left < panelX - 6 || r.top > window.innerHeight) continue;
      const leaf = el.children.length === 0 && el.textContent.trim();
      const gfx = el.tagName === "svg" || el.tagName === "CANVAS";
      if (leaf || gfx) max = Math.max(max, r.bottom);
    }
    return Math.min(window.innerHeight, max);
  }, x);
  return { x, y, width: VIEWPORT.width - x, height: Math.min(VIEWPORT.height, bottom + marginBottom) - y };
}

// Activate an inspector tab by firing the DOM click on the role=tab element —
// a plain Playwright .click() intermittently doesn't switch the panel.
async function clickTab(page, name) {
  await page.evaluate((n) => {
    const t = [...document.querySelectorAll("[role=tab]")].find((x) => x.textContent.trim() === n);
    if (t) t.click();
  }, name);
  await page.waitForTimeout(250);
}

// Activate a tab and confirm its panel rendered (retry: a trailing solve update
// can reset the tab back to "streams" right after a solve completes).
async function activatePanel(page, name, sentinelLabel) {
  for (let i = 0; i < 6; i++) {
    await clickTab(page, name);
    try {
      await page.getByLabel(sentinelLabel).first().waitFor({ state: "visible", timeout: 900 });
      return;
    } catch { await page.waitForTimeout(400); }
  }
  throw new Error(`could not activate the ${name} panel`);
}

async function shootAmmonia(browser, canvasOnly = false) {
  const page = await newPage(browser, null);
  await page.getByRole("button", { name: "Projects" }).click();
  const dialog = page.getByRole("dialog");
  await dialog.waitFor({ state: "visible", timeout: 10000 });
  await page.getByRole("button", { name: "Ammonia loop" }).click();
  // the template button closes the dialog; wait for it, Escape as a fallback
  await dialog.waitFor({ state: "detached", timeout: 5000 }).catch(async () => {
    await page.keyboard.press("Escape");
    await dialog.waitFor({ state: "detached", timeout: 5000 });
  });
  await page.waitForTimeout(400);
  await page.getByRole("radio", { name: "PFD" }).click();
  await page.getByRole("button", { name: "Solve", exact: true }).click();
  await page.waitForFunction(statusMatches("/^Solved/"), null, { timeout: 90000 });
  await page.waitForTimeout(400);
  await page.getByRole("button", { name: "Arrange" }).click().catch(() => {});
  await page.waitForTimeout(1500);
  // Apply phase coloring AFTER the solve so it uses the solved stream data — set
  // before Solve it stays gray. Toggle none→phase to force a recolor, then wait
  // until the legend's "solve to color streams" hint disappears (= streams colored).
  const colorSel = page.locator('select[title^="Color streams"]');
  await colorSel.selectOption("none");
  await page.waitForTimeout(200);
  await colorSel.selectOption("phase");
  await page.waitForFunction(
    "() => !document.querySelector('.react-flow') || " +
    "![...document.querySelectorAll('.react-flow *')].some(e => /solve to color/i.test(e.textContent || ''))",
    null, { timeout: 8000}).catch(() => {});
  await page.waitForTimeout(500);
  await page.locator(".react-flow__controls-fitview").click().catch(() => {});
  await page.waitForTimeout(700);
  // Hero shot: hide the minimap/controls/attribution, then clip tight to the
  // flowsheet's node bounding box + margin. The ammonia loop is an ⌐-shape (a
  // horizontal train with the separator/splitter dropping on the right), leaving
  // its lower-left corner empty; we tuck the phase legend into that empty corner
  // so it's kept in frame without inflating the crop.
  const clip = await page.evaluate((margin) => {
    for (const sel of [".react-flow__minimap", ".react-flow__controls", ".react-flow__attribution"]) {
      const el = document.querySelector(sel);
      if (el) el.style.display = "none";
    }
    const pane = document.querySelector(".react-flow").getBoundingClientRect();
    const nodes = [...document.querySelectorAll(".react-flow__node")].map((n) => n.getBoundingClientRect());
    const nx0 = Math.min(...nodes.map((r) => r.left)), ny0 = Math.min(...nodes.map((r) => r.top));
    const nx1 = Math.max(...nodes.map((r) => r.right)), ny1 = Math.max(...nodes.map((r) => r.bottom));
    // place the phase legend inside the empty lower-left of the flowsheet's bbox
    const legend = [...document.querySelectorAll(".react-flow__panel, .react-flow *")]
      .find((e) => e.children.length <= 6 && /two-phase|vapor/i.test(e.textContent || ""));
    if (legend) {
      legend.style.bottom = "auto";
      legend.style.right = "auto";
      legend.style.left = `${nx0 - pane.left + 8}px`;
      legend.style.top = `${ny0 - pane.top + (ny1 - ny0) * 0.72}px`;
    }
    const x0 = Math.max(pane.left, nx0 - margin);
    const y0 = Math.max(pane.top, ny0 - margin);
    const x1 = Math.min(pane.right, nx1 + margin);
    const y1 = Math.min(pane.bottom, ny1 + margin);
    return { x: x0, y: y0, width: x1 - x0, height: y1 - y0 };
  }, 60);
  await page.screenshot({ path: `${IMG}/ammonia-loop-solved.png`, clip });
  if (canvasOnly) { await page.context().close(); console.log("ammonia: 1 shot (canvas)"); return; }

  await page.getByRole("button", { name: "Cost", exact: true }).click();
  await page.waitForFunction(statusMatches("/^Costed/"), null, { timeout: 90000 });
  await clickTab(page, "Econ");
  await page.waitForTimeout(900);
  await page.screenshot({ path: `${IMG}/ammonia-econ.png`, clip: await inspectorClip(page) });

  await page.getByRole("button", { name: "Run 500 samples" }).click();
  await page.waitForFunction(statusMatches("/^Costed/"), null, { timeout: 120000 });
  await page.waitForTimeout(1200);
  await page.screenshot({ path: `${IMG}/ammonia-econ-mc.png`, clip: await inspectorClip(page) });
  await page.context().close();
  console.log("ammonia: 3 shots");
}

async function shootRigorous(browser) {
  const page = await newPage(browser, RIGOROUS_COLUMN);
  await page.getByRole("button", { name: "Solve", exact: true }).click();
  await page.waitForFunction(statusMatches("/^Solved/"), null, { timeout: 120000 });
  await page.waitForTimeout(400);
  await page.locator(".react-flow__node", { hasText: "COL" }).first().click();
  await page.waitForTimeout(300);
  await clickTab(page, "Params");
  await page.waitForTimeout(300);
  // scroll the inspector's scroll container to the bottom so both charts show
  await page.evaluate(() => {
    const t = [...document.querySelectorAll("*")].filter(
      (e) => e.children.length === 0 && e.textContent.trim() === "Liquid composition profile");
    let el = t[t.length - 1]; if (!el) return;
    for (let p = el.parentElement; p; p = p.parentElement) {
      const s = getComputedStyle(p);
      if (/(auto|scroll)/.test(s.overflowY) && p.scrollHeight > p.clientHeight) { p.scrollTop = p.scrollHeight; return; }
    }
  });
  await page.waitForTimeout(700);
  await page.screenshot({ path: `${IMG}/distillation-design-results.png`, clip: await inspectorClip(page) });
  await page.context().close();
  console.log("rigorous column: 1 shot");
}

async function shootOptimize(browser) {
  const page = await newPage(browser, FLASH_RECYCLE);
  await page.getByRole("button", { name: "Solve", exact: true }).click();
  await page.waitForFunction(statusMatches("/^Solved/"), null, { timeout: 90000 });
  await page.waitForTimeout(1200);   // let trailing solve updates settle
  await activatePanel(page, "Opt", "Sense");
  // helper: choose the first option matching a predicate (boundary streams get
  // auto-renamed on load — e.g. VAP -> PROD_VAP — and duties key as <unit>_duty)
  const optionMatching = (loc, re) =>
    loc.evaluate((s, r) => [...s.options].map((o) => o.value).find((v) => v && new RegExp(r).test(v)), re.source);

  await page.getByLabel("Sense").selectOption("min");
  await page.getByLabel("Metric type").first().selectOption("duty");
  await page.waitForTimeout(150);
  const dutyOpt = await optionMatching(page.getByLabel("Duty"), /duty|FL/i);
  await page.getByLabel("Duty").selectOption(dutyOpt);
  await page.getByRole("button", { name: "design variable" }).click();
  await page.waitForTimeout(250);
  await page.getByLabel("Unit", { exact: true }).selectOption("FL");
  await page.waitForTimeout(200);
  await page.getByLabel("Parameter", { exact: true }).selectOption("T");
  await page.getByLabel("Lower bound").fill("340");
  await page.getByLabel("Upper bound").fill("370");
  await page.getByRole("button", { name: "constraint" }).click();
  await page.waitForTimeout(250);
  await page.getByLabel("Metric type").last().selectOption("component_rate");
  await page.waitForTimeout(200);
  const vapOpt = await optionMatching(page.getByLabel("Stream").last(), /VAP/i);
  await page.getByLabel("Stream").last().selectOption(vapOpt);
  await page.getByLabel("Component").last().selectOption("n-pentane");
  await page.getByLabel("Operator").selectOption(">=");
  await page.getByLabel("Constraint value").fill("4.2");
  await page.waitForTimeout(200);
  await page.getByRole("button", { name: "Optimize" }).click();
  await page.waitForSelector("text=engine solves", { timeout: 120000 });
  await page.waitForTimeout(600);
  await page.screenshot({ path: `${IMG}/optimization-opt-panel.png`, clip: await inspectorClip(page) });
  await page.context().close();
  console.log("optimize panel: 1 shot");
}

// Optional arg selects a subset: "canvas" (hero shot only), "ammonia",
// "rigorous", or "optimize". No arg captures everything.
const only = process.argv[2] || "";
const browser = await chromium.launch({ headless: true });
try {
  if (only === "canvas") await shootAmmonia(browser, true);
  else {
    if (!only || only === "ammonia") await shootAmmonia(browser);
    if (!only || only === "rigorous") await shootRigorous(browser);
    if (!only || only === "optimize") await shootOptimize(browser);
  }
  console.log(`Done → ${IMG}`);
} finally {
  await browser.close();
}
