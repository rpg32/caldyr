# Docs image tools

The six images embedded in the tutorials (`docs/img/`) are all reproducible from
these two scripts — no manual screenshotting or cropping. Five are UI screenshots
driven with headless Playwright; one is a matplotlib data plot from the engine.

## Prerequisites

- **For the UI screenshots** (`shoot-docs.mjs`): both dev servers running, and the
  web app's Playwright/Chromium installed.
  ```bash
  # terminal 1 — engine API
  python -m uvicorn api.main:app --port 8753        # (Python 3.11 env)
  # terminal 2 — web app
  cd web && npm run dev                             # serves http://localhost:5273
  ```
- **For the reflux chart** (`reflux_chart.py`): the `caldyr` engine importable plus
  `matplotlib`. On this dev box the engine deps live in Python 3.11
  (`C:\Python311\python.exe`), not the default interpreter.

## The two commands

```bash
# 1. five UI screenshots  → docs/img/ (run from web/)
node scripts/shoot-docs.mjs            # all five; or: canvas | ammonia | rigorous | optimize

# 2. the reflux trade-off chart → docs/img/distillation-reflux-tradeoff.png
python web/scripts/reflux_chart.py     # (use the Python 3.11 that has caldyr)
```

`shoot-docs.mjs` takes an optional argument to regenerate a subset:
`canvas` (the README hero only), `ammonia`, `rigorous`, or `optimize`. No argument
captures all five.

## Which script owns which image

| Image | Script | What it shows |
|---|---|---|
| `ammonia-loop-solved.png` | `shoot-docs.mjs` (`canvas`) | Solved ammonia loop, PFD view, phase-colored — the README/tutorial hero |
| `ammonia-econ.png` | `shoot-docs.mjs` (`ammonia`) | Econ tab: LCOP/TCI/OPEX/NPV, equipment, tornado |
| `ammonia-econ-mc.png` | `shoot-docs.mjs` (`ammonia`) | Econ tab after Monte-Carlo: tornado + P10/P50/P90 histogram |
| `distillation-design-results.png` | `shoot-docs.mjs` (`rigorous`) | RigorousColumn Design-results panel + stage-profile charts |
| `optimization-opt-panel.png` | `shoot-docs.mjs` (`optimize`) | Opt builder filled in + optimization result |
| `distillation-reflux-tradeoff.png` | `reflux_chart.py` | Capital vs utilities vs reflux ratio (matplotlib) |

## Notes

- `shoot-docs.mjs` drives the real app: it loads the ammonia template (or injects a
  flowsheet via the autosave key for the cases with no gallery template), solves,
  opens the relevant tab/panel, and clips each shot to its content. The hero shot
  hides the minimap/controls, tucks the phase legend into the loop's empty corner,
  and crops to the flowsheet's bounding box + margin.
- The screenshots are deterministic except the Monte-Carlo histogram (random
  sampling) — re-running changes only that one image's bars.
- CDP/`Page.captureScreenshot` via the Chrome extension is unreliable on this
  machine; that's why these use Playwright directly (see `docs/UX_BACKLOG.md`).
