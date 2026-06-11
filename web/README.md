# web — Caldyr canvas (React + @xyflow/react)

A thin client over the engine API. **No physics in the browser** — it calls the
API, which calls the engine.

- Canvas: `@xyflow/react` (React Flow v12). Nodes = unit ops (one handle per
  Port; energy ports amber), edges = streams. Boundary feeds/products are
  UI-only nodes.
- `canvasToFlow` / `flowToCanvas` make the canvas and the `.flow` document the
  same object — flowsheet-as-canvas == flowsheet-as-code.
- Panels: unit/feed params, the solved stream table, and an economics tab with a
  capital/opex/LCOP summary and a sensitivity tornado.

## Run

Start the engine API (from the repo root), then the web dev server:

```bash
# 1) API on :8753
python -m uvicorn api.main:app --port 8753

# 2) web on :5273 (proxies /api -> :8753)
cd web && npm install && npm run dev
```

Open http://localhost:5273, click **Ammonia loop** (or **Mixer + Heater**), then
**Solve** and **Cost**. Ports are pinned (`strictPort`) to 5273/5274 and the API
to 8753 to avoid clashing with other local apps.

## What's built vs. future

Built: a single working flowsheet canvas — the M5 DoD (build + solve + cost in
the browser). Future: the BFD/PFD/P&ID three-view toggle, live plots, `.flow`
file load/save, and richer param editing.
