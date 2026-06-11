# Contributing

- Read `CLAUDE.md`, `ARCHITECTURE.md`, `docs/DATA_MODEL.md` first.
- Engine code is Python 3.11+, fully typed, `ruff`-formatted, `pytest`-tested.
- **Every unit op / property method ships with a test validated against a cited
  reference** (textbook case, DWSIM, or Aspen). Put the citation in the test.
- SI units internally; convert only at I/O edges.
- Small commits (Conventional Commits). The engine must stay runnable headless.
- The web app contains no physics — it calls the engine API.

## Quality gates (CI enforces these)

```bash
pip install -e ".[dev,api,ai,kinetics]"
python -m pytest -q              # full suite, all green
python -m ruff check engine api tests
python -m mypy engine/caldyr
cd web && npm ci && npm run build && npm test
```

Playwright e2e (`cd web && npx playwright test`) needs the API on :8753 and is
run locally before releases.

## Adding a unit operation

1. New module in `engine/caldyr/unitops/`, registered in the REGISTRY; ports
   via the common contract; typed errors for infeasible specs.
2. A sizer entry (registry in `economics/sizing.py`) + Turton constants in
   `economics/data.py` **with citations**, so `cost` works.
3. Tests citing your validation source; an example if it demonstrates
   something new.
4. The web palette, AI tools, and `.flow` schema pick it up automatically.

## License hygiene

Never copy or statically link GPL code (e.g. DWSIM internals). Interop across
process/API boundaries only.
