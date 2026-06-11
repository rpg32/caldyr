"""FastAPI bridge over the Caldyr engine (thin transport; no physics here).

Adds the engine source dir to ``sys.path`` so ``import caldyr`` works without an
editable install (mirrors the repo ``conftest.py`` / example bootstrap).
"""
import sys
from pathlib import Path

_ENGINE = Path(__file__).resolve().parent.parent / "engine"
if _ENGINE.is_dir() and str(_ENGINE) not in sys.path:
    sys.path.insert(0, str(_ENGINE))
