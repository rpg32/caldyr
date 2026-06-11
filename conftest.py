"""Make the engine importable in tests without an editable install.

`caldyr` lives under ``engine/``; add it to ``sys.path`` so ``pytest`` from the
repo root can ``import caldyr`` directly. An editable install (``pip install -e
engine``) also works and takes precedence.
"""
import sys
from pathlib import Path

ENGINE = Path(__file__).parent / "engine"
if ENGINE.is_dir() and str(ENGINE) not in sys.path:
    sys.path.insert(0, str(ENGINE))
