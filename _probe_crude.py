"""Scratch probe: flash timing on the crude slate + water; NOT a deliverable."""
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent / "engine"))

from caldyr.assay import api_to_sg, characterize_assay
from caldyr.core.components_db import resolve_component
from caldyr.thermo import make_package


def degF(f):
    return (f - 32.0) / 1.8 + 273.15


TBP = [(0, 80.0), (10, 255.0), (20, 349.0), (30, 430.0), (40, 527.0),
       (50, 635.0), (60, 751.0), (70, 915.0), (80, 1095.0),
       (90, 1277.0), (98, 1410.0)]
MW = [(0, 68.0), (10, 119.0), (20, 150.0), (30, 182.0), (40, 225.0),
      (50, 282.0), (60, 350.0), (70, 456.0), (80, 585.0),
      (90, 713.0), (98, 838.0)]
API = [(13, 63.28), (33, 54.86), (57, 45.91), (74, 38.21), (91, 26.01)]
LE = {"isobutane": 0.19, "n-butane": 0.11, "isopentane": 0.37, "n-pentane": 0.46}

t0 = time.perf_counter()
res = characterize_assay([(v, degF(f)) for v, f in TBP], kind="TBP",
                         api_gravity=48.75,
                         sg_curve=[(v, api_to_sg(a)) for v, a in API],
                         mw_curve=MW, n_cuts=12, light_ends=LE)
print(f"characterize: {time.perf_counter() - t0:.2f}s, {len(res.cuts)} cuts")

comps = [resolve_component("water")] + res.components()
ids = [c.id for c in comps]
t0 = time.perf_counter()
pp = make_package("thermo:PR", ids)
print(f"make_package({len(ids)} comps): {time.perf_counter() - t0:.2f}s")

z = res.mole_fractions()
z = {**{c: v * 0.92 for c, v in z.items()}, "water": 0.08}

# typical stage liquid: crude + a bit of water
t0 = time.perf_counter()
r = pp.flash_pt(degF(650.0), 448159.0, z)
print(f"flash_pt: {time.perf_counter() - t0:.3f}s VF={r.vapor_fraction:.3f}")

t0 = time.perf_counter()
b = pp.bubble_point(202000.0, z)
dt = time.perf_counter() - t0
print(f"bubble_point(full feed + water): {dt:.3f}s T={b.T:.1f}K y_water={b.y['water']:.4f}")

# a naphtha-ish liquid with water
zn = {ids[i]: v for i, v in enumerate([0.15, 0.01, 0.01, 0.02, 0.03, 0.30,
                                       0.25, 0.15, 0.05, 0.03])}
tot = sum(zn.values())
zn = {c: v / tot for c, v in zn.items()}
t0 = time.perf_counter()
b = pp.bubble_point(135800.0, zn)
print(f"bubble_point(naphtha+water): {time.perf_counter() - t0:.3f}s T={b.T:.1f}K")

# repeat timing (cached constants)
t0 = time.perf_counter()
for _ in range(10):
    pp.bubble_point(202000.0, z)
print(f"10x bubble_point: {time.perf_counter() - t0:.3f}s")

# heavy resid-ish liquid with trace water
zh = {c.id: f for c, f in zip(res.cuts[-4:], [0.2, 0.3, 0.3, 0.19])}
zh["water"] = 0.01
t0 = time.perf_counter()
b = pp.bubble_point(225500.0, zh)
print(f"bubble_point(resid+water): {time.perf_counter() - t0:.3f}s T={b.T:.1f}K")
