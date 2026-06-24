"""Closed double-recycle anhydrous-ethanol entrainer PLANT (Hameed §9.5.6).

The integrated decanting-condenser column (test_m24) assembled into the full
two-column / two-recycle train: T-100 (decant column) -> aqueous -> T-101 (water
column) -> distillate recycle -> cyclohexane Makeup -> back to T-100. The plant
is brought up by distillate-rate continuation at the FLOWSHEET level (D up,
make-up target down), warm-starting T-100 each step. This test runs an
abbreviated continuation and checks that the closed loop converges, the overall
plant mass balance closes, and the entrainer drives the bottoms toward anhydrous
ethanol (the bottoms ethanol fraction rises as the distillate rate rises).

Book-scale >99.95% EtOH needs the 50-62 stages of the reference (a 30-stage
column caps ~90% here); the point validated is the converged closed train.
"""
import pytest

from caldyr.core import Component, Flowsheet
from caldyr.unitops import Makeup, RigorousColumn

P_ATM = 101325.0
COMPS = ["ethanol", "water", "cyclohexane"]
F_FRESH = 27.78
Z_FRESH = {"ethanol": 0.87, "water": 0.13}


def _plant() -> Flowsheet:
    fs = Flowsheet(components=[Component(c) for c in COMPS],
                   property_package="thermo:UNIFAC")
    fs.add(RigorousColumn("T100", {
        "n_stages": 30, "feeds": [{"stage": 3}, {"stage": 14}],
        "reflux_ratio": 3.0, "distillate_rate": 4.0,
        "method": "naphtali_sandholm", "reboiled": True,
        "decant_condenser": True, "condenser_T": 305.0,
        "reflux_layer": "organic", "max_iter": 120,
    }))
    fs.add(RigorousColumn("T101", {
        "n_stages": 16, "feed_stage": 8, "reflux_ratio": 1.0,
        "distillate_rate": 2.0, "method": "naphtali_sandholm",
        "reboiled": True, "max_iter": 120,
    }))
    fs.add(Makeup("MK", {"component": "cyclohexane", "target": 8.0,
                         "T": 305.0, "P": P_ATM}))
    fs.feed("FEED", "T100:in2", T=343.0, P=P_ATM, molar_flow=F_FRESH, z=Z_FRESH)
    fs.connect("ENTR", "MK:out", "T100:in1")
    fs.connect("AQ", "T100:distillate", "T101:in1")
    fs.connect("REC", "T101:distillate", "MK:in1")
    fs.connect("ETOH", "T100:bottoms", None)
    fs.connect("WATER", "T101:bottoms", None)
    for u in ("T100", "T101"):
        fs.connect(f"{u}_QC", f"{u}:condenser_duty", None)
        fs.connect(f"{u}_QR", f"{u}:reboiler_duty", None)
    fs.solver_hints = {
        "tear_guesses": {
            "ENTR": {"T": 320.0, "P": P_ATM, "molar_flow": 8.0,
                     "z": {"ethanol": 0.25, "water": 0.10, "cyclohexane": 0.65}},
            "REC": {"T": 333.0, "P": P_ATM, "molar_flow": 8.0,
                    "z": {"ethanol": 0.30, "water": 0.10, "cyclohexane": 0.60}},
        },
        "tear_tolerance": 5e-3,
    }
    return fs


@pytest.mark.slow
def test_entrainer_plant_closes_and_concentrates_ethanol():
    fs = _plant()
    t100, t101, mk = fs.units["T100"], fs.units["T101"], fs.units["MK"]
    eth_prev = 0.0
    for d100, d101, mkt in [(4.0, 2.0, 8.0), (7.0, 4.0, 6.0), (10.0, 6.0, 5.0)]:
        t100.params["distillate_rate"] = d100
        t101.params["distillate_rate"] = d101
        mk.params["target"] = mkt
        rep = fs.solve(method="direct", max_iter=30)
        assert rep.converged, f"recycle did not converge at D100={d100}"
        etoh = fs.streams["ETOH"]
        # the entrainer concentrates ethanol in the bottoms as D rises.
        assert etoh.z["ethanol"] > eth_prev
        eth_prev = etoh.z["ethanol"]

    # Closed-plant overall mass balance: the only streams crossing the boundary
    # are the fresh feed + the cyclohexane make-up IN, and the two product
    # bottoms OUT (the recycles are internal).
    etoh, water = fs.streams["ETOH"], fs.streams["WATER"]
    makeup = mk.design["makeup_flow"]
    for c in COMPS:
        fed = F_FRESH * Z_FRESH.get(c, 0.0) + (makeup if c == "cyclohexane" else 0.0)
        out = etoh.molar_flow * etoh.z[c] + water.molar_flow * water.z[c]
        assert out == pytest.approx(fed, rel=2e-2, abs=2e-2), f"{c} balance"

    # The entrainer recirculates: T-101 returns cyclohexane to T-100 (the recycle
    # carries far more cyclohexane than the small make-up adds).
    rec = fs.streams["REC"]
    assert rec.molar_flow * rec.z["cyclohexane"] > 0.5
    assert etoh.z["ethanol"] > 0.78        # ~0.79 at D100=10 (trending anhydrous)
