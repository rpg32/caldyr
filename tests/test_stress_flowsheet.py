"""Industrial-scale stress case: a 31-unit, 3-recycle plant (three parallel
ammonia synthesis trains + shared product workup). Guards solver scaling — the
roadmap risk was "unproven beyond ~10 units / ~5 tears" — and doubles as a
regression net for cross-unit interactions at scale.

Wall-clock observations on first landing (Windows, i7-class):
sequential ~3-8 s (dominated by one-time flasher build), EO comparable.
The assertions use generous ceilings so CI variance doesn't flake.
"""
from __future__ import annotations

import time

from caldyr.core import Component, Flowsheet
from caldyr.solver import balance_report
from caldyr.unitops import (
    Compressor,
    FlashDrum,
    Heater,
    Mixer,
    Pump,
    Splitter,
    Valve,
)
from caldyr.unitops.equilibrium_reactor import EquilibriumReactor

AMMONIA_RXN = {"stoich": {"nitrogen": -1, "hydrogen": -3, "ammonia": 2}, "key": "nitrogen"}
FEED_Z = {"nitrogen": 0.2475, "hydrogen": 0.7425, "ammonia": 0.0, "argon": 0.01}


def build_plant() -> Flowsheet:
    fs = Flowsheet(
        components=[Component(c) for c in ("nitrogen", "hydrogen", "ammonia", "argon")],
        property_package="thermo:PR",
    )

    # -- three parallel synthesis trains, each with its own recycle ---------
    for t in (1, 2, 3):
        fs.add(Mixer(f"MIX{t}", {"dP": 0.0}))
        fs.add(Compressor(f"COMP{t}", {"P_out": 2.0e7, "eta": 0.78}))
        fs.add(Heater(f"PRE{t}", {"T_out": 673.15}))
        fs.add(EquilibriumReactor(f"RXN{t}", {"reaction": AMMONIA_RXN, "T": 673.15}))
        fs.add(Heater(f"COOL{t}", {"T_out": 250.0}))
        fs.add(FlashDrum(f"SEP{t}", {"T": 250.0, "P": 2.0e7}))
        fs.add(Splitter(f"SPL{t}", {"split": 0.9}))

        fs.feed(f"FEED{t}", f"MIX{t}:in1", T=300.0, P=1.0e7,
                molar_flow=100.0, z=dict(FEED_Z))
        fs.connect(f"T{t}_S1", f"MIX{t}:out", f"COMP{t}:in1")
        fs.connect(f"T{t}_S2", f"COMP{t}:out", f"PRE{t}:in1")
        fs.connect(f"T{t}_S3", f"PRE{t}:out", f"RXN{t}:in1")
        fs.connect(f"T{t}_S4", f"RXN{t}:out", f"COOL{t}:in1")
        fs.connect(f"T{t}_S5", f"COOL{t}:out", f"SEP{t}:in1")
        fs.connect(f"T{t}_VAP", f"SEP{t}:vapor", f"SPL{t}:in1")
        fs.connect(f"T{t}_RECY", f"SPL{t}:out1", f"MIX{t}:in2")   # tear t
        fs.connect(f"T{t}_PURGE", f"SPL{t}:out2", None)

    # -- shared liquid-product workup (units 22..31) -------------------------
    fs.add(Mixer("PMIX1", {"dP": 0.0}))
    fs.add(Mixer("PMIX2", {"dP": 0.0}))
    fs.connect("L1", "SEP1:liquid", "PMIX1:in1")
    fs.connect("L2", "SEP2:liquid", "PMIX1:in2")
    fs.connect("L12", "PMIX1:out", "PMIX2:in1")
    fs.connect("L3", "SEP3:liquid", "PMIX2:in2")

    fs.add(Valve("LETDOWN", {"P_out": 5.0e5}))
    fs.add(FlashDrum("DEGAS", {"P": 5.0e5}))          # adiabatic letdown flash
    fs.add(Heater("CHILL", {"T_out": 240.0}))
    fs.add(FlashDrum("POLISH", {"T": 240.0, "P": 5.0e5}))
    fs.add(Pump("PROD_PUMP", {"P_out": 2.0e6, "eta": 0.7}))
    fs.add(Compressor("GAS_COMP", {"P_out": 2.0e6, "eta": 0.75}))
    fs.add(Heater("GAS_COOL", {"T_out": 310.0}))

    fs.connect("W1", "PMIX2:out", "LETDOWN:in1")
    fs.connect("W2", "LETDOWN:out", "DEGAS:in1")
    fs.connect("W3", "DEGAS:liquid", "CHILL:in1")
    fs.connect("W4", "CHILL:out", "POLISH:in1")
    fs.connect("NH3_PROD", "POLISH:liquid", "PROD_PUMP:in1")
    fs.connect("NH3_HP", "PROD_PUMP:out", None)
    fs.connect("OFFGAS1", "DEGAS:vapor", "GAS_COMP:in1")
    fs.connect("OFFGAS2", "GAS_COMP:out", "GAS_COOL:in1")
    fs.connect("FUELGAS", "GAS_COOL:out", None)
    fs.connect("OFFGAS3", "POLISH:vapor", None)

    return fs


def test_plant_has_industrial_scale():
    fs = build_plant()
    assert len(fs.units) == 30
    assert len(fs.connections) >= 36


def test_sequential_solves_three_recycles_at_scale():
    fs = build_plant()
    t0 = time.perf_counter()
    report = fs.solve(backend="sequential", tol=1e-6)
    wall = time.perf_counter() - t0
    assert report.converged, report.messages
    assert sorted(report.tear_streams) == ["T1_RECY", "T2_RECY", "T3_RECY"]
    assert wall < 120.0, f"sequential solve took {wall:.1f}s"

    # mass closes plant-wide; every train actually made ammonia
    bal = balance_report(fs)
    assert bal["overall"]["mass_rel"] < 1e-5
    prod = fs.streams["NH3_PROD"]
    assert prod.molar_flow is not None and prod.molar_flow > 60.0
    assert prod.normalized_z()["ammonia"] > 0.95


def test_equation_oriented_agrees_at_scale():
    fs_sm = build_plant()
    fs_sm.solve(backend="sequential", tol=1e-8)

    fs_eo = build_plant()
    t0 = time.perf_counter()
    report = fs_eo.solve(backend="equation_oriented", tol=1e-8)
    wall = time.perf_counter() - t0
    assert report.converged, report.messages
    assert wall < 300.0, f"EO solve took {wall:.1f}s"

    for sid in ("NH3_PROD", "FUELGAS", "T1_RECY", "T3_RECY"):
        a, b = fs_sm.streams[sid], fs_eo.streams[sid]
        assert a.molar_flow is not None and b.molar_flow is not None
        assert abs(a.molar_flow - b.molar_flow) < 1e-3 * max(1.0, a.molar_flow)
        assert a.T is not None and b.T is not None and abs(a.T - b.T) < 0.01
