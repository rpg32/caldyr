"""M21 tests: the multi-stream (LNG / plate-fin) heat exchanger.

The LNG operation (Hameed 2025, *Chemical Process Simulations using Aspen
HYSYS*, sec. 9.5.2) "solves heat and material balances for multi-stream heat
exchangers": an arbitrary number of passes share one core, each pass's hot/cold
role is decided by the solution (not a label), and the **Weighted** method
splits the heating curve into intervals to find the minimum internal approach
and the conductance.

The book embeds its LNG-100 inside a full turbo-expander LPG plant (expander +
separators + recycle + a recovery column) and reports no standalone exchanger
numbers, so the validation here is on the physics the operation must obey:

* the two-pass form reproduces the simple :class:`HeatExchanger` exactly (same
  duty, same outlet temperatures, and a required ``UA`` equal to ``Q/LMTD``);
* the overall energy balance closes to machine precision over N passes;
* the global **minimum-approach** and **UA** specs are met by the iterative
  solve (the "multiple unknowns" case);
* a phase change on a pass is captured by the zone (weighted) analysis;
* an infeasible temperature cross is rejected with a typed error rather than a
  silently wrong answer.
"""
import math

import pytest

from caldyr.core import Component, Flowsheet, Stream
from caldyr.economics.costing import cost_equipment
from caldyr.economics.sizing import size_flowsheet
from caldyr.io import from_dict, to_dict
from caldyr.thermo import make_package
from caldyr.unitops import (
    HeatExchanger,
    MultiStreamExchanger,
    MultiStreamExchangerError,
)

_COMPS = ["methane", "ethane", "propane"]
_Z = {"methane": 0.9, "ethane": 0.07, "propane": 0.03}
_P = 3.0e6


def _pkg(comps=_COMPS):
    return make_package("thermo:PR", comps)


def _stream(T, n=1.0, P=_P, z=_Z, comps=_COMPS):
    return Stream(id="s", components=comps, T=T, P=P, molar_flow=n, z=z)


# -- two-pass form reproduces the simple HeatExchanger ------------------------
def test_two_pass_matches_heat_exchanger():
    pp = _pkg()
    hot, cold = _stream(350.0), _stream(250.0)
    Q = 1500.0

    hx = HeatExchanger("HX", {"duty": Q})
    ref = hx.solve({"hot_in": hot, "cold_in": cold}, pp)

    ms = MultiStreamExchanger("MS", {"passes": [{"duty": -Q}, {}]})
    out = ms.solve({"pass1_in": hot, "pass2_in": cold}, pp)

    # identical outlet temperatures
    assert abs(out["pass1_out"].T - ref["hot_out"].T) < 1e-6
    assert abs(out["pass2_out"].T - ref["cold_out"].T) < 1e-6
    # energy balance: exactly equal-and-opposite duties
    assert abs(sum(ms.design["pass_duties"])) < 1e-9
    assert abs(ms.design["pass_duties"][0] + Q) < 1e-9

    # required UA ~ Q/LMTD (the zone analysis refines the single-LMTD value by
    # capturing the real-gas Cp variation, so allow a small tolerance)
    dt1 = hot.T - ref["cold_out"].T
    dt2 = ref["hot_out"].T - cold.T
    ua_ref = Q / HeatExchanger.lmtd(dt1, dt2)
    assert abs(ms.design["UA"] - ua_ref) / ua_ref < 0.02
    # MITA equals the smaller end approach of a counter-current 2-stream unit
    assert abs(ms.design["min_approach"] - min(dt1, dt2)) < 0.2


# -- N-pass energy balance ----------------------------------------------------
def test_three_pass_energy_balance_closes():
    pp = _pkg()
    # one hot pass splits its heat to two cold passes; the third is free.
    ms = MultiStreamExchanger(
        "M3", {"passes": [{"T_out": 260.0}, {"T_out": 320.0}, {}]})
    ms.solve(
        {"pass1_in": _stream(380.0), "pass2_in": _stream(250.0),
         "pass3_in": _stream(260.0)}, pp)
    d = ms.design
    assert abs(sum(d["pass_duties"])) < 1e-6          # Sum Q = 0 (adiabatic)
    assert d["pass_duties"][0] < 0                     # pass 1 is the hot one
    assert d["pass_duties"][1] > 0 and d["pass_duties"][2] > 0
    assert d["min_approach"] > 0
    assert math.isfinite(d["UA"]) and d["UA"] > 0
    # hot duty released equals cold duty absorbed
    assert abs(d["hot_duty"] - d["cold_duty"]) < 1e-6


def test_hot_cold_designation_is_by_sign_not_label():
    # A pass given a T_out ABOVE its inlet is cold; below is hot — the unit
    # never asks for a hot/cold label (as HYSYS notes, a mislabelled pass still
    # solves). Pass 1 here is heated, pass 2 cooled, decided by the solution.
    pp = _pkg()
    ms = MultiStreamExchanger("MD", {"passes": [{"T_out": 330.0}, {}]})
    out = ms.solve({"pass1_in": _stream(250.0), "pass2_in": _stream(380.0)}, pp)
    assert out["pass1_out"].T > 250.0                  # pass 1 heated -> cold
    assert out["pass2_out"].T < 380.0                  # pass 2 cooled -> hot
    assert ms.design["pass_duties"][0] > 0
    assert ms.design["pass_duties"][1] < 0


# -- global specs (the iterative, multiple-unknown case) ----------------------
def test_global_min_approach_spec_is_met():
    pp = _pkg()
    ms = MultiStreamExchanger("MA", {"passes": [{}, {}], "min_approach": 15.0})
    ms.solve({"pass1_in": _stream(380.0), "pass2_in": _stream(250.0)}, pp)
    assert abs(ms.design["min_approach"] - 15.0) < 0.05
    assert abs(sum(ms.design["pass_duties"])) < 1e-6


def test_global_ua_spec_is_met():
    pp = _pkg()
    ms = MultiStreamExchanger("MU", {"passes": [{}, {}], "ua": 50.0})
    ms.solve({"pass1_in": _stream(380.0), "pass2_in": _stream(250.0)}, pp)
    assert abs(ms.design["UA"] - 50.0) < 0.1
    assert ms.design["min_approach"] > 0


def test_two_global_specs_rejected():
    pp = _pkg()
    ms = MultiStreamExchanger(
        "MB", {"passes": [{}, {}], "min_approach": 10.0, "ua": 50.0})
    with pytest.raises(MultiStreamExchangerError, match="at most one global"):
        ms.solve({"pass1_in": _stream(380.0), "pass2_in": _stream(250.0)}, pp)


def test_unreachable_min_approach_raises():
    # The cold inlet sits a hair below the hot outlet end, so the streams cannot
    # be held 30 K apart — the global solve fails cleanly, not silently.
    pp = _pkg()
    ms = MultiStreamExchanger(
        "MR", {"passes": [{"T_out": 255.0}, {}, {}], "min_approach": 30.0})
    with pytest.raises(MultiStreamExchangerError, match="could not satisfy"):
        ms.solve({"pass1_in": _stream(380.0), "pass2_in": _stream(250.0),
                  "pass3_in": _stream(255.0)}, pp)


# -- phase change captured by the zone analysis -------------------------------
def test_phase_change_is_captured():
    pp = _pkg()
    # a heavier (propane/ethane-rich) gas partially condenses as it cools at
    # 30 bar; the composite curve must bend through the dew point.
    zr = {"methane": 0.45, "ethane": 0.30, "propane": 0.25}
    hot = _stream(340.0, n=1.0, z=zr)
    cold = _stream(255.0, n=1.2, z=zr)
    ms = MultiStreamExchanger("PC", {"passes": [{"T_out": 268.0}, {}]})
    out = ms.solve({"pass1_in": hot, "pass2_in": cold}, pp)
    assert 0.0 < out["pass1_out"].vapor_fraction < 1.0   # really two-phase
    d = ms.design
    assert d["min_approach"] > 0
    assert math.isfinite(d["UA"]) and d["UA"] > 0
    assert abs(d["hot_duty"] - d["cold_duty"]) < 1e-6
    # the hot composite has more than two points (the phase-change kink)
    assert len(d["hot_composite"]) > 2


# -- heat leak ----------------------------------------------------------------
def test_heat_loss_shifts_the_energy_balance():
    pp = _pkg()
    loss = 200.0
    ms = MultiStreamExchanger(
        "HL", {"passes": [{"duty": -2000.0}, {}], "heat_loss": loss})
    ms.solve({"pass1_in": _stream(350.0), "pass2_in": _stream(250.0)}, pp)
    # Sum of pass duties = -heat_loss (the core leaks `loss` W to ambient)
    assert abs(sum(ms.design["pass_duties"]) + loss) < 1e-6


# -- feasibility / error paths ------------------------------------------------
def test_temperature_cross_raises():
    pp = _pkg()
    # cool the hot pass below the cold inlet -> a genuine composite cross
    ms = MultiStreamExchanger("CX", {"passes": [{"duty": -6500.0}, {}]})
    with pytest.raises(MultiStreamExchangerError, match="cross"):
        ms.solve({"pass1_in": _stream(350.0), "pass2_in": _stream(250.0)}, pp)


def test_extreme_duty_gives_typed_error_not_raw_thermo():
    pp = _pkg()
    ms = MultiStreamExchanger("XT", {"passes": [{"duty": -1e5}, {}]})
    with pytest.raises(MultiStreamExchangerError):
        ms.solve({"pass1_in": _stream(350.0), "pass2_in": _stream(250.0)}, pp)


def test_wrong_free_pass_count_raises():
    pp = _pkg()
    ins = {"pass1_in": _stream(380.0), "pass2_in": _stream(250.0),
           "pass3_in": _stream(260.0)}
    # all fixed (zero free) with no global spec -> under-determined
    ms = MultiStreamExchanger(
        "F0", {"passes": [{"T_out": 300}, {"T_out": 310}, {"T_out": 320}]})
    with pytest.raises(MultiStreamExchangerError, match="exactly one pass"):
        ms.solve(ins, pp)
    # all free (three free) with no global spec
    ms = MultiStreamExchanger("F3", {"passes": [{}, {}, {}]})
    with pytest.raises(MultiStreamExchangerError, match="exactly one pass"):
        ms.solve(ins, pp)


def test_too_few_passes_rejected():
    with pytest.raises(MultiStreamExchangerError, match="at least two"):
        MultiStreamExchanger("P1", {"passes": [{}]}).define_ports()


def test_pass_with_both_t_out_and_duty_rejected():
    pp = _pkg()
    ms = MultiStreamExchanger(
        "BD", {"passes": [{"T_out": 300, "duty": -1000}, {}]})
    with pytest.raises(MultiStreamExchangerError, match="both"):
        ms.solve({"pass1_in": _stream(380.0), "pass2_in": _stream(250.0)}, pp)


# -- full flowsheet: solve, round-trip, size, cost ----------------------------
def _lng_fs():
    fs = Flowsheet(components=[Component(c) for c in _COMPS],
                   property_package="thermo:PR")
    fs.add(MultiStreamExchanger(
        "LNG", {"passes": [{"T_out": 260.0, "dP": 20000.0},
                           {"T_out": 320.0}, {}]}))
    fs.feed("H1", "LNG:pass1_in", T=380.0, P=_P, molar_flow=1.0, z=_Z)
    fs.feed("C1", "LNG:pass2_in", T=250.0, P=_P, molar_flow=1.0, z=_Z)
    fs.feed("C2", "LNG:pass3_in", T=260.0, P=_P, molar_flow=1.0, z=_Z)
    for p in ("pass1_out", "pass2_out", "pass3_out"):
        fs.connect(p.upper(), f"LNG:{p}", None)
    return fs


def test_flowsheet_solves_and_applies_pressure_drop():
    fs = _lng_fs()
    assert fs.solve().converged
    d = fs.units["LNG"].design
    assert abs(sum(d["pass_duties"])) < 1e-6
    assert d["min_approach"] > 0
    # the dP on pass 1 is applied
    assert abs(fs.streams["PASS1_OUT"].P - (_P - 20000.0)) < 1e-6


def test_flowsheet_round_trips_and_resolves():
    fs = _lng_fs()
    fs.solve()
    ua = fs.units["LNG"].design["UA"]
    fs2 = from_dict(to_dict(fs))
    fs2.solve()
    assert abs(fs2.units["LNG"].design["UA"] - ua) < 1e-6


def test_sizing_and_costing():
    fs = _lng_fs()
    report = fs.solve()
    sizes = size_flowsheet(fs, report, fs.property_package)
    assert len(sizes) == 1
    s = sizes[0]
    assert s.attribute_name == "area_m2" and s.attribute > 0
    cost = cost_equipment(s)
    assert cost.bare_module > 0
