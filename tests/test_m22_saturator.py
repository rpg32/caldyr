"""M22 tests: the gas Saturator (P9 — the extractor half is covered by the
existing ExtractionColumn, validated in test_m13_specialty_distillation).

The Saturator is the open analogue of HYSYS's Stream Saturator (Hameed 2025
sec. 10.4, humidifying the acid-gas and air feeds of a sulfur-recovery unit). A
gas is loaded with a condensable component (water by default) up to its
saturation partial pressure at the operating T, P.

The headline validation uses **n-hexane**, for which the cubic EOS predicts the
vapour pressure to a few percent (water VLE under PR is poor — the same caveat
the steam-table / humidity work documents — so the water cases here assert the
internal saturation consistency and the mass/energy balances, not an absolute
psychrometric number).
"""
import pytest

from caldyr.core import Component, Flowsheet, Stream
from caldyr.economics.costing import cost_equipment
from caldyr.economics.sizing import size_flowsheet
from caldyr.io import from_dict, to_dict
from caldyr.thermo import make_package
from caldyr.unitops import Saturator, SaturatorError

P_ATM = 101325.0


def _dry_n2_hexane():
    comps = ["nitrogen", "hexane"]
    pp = make_package("thermo:PR", comps)
    gas = Stream(id="g", components=comps, T=300.0, P=P_ATM, molar_flow=1.0,
                 z={"nitrogen": 1.0, "hexane": 0.0})
    return comps, pp, gas


def _hexane_psat(T):
    from thermo import Chemical
    return Chemical("hexane", T=T, P=P_ATM).Psat


# -- saturation thermodynamics ------------------------------------------------
def test_saturates_to_the_vapor_pressure():
    comps, pp, gas = _dry_n2_hexane()
    sat = Saturator("SAT", {"saturant": "hexane"})
    out = sat.solve({"gas_in": gas}, pp)
    g = out["gas_out"]
    # the hexane partial pressure equals its (PR) saturation pressure; PR is
    # within a few % of the literature vapour pressure at 300 K
    p_hexane = g.z["hexane"] * g.P
    assert abs(p_hexane - _hexane_psat(300.0)) / _hexane_psat(300.0) < 0.05
    assert abs(sat.design["relative_humidity_achieved"] - 1.0) < 1e-9
    # the gas leaves as a saturated vapour
    assert g.vapor_fraction == 1.0


def test_relative_humidity_scales_the_loading():
    comps, pp, gas = _dry_n2_hexane()
    full = Saturator("F", {"saturant": "hexane"}).solve({"gas_in": gas}, pp)
    half = Saturator("H", {"saturant": "hexane", "relative_humidity": 0.5}
                     ).solve({"gas_in": gas}, pp)
    # half the relative humidity -> half the saturant partial pressure
    assert abs(half["gas_out"].z["hexane"]
               - 0.5 * full["gas_out"].z["hexane"]) < 1e-6


# -- mass / energy balances ---------------------------------------------------
def test_dry_gas_is_conserved_and_balances_close():
    comps, pp, gas = _dry_n2_hexane()
    sat = Saturator("SAT", {"saturant": "hexane"})
    out = sat.solve({"gas_in": gas}, pp)
    g = out["gas_out"]
    # nitrogen (the dry gas) passes through unchanged
    assert abs(g.molar_flow * g.z["nitrogen"] - 1.0) < 1e-9
    # hexane balance: added saturant = hexane leaving in the gas
    added = sat.design["saturant_added"]
    assert abs(g.molar_flow * g.z["hexane"] - added) < 1e-9
    # isothermal saturation absorbs the latent heat -> positive duty (heat in)
    assert out["duty"].duty > 0.0


def test_wired_supply_matches_auto_supply():
    comps, pp, gas = _dry_n2_hexane()
    auto = Saturator("A", {"saturant": "hexane"})
    a = auto.solve({"gas_in": gas}, pp)
    need = auto.design["saturant_added"]
    # supply exactly the needed saturant as a liquid at the same temperature
    water = Stream(id="w", components=comps, T=300.0, P=P_ATM, molar_flow=need,
                   z={"nitrogen": 0.0, "hexane": 1.0})
    wired = Saturator("W", {"saturant": "hexane"})
    w = wired.solve({"gas_in": gas, "water_in": water}, pp)
    assert abs(w["gas_out"].z["hexane"] - a["gas_out"].z["hexane"]) < 1e-9
    assert abs(wired.design["duty"] - auto.design["duty"]) < 1e-6
    assert w["liquid_out"].molar_flow < 1e-12


def test_excess_supply_drains_to_liquid():
    comps, pp, gas = _dry_n2_hexane()
    base = Saturator("B", {"saturant": "hexane"})
    base.solve({"gas_in": gas}, pp)
    need = base.design["saturant_added"]
    water = Stream(id="w", components=comps, T=300.0, P=P_ATM,
                   molar_flow=need + 0.1, z={"nitrogen": 0.0, "hexane": 1.0})
    sat = Saturator("S", {"saturant": "hexane"})
    out = sat.solve({"gas_in": gas, "water_in": water}, pp)
    assert abs(out["liquid_out"].molar_flow - 0.1) < 1e-6   # surplus drains
    assert abs(sat.design["relative_humidity_achieved"] - 1.0) < 1e-9


def test_short_supply_leaves_gas_sub_saturated():
    comps, pp, gas = _dry_n2_hexane()
    water = Stream(id="w", components=comps, T=300.0, P=P_ATM, molar_flow=0.1,
                   z={"nitrogen": 0.0, "hexane": 1.0})
    sat = Saturator("S", {"saturant": "hexane"})
    out = sat.solve({"gas_in": gas, "water_in": water}, pp)
    assert sat.design["sub_saturated"] is True
    assert sat.design["relative_humidity_achieved"] < 1.0
    assert out["liquid_out"].molar_flow < 1e-12           # all of it evaporates
    # the 0.1 mol supplied all ends up in the gas
    assert abs(out["gas_out"].molar_flow * out["gas_out"].z["hexane"] - 0.1) < 1e-9


# -- the book's use: humidifying an acid-gas / air stream with water ----------
def test_water_saturation_of_air_is_internally_consistent():
    # PR's water vapour pressure is poor (documented), so assert the saturation
    # is self-consistent (the gas's water partial pressure equals PR's own water
    # saturation) and the balances close, not an absolute steam-table number.
    comps = ["nitrogen", "oxygen", "water"]
    pp = make_package("thermo:PR", comps)
    air = Stream(id="air", components=comps, T=320.0, P=P_ATM, molar_flow=10.0,
                 z={"nitrogen": 0.79, "oxygen": 0.21, "water": 0.0})
    sat = Saturator("SAT", {"saturant": "water"})
    out = sat.solve({"gas_in": air}, pp)
    g = out["gas_out"]
    # water partial pressure == the saturation value the unit read for water
    assert abs(g.z["water"] - sat.design["y_saturation"]) < 1e-9
    # dry air conserved, hexane-free
    assert abs(g.molar_flow * g.z["nitrogen"] - 7.9) < 1e-9
    assert abs(g.molar_flow * g.z["oxygen"] - 2.1) < 1e-9


# -- error paths --------------------------------------------------------------
def test_saturant_must_be_a_component():
    comps, pp, gas = _dry_n2_hexane()
    with pytest.raises(SaturatorError, match="not in the component"):
        Saturator("S", {"saturant": "water"}).solve({"gas_in": gas}, pp)


def test_bad_relative_humidity_rejected():
    comps, pp, gas = _dry_n2_hexane()
    with pytest.raises(SaturatorError, match="relative_humidity"):
        Saturator("S", {"saturant": "hexane", "relative_humidity": 1.5}
                  ).solve({"gas_in": gas}, pp)


def test_non_condensable_saturant_rejected():
    # nitrogen cannot condense at 300 K / 1 atm, so it cannot saturate the gas
    comps, pp, gas = _dry_n2_hexane()
    with pytest.raises(SaturatorError, match="does not condense"):
        Saturator("S", {"saturant": "nitrogen"}).solve({"gas_in": gas}, pp)


# -- full flowsheet: solve, round-trip, size, cost ----------------------------
def _saturator_fs():
    comps = ["nitrogen", "hexane"]
    fs = Flowsheet(components=[Component(c) for c in comps],
                   property_package="thermo:PR")
    fs.add(Saturator("SAT", {"saturant": "hexane", "dP": 5000.0}))
    fs.feed("GAS", "SAT:gas_in", T=300.0, P=P_ATM, molar_flow=1.0,
            z={"nitrogen": 1.0, "hexane": 0.0})
    fs.connect("WETGAS", "SAT:gas_out", None)
    fs.connect("DRAIN", "SAT:liquid_out", None)
    fs.connect("Q", "SAT:duty", None)
    return fs


def test_flowsheet_solves_round_trips_and_costs():
    fs = _saturator_fs()
    assert fs.solve().converged
    assert fs.streams["WETGAS"].z["hexane"] > 0.0
    assert abs(fs.streams["WETGAS"].P - (P_ATM - 5000.0)) < 1e-6

    fs2 = from_dict(to_dict(fs))
    report = fs2.solve()
    assert abs(fs2.streams["WETGAS"].z["hexane"]
               - fs.streams["WETGAS"].z["hexane"]) < 1e-9

    pp = make_package(fs2.property_package, fs2.component_ids)
    sizes = size_flowsheet(fs2, report, pp)
    sat = next(s for s in sizes if s.unit_id == "SAT")
    assert sat.attribute > 0
    assert cost_equipment(sat).bare_module > 0
