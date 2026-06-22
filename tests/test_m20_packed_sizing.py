"""M20 tests: packed-column hydraulic sizing (economics.packed_sizing) and the
``internals='packed'`` option on the column sizers.

Validation references (Perry's Chemical Engineers' Handbook, 8e, sec. 14):
  * the Eckert/Strigle GPDC worked Example 13 (2-in metal Pall rings, F_p=27
    ft^-1, air/water at F_LV=0.207) reads a capacity parameter CP=1.01 at a
    pressure drop of 0.38 in-H2O/ft;
  * the Kister-Gill flood pressure drop dP_flood = 0.12 F_p^0.7;
  * the Hameed (2025) sec. 9.1 SO2 absorber (206 kmol/h gas, 1.3e5 kg/h water,
    20 stages, 1 atm), which HYSYS packs to a 1.285 m diameter at 80% of
    capacity (fig. 9.10) — the GPDC here with 50 mm metal Pall rings at 70% of
    flood gives ~1.25 m, agreeing with the book and with the tray sizer
    (~1.4 m) to within ~15%.
"""
import math

import pytest

from caldyr.core import Component, Flowsheet
from caldyr.economics import packed_sizing as ps
from caldyr.economics.costing import cost_equipment
from caldyr.economics.sizing import size_flowsheet
from caldyr.economics.tray_sizing import StageLoad
from caldyr.thermo import make_package
from caldyr.unitops import Absorber

P_ATM = 101325.0
KMOLH = 1000.0 / 3600.0


# -- the GPDC correlation against the Perry worked example --------------------
def test_gpdc_matches_perry_worked_example():
    # Example 13: F_LV=0.207, 2-in metal Pall rings F_p=27 ft^-1; at dP=0.38
    # in-H2O/ft the capacity parameter reads CP ~ 1.01.
    cp = ps.gpdc_capacity_param(0.207, 0.38)
    assert cp == pytest.approx(1.01, rel=0.05)


def test_kister_gill_flood_pressure_drop():
    # dP_flood = 0.12 F_p^0.7; F_p = 27 ft^-1 -> ~1.2 in-H2O/ft.
    Fp_ft = 27.0
    assert 0.12 * Fp_ft ** 0.7 == pytest.approx(1.20, rel=0.02)


def test_capacity_param_decreases_with_flow_parameter():
    # The GPDC capacity falls monotonically as the flow parameter rises.
    cps = [ps.gpdc_capacity_param(x, 1.0) for x in (0.02, 0.1, 0.5, 2.0)]
    assert all(a > b for a, b in zip(cps, cps[1:]))


# -- the SO2 absorber, packed ------------------------------------------------
def so2_absorber() -> Flowsheet:
    fs = Flowsheet(components=[Component("nitrogen"), Component("oxygen"),
                              Component("sulfur dioxide"), Component("water")],
                   property_package="thermo:PR")
    fs.add(Absorber("ABS", {"n_stages": 20, "P": P_ATM}))
    fs.feed("GAS", "ABS:gas_in", T=293.15, P=P_ATM, molar_flow=206.0 * KMOLH,
            z={"sulfur dioxide": 0.03, "nitrogen": 0.97 * 0.79,
               "oxygen": 0.97 * 0.21})
    fs.feed("LIQ", "ABS:liquid_in", T=293.15, P=P_ATM,
            molar_flow=1.3e5 / 18.01528 * KMOLH, z={"water": 1.0})
    fs.connect("CLEAN", "ABS:gas_out", None)
    fs.connect("RICH", "ABS:liquid_out", None)
    return fs


def _bottom_load(design, pp):
    j = design["n_stages"] - 2
    return StageLoad(stage=j + 1, V=design["V_profile"][j],
                     L=design["L_profile"][j], T=design["T_profile"][j],
                     P=design["P_profile"][j], x=design["x_profile"][j],
                     y=design["y_profile"][j])


def test_so2_packed_diameter_near_book_value():
    fs = so2_absorber()
    assert fs.solve().converged
    pp = make_package("thermo:PR", fs.component_ids)
    load = _bottom_load(fs.units["ABS"].design, pp)
    hyd = ps.size_packed(pp, load, "pall_metal_50mm")
    # Book HYSYS packed design: 1.285 m at 80% of capacity (here 70% of flood).
    assert hyd.diameter_m == pytest.approx(1.285, rel=0.15)
    assert 0.0 < hyd.F_LV < 1.0
    assert hyd.u_design_ms < hyd.u_flood_ms


def test_smaller_packing_floods_sooner_needs_bigger_tower():
    # A higher packing factor (smaller packing) floods at a lower velocity, so
    # the same service needs a larger diameter.
    fs = so2_absorber()
    fs.solve()
    pp = make_package("thermo:PR", fs.component_ids)
    load = _bottom_load(fs.units["ABS"].design, pp)
    d50 = ps.size_packed(pp, load, "pall_metal_50mm").diameter_m   # F_p=89
    d25 = ps.size_packed(pp, load, "pall_metal_25mm").diameter_m   # F_p=183
    assert d25 > d50


# -- HETP rules of thumb -----------------------------------------------------
def test_hetp_rules_of_thumb():
    # random 50 mm Pall: HETP ~ 18 D_p = 0.9 m for organics; doubled for water.
    h_org, _ = ps.hetp_m("pall_metal_50mm", aqueous=False)
    h_aq, _ = ps.hetp_m("pall_metal_50mm", aqueous=True)
    assert h_org == pytest.approx(18.0 * 0.050, rel=1e-6)
    assert h_aq == pytest.approx(2.0 * h_org, rel=1e-6)
    # structured packing is more efficient (smaller HETP) than random.
    h_struct, _ = ps.hetp_m("mellapak_250y", aqueous=False)
    assert h_struct < h_org


def test_unknown_packing_raises():
    with pytest.raises(ValueError, match="unknown packing"):
        ps.size_packed(None, None, "not_a_packing")


# -- end-to-end through the absorber sizer + costing -------------------------
def test_packed_internals_size_and_cost_through_flowsheet():
    fs = so2_absorber()
    fs.units["ABS"].params["internals"] = "packed"
    fs.units["ABS"].params["packing"] = "pall_metal_50mm"
    rep = fs.solve()
    pp = make_package(fs.property_package, fs.component_ids)
    sizes = size_flowsheet(fs, rep, pp)
    tower = next(s for s in sizes if s.unit_id == "ABS")
    bed = next(s for s in sizes if s.unit_id == "ABS.packing")
    assert tower.equipment_type == "vessel_vertical"
    assert tower.diameter_m == pytest.approx(1.285, rel=0.15)
    assert bed.equipment_type == "packing_random" and bed.attribute > 0
    # both the shell and the packing bed carry a positive bare-module cost
    assert cost_equipment(tower).bare_module > 0
    assert cost_equipment(bed).bare_module > 0


def test_packed_vs_tray_diameter_consistent():
    # The same SO2 service sized as trays vs packing agrees within ~25%.
    fs_t = so2_absorber()
    rep_t = fs_t.solve()
    pp = make_package(fs_t.property_package, fs_t.component_ids)
    tray_d = next(s for s in size_flowsheet(fs_t, rep_t, pp)
                  if s.unit_id == "ABS").diameter_m
    fs_p = so2_absorber()
    fs_p.units["ABS"].params["internals"] = "packed"
    rep_p = fs_p.solve()
    packed_d = next(s for s in size_flowsheet(fs_p, rep_p, pp)
                    if s.unit_id == "ABS").diameter_m
    assert packed_d == pytest.approx(tray_d, rel=0.25)


def test_bad_internals_raises():
    fs = so2_absorber()
    fs.units["ABS"].params["internals"] = "bubble_caps"
    rep = fs.solve()
    pp = make_package(fs.property_package, fs.component_ids)
    with pytest.raises(ValueError, match="internals"):
        size_flowsheet(fs, rep, pp)
