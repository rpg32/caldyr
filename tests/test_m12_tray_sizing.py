"""M12 tests: tray hydraulic sizing (Fair's flooding correlation).

The column diameter now comes from Fair's (1961) flooding capacity chart in
the Lygeros & Magoulas closed form quoted by Seader, Henley & Roper,
*Separation Process Principles* 3e, sec. 6.6 (eq. 6-42), at 80% of flooding
on the net (downcomer-corrected) area, evaluated at the governing tray of
the converged stage profiles (top-most and bottom-most trays).

Validation:
  * the correlation itself reproduces Fair's chart values (C_sb ~ 0.09-0.11
    m/s at F_LV = 0.1 and 24-in spacing — Seader 3e fig. 6.24, Wankat 3e
    fig. 10-25) and its monotonic structure;
  * the book absorber of Hameed 2025 sec. 9.1 (206 kmol/h gas, 1.3e5 kg/h
    water, 20 stages at 1 atm) sizes to D ~ 1.4 m vs the book's HYSYS
    packed-column design of **1.285 m** at 80% of max capacity (fig. 9.10)
    — tray vs packed diameters for the same service within ~10%;
  * the replacement stays within sanity of the old fixed-F-factor heuristic
    (1.2 Pa^0.5) on the benzene/toluene benchmark: 2.0 m vs 2.5 m (the old
    rule was conservative; 80% of Fair flooding corresponds to F_s ~ 2.2
    Pa^0.5 here, squarely in the normal sieve-tray range).
"""
import math

import pytest

from caldyr.core import Component, Flowsheet
from caldyr.economics import tray_sizing
from caldyr.economics.sizing import SizingOptions, size_flowsheet
from caldyr.thermo import make_package
from caldyr.unitops import Absorber, RigorousColumn

P_ATM = 101325.0
KMOLH = 1000.0 / 3600.0


# -- 1. the correlation itself -------------------------------------------------
def test_fair_capacity_matches_chart_anchor():
    """Fair's chart at F_LV = 0.1 and 24-in (0.61 m) spacing reads
    C_sb,f ~ 0.30-0.36 ft/s = 0.09-0.11 m/s (Seader 3e fig. 6.24; Wankat 3e
    fig. 10-25). The Lygeros & Magoulas fit must land in that window."""
    c = tray_sizing.fair_capacity(0.1, 0.61)
    assert 0.08 < c < 0.12


def test_fair_capacity_monotone():
    # More liquid load (higher F_LV) -> lower capacity.
    cs = [tray_sizing.fair_capacity(f, 0.61) for f in (0.01, 0.1, 0.5, 1.0)]
    assert all(a > b for a, b in zip(cs, cs[1:]))
    # Wider tray spacing -> higher capacity.
    spacings = [tray_sizing.fair_capacity(0.1, s) for s in (0.3, 0.45, 0.61)]
    assert all(a < b for a, b in zip(spacings, spacings[1:]))


def test_size_tray_responds_to_load():
    """Doubling the vapor load must grow the diameter by ~sqrt(2)."""
    pp = make_package("thermo:PR", ["benzene", "toluene"])
    x = {"benzene": 0.5, "toluene": 0.5}
    res = pp.bubble_point(P_ATM, x)
    assert res.y is not None
    base = tray_sizing.size_tray(pp, tray_sizing.StageLoad(
        stage=1, V=50.0, L=40.0, T=res.T, P=P_ATM, x=x, y=res.y))
    double = tray_sizing.size_tray(pp, tray_sizing.StageLoad(
        stage=1, V=100.0, L=80.0, T=res.T, P=P_ATM, x=x, y=res.y))
    assert double.diameter_m == pytest.approx(
        base.diameter_m * math.sqrt(2.0), rel=1e-6)
    assert 0.0 < base.u_design_ms == base.flood_fraction * base.u_flood_ms
    assert base.rho_L > base.rho_V > 0.0


# -- 2. the book absorber's tray design (Hameed 2025 sec. 9.1.3) -----------------
def test_book_so2_absorber_diameter_near_book_packed_design():
    """Book fig. 9.10: HYSYS packs the 20-stage SO2 absorber into a 1.285 m
    column at 80% of maximum capacity. Sizing the same service as a sieve
    tray column at 80% of Fair flooding must land in the same hydraulic
    ballpark (achieved: ~1.40 m, +9%; asserted within 25% — tray and packed
    capacities for one service genuinely differ by that order)."""
    fs = Flowsheet(components=[Component("water"), Component("nitrogen"),
                               Component("oxygen"),
                               Component("sulfur dioxide")],
                   property_package="thermo:PR")
    fs.add(Absorber("ABS", {"n_stages": 20, "P": P_ATM}))
    fs.feed("GAS", "ABS:gas_in", T=293.15, P=P_ATM, molar_flow=206.0 * KMOLH,
            z={"sulfur dioxide": 0.03, "nitrogen": 0.97 * 0.79,
               "oxygen": 0.97 * 0.21})
    fs.feed("WATER", "ABS:liquid_in", T=293.15, P=P_ATM,
            molar_flow=1.3e5 / 18.01528 * KMOLH, z={"water": 1.0})
    fs.connect("GASOUT", "ABS:vapor_out", None)
    fs.connect("LIQOUT", "ABS:liquid_out", None)
    rep = fs.solve()
    assert rep.converged

    pp = make_package("thermo:PR", fs.component_ids)
    sizes = size_flowsheet(fs, rep, pp, SizingOptions())
    tower = next(s for s in sizes if s.unit_id == "ABS")
    assert tower.diameter_m == pytest.approx(1.285, rel=0.25)
    # The high liquid load of a water wash pushes F_LV toward the top of the
    # Fair chart — recorded in the notes for auditability.
    assert any("Fair flooding" in n for n in tower.notes)


# -- 3. the replacement stays within sanity of the old heuristic ------------------
def bt_column() -> Flowsheet:
    fs = Flowsheet(components=[Component("benzene"), Component("toluene")],
                   property_package="thermo:PR")
    fs.add(RigorousColumn("COL", {"n_stages": 16, "feed_stage": 8,
                                  "reflux_ratio": 1.545,
                                  "distillate_rate": 50.0, "P": P_ATM}))
    fs.feed("FEED", "COL:in1", T=365.0, P=P_ATM, molar_flow=100.0,
            z={"benzene": 0.5, "toluene": 0.5})
    for sid, port in (("DIST", "distillate"), ("BOT", "bottoms"),
                      ("QC", "condenser_duty"), ("QR", "reboiler_duty")):
        fs.connect(sid, f"COL:{port}", None)
    return fs


def test_distillation_diameter_within_sanity_of_old_heuristic():
    """The old sizer used a fixed F-factor of 1.2 Pa^0.5 on the overhead
    vapor (D = 2.53 m on this benchmark). Fair flooding at 80% gives 2.01 m
    (-20%): same order, the old rule was conservative. Assert the new result
    sits within the +-40% sanity band of the old one and that the governing
    tray is the bottom one (higher vapor density *and* traffic here)."""
    fs = bt_column()
    rep = fs.solve()
    pp = make_package("thermo:PR", ["benzene", "toluene"])
    sizes = size_flowsheet(fs, rep, pp, SizingOptions())
    tower = next(s for s in sizes if s.unit_id == "COL")
    d_old = 2.53
    assert 0.6 * d_old < tower.diameter_m < 1.4 * d_old
    assert tower.diameter_m == pytest.approx(2.01, abs=0.05)

    # Same diameter computed straight from the profiles: the governing tray
    # of {top tray, bottom tray} is the bottom one.
    design = fs.units["COL"].design
    n = design["n_stages"]
    hyd = tray_sizing.governing_tray(
        pp, tray_sizing.loads_from_profiles(design, [1, n - 2]))
    assert hyd.stage == n - 1                 # 1-based bottom tray
    assert hyd.diameter_m == pytest.approx(tower.diameter_m, rel=1e-9)


def test_flood_fraction_option_scales_the_tower():
    """Halving the design flood fraction must grow the area by 2x (diameter
    by sqrt(2)) — the option is honored end to end."""
    fs = bt_column()
    rep = fs.solve()
    pp = make_package("thermo:PR", ["benzene", "toluene"])
    d80 = next(s for s in size_flowsheet(fs, rep, pp, SizingOptions())
               if s.unit_id == "COL").diameter_m
    d40 = next(s for s in size_flowsheet(
        fs, rep, pp, SizingOptions(tray_flood_fraction=0.4))
        if s.unit_id == "COL").diameter_m
    assert d40 == pytest.approx(d80 * math.sqrt(2.0), rel=1e-9)


# -- 5. tray internals rating (P10): pressure drop, weeping, downcomer backup --
def _bt_load(V=100.0, L=80.0):
    return tray_sizing.StageLoad(
        stage=10, V=V, L=L, T=370.0, P=1.2e5,
        x={"benzene": 0.5, "toluene": 0.5},
        y={"benzene": 0.7, "toluene": 0.3})


def test_tray_rating_pressure_drop_is_physical():
    """A sieve tray's total pressure drop is a few tenths of a kPa to ~1 kPa
    (Seader 3e sec. 6.6; Kister ch. 6): dry head + clear-liquid head, each tens
    of mm of liquid."""
    pp = make_package("thermo:PR", ["benzene", "toluene"])
    hyd = tray_sizing.size_tray(pp, _bt_load(), flood_fraction=0.8)
    r = hyd.rating
    assert r is not None
    assert 0.02 < r.head_dry_m < 0.15           # dry head 20-150 mm
    assert 0.0 < r.weir_crest_m < 0.05          # weir crest a few cm
    assert r.head_total_m == pytest.approx(r.head_dry_m + r.head_liquid_m)
    assert 300.0 < r.delta_P_Pa < 1500.0        # ~0.3-1.5 kPa per tray
    # delta_P is the static head of the total froth height
    assert r.delta_P_Pa == pytest.approx(hyd.rho_L * 9.80665 * r.head_total_m)


def test_tray_rating_matches_the_hand_correlations():
    """Recompute h_d and the Francis weir crest from first principles and match
    the rating — the formulas are the citation."""
    pp = make_package("thermo:PR", ["benzene", "toluene"])
    load = _bt_load()
    hyd = tray_sizing.size_tray(pp, load, flood_fraction=0.8,
                                downcomer_frac=0.10)
    r = hyd.rating
    a_active = hyd.area_m2 * (1.0 - 2.0 * 0.10)
    a_hole = 0.10 * a_active
    q_v = load.V * pp.volume_vapor(load.T, load.P, load.y)
    q_l = load.L * pp.volume_liquid(load.T, load.P, load.x)
    u_h = q_v / a_hole
    h_d = u_h ** 2 / (2.0 * 9.80665) * (hyd.rho_V / hyd.rho_L) / 0.73 ** 2
    L_w = 0.73 * hyd.diameter_m * r.n_passes
    how = 0.664 * (q_l / L_w) ** (2.0 / 3.0)
    assert r.head_dry_m == pytest.approx(h_d, rel=1e-9)
    assert r.weir_crest_m == pytest.approx(how, rel=1e-9)


def test_low_vapor_rate_weeps():
    """At a low vapour rate the dynamic + surface-tension head can no longer
    hold the liquid on the tray (h_d + h_sigma < h_w + how) -> weeping."""
    pp = make_package("thermo:PR", ["benzene", "toluene"])
    high = tray_sizing.size_tray(pp, _bt_load(V=120.0, L=60.0), flood_fraction=0.85)
    low = tray_sizing.size_tray(pp, _bt_load(V=12.0, L=60.0), flood_fraction=0.85)
    assert not high.rating.weeps                 # comfortably above the weep point
    assert low.rating.weeps                       # starved of vapour -> weeps
    assert low.rating.weep_margin_m < 0.0


def test_n_passes_increases_with_diameter():
    assert tray_sizing._n_passes(1.5) == 1
    assert tray_sizing._n_passes(3.0) == 2
    assert tray_sizing._n_passes(4.5) == 3
    assert tray_sizing._n_passes(6.0) == 4


def test_column_sizer_surfaces_the_internals_note():
    """The internals rating rides through to the equipment-size notes so the
    tray pressure drop and operability flags are visible after costing."""
    fs = bt_column()
    rep = fs.solve()
    pp = make_package("thermo:PR", ["benzene", "toluene"])
    tower = next(s for s in size_flowsheet(fs, rep, pp, SizingOptions())
                 if s.unit_id == "COL")
    joined = " ".join(tower.notes)
    assert "internals:" in joined and "tray dP" in joined
