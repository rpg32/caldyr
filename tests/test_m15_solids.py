"""M15 tests: solids operations — Cyclone, RotaryVacuumFilter, BaghouseFilter
(Hameed, *Chemical Process Simulations using Aspen HYSYS*, Wiley 2025, ch. 12).

The v1 particle model (see caldyr/unitops/solids.py): the solid is an ordinary
stream component; the particle size distribution is a *unit* param (``psd``),
exactly how the book's worked example enters it. PSDs do not propagate
between units.

References & measured deltas
----------------------------
* **Hameed sec. 12.1 (cyclone)** — 150 kmol/h of air at 45 C / 150 kPa with
  5 wt% carbon dust (HYSYS particle density 1642 kg/m^3); the HYSYS design
  mode hits exactly the 95% efficiency spec (Leith/Licht) with:
  High Output (Stairmand HT ratios, Fig. 12.6): D = 2.938 m, 2 cyclones
  (single 1 mm particle); D = 2.881 m, 3 cyclones with the discrete PSD
  {1.5 mm: 40%, 2.0: 35%, 2.5: 15%, 3.0: 10%} (Fig. 12.9); High Efficiency
  (Stairmand HE, Fig. 12.5): D = 2.985 m, 34 cyclones.
  This engine's **Lapple** model evaluated at the book's geometries gives
  overall efficiencies of 96.4% / 98.4% / 92.1% vs the book's 95% design
  point — deltas of +1.4 / +3.4 / -2.9 percentage points (different grade
  model: Lapple vs HYSYS Leith/Licht, which the book itself notes predicts
  *higher* efficiency than Lapple).
* **Cooper & Alley, Air Pollution Control: A Design Approach, 4e, ch. 4** —
  Lapple cut diameter (Eq. 4.7), grade efficiency (Eq. 4.9), effective turns
  (Eq. 4.4), Shepherd-Lapple dP = 0.5 rho v^2 N_H with N_H = 16HW/De^2
  (Eqs. 4.10-4.11), standard geometries (Table 4.2). The classic hand case —
  standard Lapple proportions, D = 1 m, air at ~350 K / 1 atm, v_i = 15 m/s,
  rho_p = 1600 kg/m^3 — gives d50 = 7.2 um and dP ~ 0.9 kPa; reproduced
  below to 0.1% with an independent in-test hand calculation.
* **Hameed sec. 12.2 (rotary vacuum filter)** — 4000 kg/h slurry, 30 wt%
  calcium in water at 25 C / 160 kPa, dP = 10 kPa, 0.8 m drum: HYSYS reports
  area 2.4368 m^2 and width 0.96957 m (Fig. 12.13). HYSYS computes the cake
  resistance internally from particle size/sphericity and the book does not
  publish it, so alpha is *backed out* (1.032e14 m/kg) — the test then
  reproduces the book's area AND the drum geometry chain W = A/(2 pi R)
  exactly. The filtration physics itself is checked against McCabe, Smith &
  Harriott, *Unit Operations*, 7e, ch. 29 (continuous rotary filtration,
  negligible medium resistance) on the data of the book's Exercise 12.2
  (= the McCabe CaCO3 worked case): independent hand recomputation matches
  the unit to 1e-9; the absolute area differs ~ +27% from a real-water-
  density hand value because PR overpredicts liquid water molar volume
  (known cubic-EOS bias, documented below).
* **Hameed sec. 12.3 (baghouse)** — 60 kmol/h air + 5 mol% sulfur at 35 C /
  180 kPa. Air-to-cloth sizing A = Q/v_f and the filter-drag dP model
  dP = S_E v + K2 c v^2 t are from Cooper & Alley 4e ch. 6 / EPA Air
  Pollution Control Cost Manual 6e sec. 6 ch. 1. The book's HYSYS filtration
  time at a 2 kPa dirty-bag dP is 16436:54:34 = 5.92e7 s; HYSYS's internal
  S_E/K2/face-velocity defaults are not published, so only the *structure*
  (time linear in dP, the book's Fig. 12.18) is asserted, not the number —
  with our mid-range defaults the same state gives ~2.9e3 s.
* **Hameed Exercise 12.3 / sec. 12.3.2 steps 1-14 (cyclone + baghouse
  train)** — 10,000 kg/h air + 2000 kg/h kaolin at 30 C / 5 bar, kaolin PSD
  log-probability (mean 8.343 um, sigma 5.227) discretized by HYSYS into the
  9 bins of Fig. 12.22; a 75%-efficient cyclone leaves 500 kg/h in the gas
  (book p. 409) and the baghouse then removes ~100%. Kaolin has no critical
  constants in the databank, so carbon stands in as the dust component (same
  mass rates via its own MW; the PSD and particle density are the book's
  kaolin values) — with the cyclone body sized to 75% the train reproduces
  the book's 500 kg/h to <0.1% and emits 0.5 kg/h (book: 0.0; our default
  baghouse efficiency is 99.9%, not 100%).
"""
import math

import pytest

from caldyr.core import Component, Flowsheet
from caldyr.economics.costing import cost_equipment
from caldyr.economics.sizing import size_flowsheet
from caldyr.io import from_dict, to_dict
from caldyr.thermo import make_package
from caldyr.unitops import BaghouseFilter, Cyclone, RotaryVacuumFilter, SolidsOperationError
from caldyr.unitops.solids import (
    CYCLONE_GEOMETRIES,
    lapple_d50,
    lapple_grade_efficiency,
)

MW_C = 0.0120107        # kg/mol carbon
MW_CA = 0.040078        # kg/mol calcium
MW_W = 0.01801528       # kg/mol water
MW_S = 0.032065         # kg/mol sulfur (databank value)
MW_AIR = 0.028851       # kg/mol 79/21 N2/O2

# 5 wt% carbon in 79/21 air (the book's sec. 12.1 feed), as mole fractions.
_mc, _mair = 0.05 / MW_C, 0.95 / MW_AIR
X_C = _mc / (_mc + _mair)
Z_AIR_C = {"nitrogen": (1 - X_C) * 0.79, "oxygen": (1 - X_C) * 0.21, "carbon": X_C}


def solid_mass_rate(stream, comp: str, mw: float) -> float:
    """kg/s of one component in a solved stream."""
    return stream.molar_flow * stream.z.get(comp, 0.0) * mw


def cyclone_fs(params: dict, *, T=318.15, P=150e3, n=150e3 / 3600.0,
               z=None) -> Flowsheet:
    fs = Flowsheet(
        components=[Component("nitrogen"), Component("oxygen"), Component("carbon")],
        property_package="thermo:PR")
    fs.add(Cyclone("CYC", params))
    fs.feed("F", "CYC:gas_in", T=T, P=P, molar_flow=n, z=z or dict(Z_AIR_C))
    fs.connect("G", "CYC:gas_out", None)
    fs.connect("S", "CYC:solids_out", None)
    return fs


PSD_1MM = [{"d_microns": 1000.0, "mass_frac": 1.0}]
PSD_BOOK = [{"d_microns": 1500.0, "mass_frac": 0.40},
            {"d_microns": 2000.0, "mass_frac": 0.35},
            {"d_microns": 2500.0, "mass_frac": 0.15},
            {"d_microns": 3000.0, "mass_frac": 0.10}]


# -- cyclone: book sec. 12.1 reproductions -------------------------------------
def test_book_cyclone_high_output_single_particle():
    """Book Fig. 12.6: High Output (Stairmand HT), D = 2.938 m, 2 cyclones,
    single 1 mm carbon particle, designed by HYSYS to exactly 95% efficiency.
    Our Lapple grade model at that geometry: 96.4% (delta +1.4 points)."""
    fs = cyclone_fs({"solids": "carbon", "particle_density": 1642.0,
                     "geometry": "Stairmand_HT", "body_diameter": 2.938,
                     "n_cyclones": 2, "psd": PSD_1MM})
    fs.solve()
    d = fs.units["CYC"].design
    assert d["overall_efficiency"] >= 0.95                  # the book's spec
    assert d["overall_efficiency"] == pytest.approx(0.9635, abs=2e-3)
    # exact component mass balance around the split
    feed, gas, sol = fs.streams["F"], fs.streams["G"], fs.streams["S"]
    for comp in ("nitrogen", "oxygen", "carbon"):
        n_in = feed.molar_flow * feed.z[comp]
        n_out = (gas.molar_flow * gas.z.get(comp, 0.0)
                 + sol.molar_flow * sol.z.get(comp, 0.0))
        assert n_out == pytest.approx(n_in, rel=1e-12)
    # all the carrier gas leaves with the gas; the solids stream is pure dust
    assert sol.z.get("nitrogen", 0.0) == 0.0
    assert sol.z.get("carbon") == pytest.approx(1.0)
    # the split is the efficiency
    cap = solid_mass_rate(sol, "carbon", MW_C)
    assert cap == pytest.approx(d["solids_captured_kg_s"], rel=1e-12)
    assert cap / d["solids_in_kg_s"] == pytest.approx(d["overall_efficiency"], rel=1e-12)


def test_book_cyclone_high_output_psd():
    """Book Fig. 12.9: the same problem with the discrete PSD (min 1 mm) —
    HYSYS picks D = 2.881 m, 3 cyclones for 95%. Lapple there: 98.4%
    (delta +3.4 points; every bin is far above the ~230 um cut size)."""
    fs = cyclone_fs({"solids": "carbon", "particle_density": 1642.0,
                     "geometry": "Stairmand_HT", "body_diameter": 2.881,
                     "n_cyclones": 3, "psd": PSD_BOOK})
    fs.solve()
    d = fs.units["CYC"].design
    assert d["overall_efficiency"] >= 0.95
    assert d["overall_efficiency"] == pytest.approx(0.9842, abs=2e-3)
    # grade efficiency is monotone in particle size (psd bins are sorted)
    effs = [g["efficiency"] for g in d["grade"]]
    assert effs == sorted(effs)
    assert all(0.0 < e < 1.0 for e in effs)


def test_book_cyclone_high_efficiency_geometry():
    """Book Fig. 12.5: High Efficiency (Stairmand HE), D = 2.985 m, 34
    cyclones for the 95% spec. Lapple at that geometry: 92.1% (delta -2.9
    points) — the book itself notes Lapple 'greatly underestimates' vs
    Leith/Licht (sec. 12.1.2), so a small shortfall is the expected sign."""
    fs = cyclone_fs({"solids": "carbon", "particle_density": 1642.0,
                     "geometry": "Stairmand_HE", "body_diameter": 2.985,
                     "n_cyclones": 34, "psd": PSD_1MM})
    fs.solve()
    d = fs.units["CYC"].design
    assert d["overall_efficiency"] == pytest.approx(0.9214, abs=2e-3)
    assert d["overall_efficiency"] == pytest.approx(0.95, abs=0.04)


# -- cyclone: independent textbook hand calculation ----------------------------
def test_hand_lapple_cut_diameter_cooper_alley():
    """Classic Lapple design point (Cooper & Alley 4e ch. 4): standard Lapple
    proportions, D = 1.0 m, air at 350 K / 1 atm, v_i = 15 m/s, rho_p = 1600
    kg/m^3 -> d50 = 7.2 um, dP = N_H velocity heads ~ 0.90 kPa. The hand
    numbers below are recomputed independently from the cited formulas."""
    pp = make_package("thermo:PR", ["nitrogen", "oxygen", "carbon"])
    vm = pp.volume(350.0, 101325.0, {"nitrogen": 0.79, "oxygen": 0.21})
    geo = CYCLONE_GEOMETRIES["Lapple"]
    n_gas = 15.0 * geo.H * geo.W / vm           # Q = v_i * H * W for D = 1 m
    xc = 1e-3
    psd = [{"d_microns": dm, "mass_frac": w} for dm, w in
           [(1, 0.03), (5, 0.20), (10, 0.15), (20, 0.20),
            (30, 0.16), (50, 0.16), (100, 0.10)]]
    fs = cyclone_fs({"solids": "carbon", "particle_density": 1600.0,
                     "geometry": "Lapple", "body_diameter": 1.0, "psd": psd},
                    T=350.0, P=101325.0, n=n_gas / (1 - xc),
                    z={"nitrogen": 0.79 * (1 - xc), "oxygen": 0.21 * (1 - xc),
                       "carbon": xc})
    fs.solve()
    d = fs.units["CYC"].design

    # Lapple standard proportions: N_e = (2 + 2/2)/0.5 = 6 turns, N_H = 8.
    assert d["Ne_turns"] == pytest.approx(6.0)
    assert d["NH_velocity_heads"] == pytest.approx(8.0)
    assert d["inlet_velocity_m_s"] == pytest.approx(15.0, rel=1e-9)
    # transport properties vs handbook air values (Perry's 8e: mu ~ 2.08e-5
    # Pa s at 350 K; ideal-gas density P MW / R T = 1.005 kg/m^3)
    assert d["gas_viscosity_Pa_s"] == pytest.approx(2.08e-5, rel=2e-2)
    rho_ideal = 101325.0 * MW_AIR / (8.314462618 * 350.0)
    assert d["gas_density_kg_m3"] == pytest.approx(rho_ideal, rel=5e-3)

    # the cut diameter, recomputed by hand from Cooper & Alley Eq. 4.7
    mu, rho = d["gas_viscosity_Pa_s"], d["gas_density_kg_m3"]
    d50_hand = math.sqrt(9.0 * mu * 0.25 / (2.0 * math.pi * 6.0 * 15.0 * (1600.0 - rho)))
    assert d["d50_m"] == pytest.approx(d50_hand, rel=1e-9)
    assert d["d50_microns"] == pytest.approx(7.2, abs=0.15)   # the textbook answer

    # Shepherd-Lapple pressure drop, by hand: 0.5 rho v^2 N_H ~ 904 Pa
    assert d["dP_Pa"] == pytest.approx(0.5 * rho * 15.0 ** 2 * 8.0, rel=1e-9)
    assert d["dP_Pa"] == pytest.approx(904.0, rel=1e-2)
    assert fs.streams["G"].P == pytest.approx(101325.0 - d["dP_Pa"], rel=1e-12)

    # the overall efficiency is the mass-weighted grade sum, by hand
    eta_hand = sum(b["mass_frac"]
                   * (1.0 / (1.0 + (d50_hand / (b["d_microns"] * 1e-6)) ** 2))
                   for b in psd)
    assert d["overall_efficiency"] == pytest.approx(eta_hand, rel=1e-9)


def test_lapple_grade_efficiency_properties():
    # exactly 50% at the cut size; monotone in d; -> 0 small, -> 1 large
    assert lapple_grade_efficiency(5e-6, 5e-6) == pytest.approx(0.5)
    es = [lapple_grade_efficiency(5e-6, d * 1e-6) for d in (0.5, 1, 2, 5, 10, 50, 500)]
    assert es == sorted(es)
    assert es[0] < 0.02 and es[-1] > 0.999
    # d50 shrinks (better capture) with faster inlet velocity and denser dust
    base = lapple_d50(1.8e-5, 0.25, 6.0, 15.0, 2000.0, 1.2)
    assert lapple_d50(1.8e-5, 0.25, 6.0, 30.0, 2000.0, 1.2) < base
    assert lapple_d50(1.8e-5, 0.25, 6.0, 15.0, 4000.0, 1.2) < base


def test_cyclone_velocity_sizing_mode():
    """Sizing from a design inlet velocity: D = sqrt(Q_each / (v H_r W_r))."""
    fs = cyclone_fs({"solids": "carbon", "particle_density": 1642.0,
                     "inlet_velocity": 18.0, "n_cyclones": 2, "psd": PSD_1MM})
    fs.solve()
    d = fs.units["CYC"].design
    assert d["inlet_velocity_m_s"] == pytest.approx(18.0, rel=1e-9)
    geo = CYCLONE_GEOMETRIES["Lapple"]
    d_hand = math.sqrt(d["Q_per_cyclone_m3_s"] / (18.0 * geo.H * geo.W))
    assert d["body_diameter_m"] == pytest.approx(d_hand, rel=1e-9)


def test_cyclone_typed_errors():
    ok = {"solids": "carbon", "particle_density": 1642.0,
          "body_diameter": 1.0, "psd": PSD_1MM}
    cases = [
        (dict(ok, solids=None), "solids"),
        (dict(ok, solids="quartz"), "not in the flowsheet"),
        (dict(ok, psd=None), "psd"),
        (dict(ok, psd=[]), "psd"),
        (dict(ok, psd=[{"d_microns": -5.0, "mass_frac": 1.0}]), "non-positive"),
        (dict(ok, psd=[{"d_microns": 5.0, "mass_frac": -0.1}]), "negative"),
        (dict(ok, psd=[{"d_microns": 5.0, "mass_frac": 0.5}]), "sum"),
        (dict(ok, psd=[{"diameter": 5.0}]), "psd bin"),
        (dict(ok, particle_density=None), "particle_density"),
        (dict(ok, particle_density=0.5), "exceed the gas density"),
        (dict(ok, geometry="Swift"), "unknown geometry"),
        (dict(ok, n_cyclones=0), "n_cyclones"),
        ({"solids": "carbon", "particle_density": 1642.0, "psd": PSD_1MM},
         "body_diameter"),
    ]
    for params, match in cases:
        with pytest.raises(SolidsOperationError, match=match):
            cyclone_fs(params).solve()
    assert issubclass(SolidsOperationError, ValueError)     # engine convention


def test_cyclone_pressure_drop_exhausting_inlet_raises():
    # 200 m/s inlet velocity at 1 atm: dP = 0.5 rho v^2 N_H ~ 1.6e5 Pa > P.
    fs = cyclone_fs({"solids": "carbon", "particle_density": 1642.0,
                     "inlet_velocity": 200.0, "psd": PSD_1MM}, P=101325.0)
    with pytest.raises(SolidsOperationError, match="pressure drop"):
        fs.solve()


def test_cyclone_all_solid_feed_raises():
    fs = cyclone_fs({"solids": "carbon", "particle_density": 1642.0,
                     "body_diameter": 1.0, "psd": PSD_1MM},
                    z={"nitrogen": 0.0, "oxygen": 0.0, "carbon": 1.0})
    with pytest.raises(SolidsOperationError, match="no gas"):
        fs.solve()


def test_cyclone_flow_roundtrip_and_costing():
    fs = cyclone_fs({"solids": "carbon", "particle_density": 1642.0,
                     "geometry": "Stairmand_HT", "body_diameter": 2.938,
                     "n_cyclones": 2, "psd": PSD_BOOK})
    rep = fs.solve()
    doc = to_dict(fs)
    assert doc == to_dict(from_dict(doc))
    fs2 = from_dict(doc)
    fs2.solve()
    assert (fs2.units["CYC"].design["overall_efficiency"]
            == pytest.approx(fs.units["CYC"].design["overall_efficiency"], rel=1e-12))

    pp = make_package(fs.property_package, fs.component_ids)
    sizes = size_flowsheet(fs, rep, pp)
    assert len(sizes) == 1
    size = sizes[0]
    assert size.equipment_type == "cyclone"
    assert size.quantity == 2                      # parallel cyclones
    assert size.attribute == pytest.approx(
        fs.units["CYC"].design["Q_per_cyclone_m3_s"], rel=1e-12)
    cost = cost_equipment(size)
    assert math.isfinite(cost.bare_module) and cost.bare_module > 0.0
    assert cost.bare_module == pytest.approx(2.0 * cost.purchased, rel=1e-12)  # Fbm=2


# -- rotary vacuum filter ------------------------------------------------------
def rvf_fs(params: dict, *, T=298.15, P=160e3, m_solids=1200.0, m_water=2800.0,
           solid="calcium", mw_solid=MW_CA) -> Flowsheet:
    """A solids-in-water slurry feed (mass rates in kg/h)."""
    n_s = (m_solids / 3600.0) / mw_solid
    n_w = (m_water / 3600.0) / MW_W
    n = n_s + n_w
    fs = Flowsheet(components=[Component("water"), Component(solid)],
                   property_package="thermo:PR")
    fs.add(RotaryVacuumFilter("RVF", params))
    fs.feed("F", "RVF:slurry_in", T=T, P=P, molar_flow=n,
            z={"water": n_w / n, solid: n_s / n})
    fs.connect("L", "RVF:filtrate_out", None)
    fs.connect("C", "RVF:cake_out", None)
    return fs


BOOK_RVF = {"solids": "calcium", "pressure_drop": 10e3, "cycle_time_s": 600.0,
            "submergence": 0.20, "cake_moisture": 0.5, "drum_radius_m": 0.4,
            # HYSYS computes alpha internally from particle size/sphericity and
            # the book does not publish it; this value is backed out so the
            # area matches the book's 2.4368 m^2 (see module docstring).
            "alpha": 1.032e14}


def test_book_rvf_area_and_drum_width():
    """Book Fig. 12.13: radius 0.4 m -> area 2.4368 m^2, width 0.96957 m.
    With the backed-out alpha the area matches; the drum geometry chain
    A = 2 pi R W is then an independent reproduction of the book's width."""
    fs = rvf_fs(BOOK_RVF)
    fs.solve()
    d = fs.units["RVF"].design
    assert d["area_m2"] == pytest.approx(2.4368, rel=2e-3)
    assert d["drum_width_m"] == pytest.approx(0.96957, rel=2e-3)
    assert d["area_m2"] == pytest.approx(
        2.0 * math.pi * 0.4 * d["drum_width_m"], rel=1e-12)


def test_book_rvf_mass_balance_and_cake_moisture():
    fs = rvf_fs(BOOK_RVF)
    fs.solve()
    feed, filt, cake = fs.streams["F"], fs.streams["L"], fs.streams["C"]
    for comp in ("water", "calcium"):
        n_in = feed.molar_flow * feed.z[comp]
        n_out = (filt.molar_flow * filt.z.get(comp, 0.0)
                 + cake.molar_flow * cake.z.get(comp, 0.0))
        assert n_out == pytest.approx(n_in, rel=1e-12)
    # 100% capture (the book's HYSYS assumption): no solids in the filtrate
    assert solid_mass_rate(filt, "calcium", MW_CA) == pytest.approx(0.0, abs=1e-15)
    # 50% wet-basis moisture: cake liquid mass equals cake solids mass
    m_cake_s = solid_mass_rate(cake, "calcium", MW_CA)
    m_cake_w = solid_mass_rate(cake, "water", MW_W)
    assert m_cake_s == pytest.approx(1200.0 / 3600.0, rel=1e-12)
    assert m_cake_w == pytest.approx(m_cake_s, rel=1e-12)
    # filtrate pressure sits on the vacuum side of the cake
    assert filt.P == pytest.approx(160e3 - 10e3)


def test_rvf_partial_capture():
    fs = rvf_fs(dict(BOOK_RVF, solids_capture=0.95))
    fs.solve()
    feed, filt = fs.streams["F"], fs.streams["L"]
    assert (solid_mass_rate(filt, "calcium", MW_CA)
            == pytest.approx(0.05 * solid_mass_rate(feed, "calcium", MW_CA),
                             rel=1e-12))


def test_mccabe_continuous_filtration_hand_calc():
    """Book Exercise 12.2 (= the McCabe 7e ch. 29 continuous-filtration CaCO3
    case): aqueous slurry of 14.7 lb solids per ft^3 of water (236 kg/m^3),
    dP = 20 in Hg, 50% wet-basis cake moisture, 100 gal/min of slurry, 5 min
    cycle, 30% submergence, 20 C / 2 bar. alpha at that dP from the McCabe
    Example 29.1 fit, 1.91e11 ft/lb = 1.283e11 m/kg (cited from memory —
    order-of-magnitude CaCO3 cake; the *equation* check below is exact).

    The in-test hand calculation recomputes c and the flux from the unit's
    own outlet streams (same property package), so it validates the
    McCabe flux equation wiring to 1e-9. The absolute area carries the PR
    liquid-water molar-volume bias (~+17% vs real water density, hence
    ~+27% on area vs a real-density hand value of ~0.25 m^2) — documented,
    not hidden: use coolprop-quality liquid densities for absolute numbers.
    """
    rho_w, rho_s = 998.2, 2800.0
    m_per_m3w = 14.7 * 0.45359237 / 0.0283168       # 235.5 kg solids/m^3 water
    q_slurry = 100.0 * 3.785411784e-3 / 60.0        # m^3/s
    q_water = q_slurry / (1.0 + m_per_m3w / rho_s)  # m^3 water /s
    m_w, m_s = q_water * rho_w * 3600.0, q_water * m_per_m3w * 3600.0  # kg/h
    dp = 20.0 * 3386.389                            # 20 in Hg, Pa
    alpha = 1.91e11 * 0.3048 / 0.45359237           # ft/lb -> m/kg

    fs = rvf_fs({"solids": "calcium", "pressure_drop": dp, "cycle_time_s": 300.0,
                 "submergence": 0.30, "cake_moisture": 0.5, "alpha": alpha},
                T=293.15, P=2e5, m_solids=m_s, m_water=m_w)
    fs.solve()
    d = fs.units["RVF"].design

    # independent recomputation of the McCabe flux equation from the outlet
    # streams: c = cake solids / filtrate volume; flux = sqrt(2 c dP f/(a mu tc))
    cake = fs.streams["C"]
    m_cake_solids = solid_mass_rate(cake, "calcium", MW_CA)
    q_f = d["filtrate_m3_s"]
    c_hand = m_cake_solids / q_f
    flux_hand = math.sqrt(2.0 * c_hand * dp * 0.30
                          / (alpha * d["filtrate_viscosity_Pa_s"] * 300.0))
    assert d["c_solids_kg_m3_filtrate"] == pytest.approx(c_hand, rel=1e-9)
    assert d["flux_m3_m2_s"] == pytest.approx(flux_hand, rel=1e-9)
    assert d["area_m2"] == pytest.approx(q_f / flux_hand, rel=1e-9)

    # databank filtrate viscosity ~ water at 20 C (Crane: 1.0 cP)
    assert d["filtrate_viscosity_Pa_s"] == pytest.approx(1.0e-3, rel=5e-2)

    # vs the all-real-properties hand value (rho_w = 998 kg/m^3): ~0.25 m^2.
    # PR's liquid-water volume bias inflates Qf and the area by ~27%.
    q_f_real = (m_w - m_s) / 3600.0 / rho_w         # 50% moisture takes m_s of water
    c_real = (m_s / 3600.0) / q_f_real
    a_real = q_f_real / math.sqrt(2.0 * c_real * dp * 0.30 / (alpha * 1.0e-3 * 300.0))
    assert a_real == pytest.approx(0.25, abs=0.01)
    assert d["area_m2"] == pytest.approx(a_real, rel=0.35)


def test_rvf_area_scales_with_sqrt_alpha():
    a1 = rvf_fs(BOOK_RVF)
    a2 = rvf_fs(dict(BOOK_RVF, alpha=4.0 * BOOK_RVF["alpha"]))
    a1.solve()
    a2.solve()
    assert (a2.units["RVF"].design["area_m2"]
            == pytest.approx(2.0 * a1.units["RVF"].design["area_m2"], rel=1e-9))


def test_rvf_typed_errors():
    cases = [
        ({k: v for k, v in BOOK_RVF.items() if k != "alpha"}, "alpha"),
        ({k: v for k, v in BOOK_RVF.items() if k != "pressure_drop"}, "pressure_drop"),
        ({k: v for k, v in BOOK_RVF.items() if k != "cycle_time_s"}, "cycle_time_s"),
        (dict(BOOK_RVF, solids_capture=0.0), "solids_capture"),
        (dict(BOOK_RVF, solids_capture=1.5), "solids_capture"),
        (dict(BOOK_RVF, cake_moisture=1.0), "cake_moisture"),
        (dict(BOOK_RVF, cake_moisture=0.99), "reduce cake_moisture"),
        (dict(BOOK_RVF, submergence=0.0), "submergence"),
        (dict(BOOK_RVF, pressure_drop=200e3), "inlet pressure"),
        (dict(BOOK_RVF, solids="benzene"), "not in the flowsheet"),
    ]
    for params, match in cases:
        with pytest.raises(SolidsOperationError, match=match):
            rvf_fs(params).solve()


def test_rvf_no_solids_in_feed_raises():
    fs = rvf_fs(BOOK_RVF, m_solids=0.0)
    with pytest.raises(SolidsOperationError, match="nothing to filter"):
        fs.solve()


def test_rvf_costing():
    fs = rvf_fs(BOOK_RVF)
    rep = fs.solve()
    pp = make_package(fs.property_package, fs.component_ids)
    sizes = size_flowsheet(fs, rep, pp)
    assert len(sizes) == 1 and sizes[0].equipment_type == "rotary_vacuum_filter"
    assert sizes[0].attribute == pytest.approx(2.4368, rel=2e-3)
    cost = cost_equipment(sizes[0])
    assert math.isfinite(cost.bare_module) and cost.bare_module > 0.0


# -- baghouse filter ------------------------------------------------------------
def baghouse_fs(params: dict, *, T=308.15, P=180e3, n=60e3 / 3600.0,
                z=None) -> Flowsheet:
    """The book's sec. 12.3 feed: 60 kmol/h air + 5 mol% sulfur, 35 C/180 kPa."""
    fs = Flowsheet(
        components=[Component("nitrogen"), Component("oxygen"), Component("sulfur")],
        property_package="thermo:PR")
    fs.add(BaghouseFilter("BH", params))
    fs.feed("F", "BH:gas_in", T=T, P=P, molar_flow=n,
            z=z or {"nitrogen": 0.95 * 0.79, "oxygen": 0.95 * 0.21, "sulfur": 0.05})
    fs.connect("G", "BH:gas_out", None)
    fs.connect("S", "BH:solids_out", None)
    return fs


def test_book_baghouse_sizing_and_capture():
    fs = baghouse_fs({"solids": "sulfur"})
    fs.solve()
    d = fs.units["BH"].design
    # gas flow ~ ideal: Q = n_gas R T / P; cloth area = Q / v_f at 0.01 m/s
    n_gas = 0.95 * 60e3 / 3600.0
    q_ideal = n_gas * 8.314462618 * 308.15 / 180e3
    assert d["Q_m3_s"] == pytest.approx(q_ideal, rel=5e-3)
    assert d["cloth_area_m2"] == pytest.approx(d["Q_m3_s"] / 0.01, rel=1e-12)
    assert d["n_bags"] == math.ceil(d["cloth_area_m2"] / (math.pi * 0.15 * 3.0))
    # inlet dust loading by hand: 3 kmol/h of S = 0.0267 kg/s over Q
    m_s = (0.05 * 60e3 / 3600.0) * MW_S
    assert d["dust_loading_kg_m3"] == pytest.approx(m_s / d["Q_m3_s"], rel=1e-12)
    # default 99.9% capture: 0.1% of the sulfur escapes with the gas
    feed, gas, sol = fs.streams["F"], fs.streams["G"], fs.streams["S"]
    assert (solid_mass_rate(gas, "sulfur", MW_S)
            == pytest.approx(1e-3 * solid_mass_rate(feed, "sulfur", MW_S), rel=1e-12))
    for comp in ("nitrogen", "oxygen", "sulfur"):
        n_in = feed.molar_flow * feed.z[comp]
        n_out = (gas.molar_flow * gas.z.get(comp, 0.0)
                 + sol.molar_flow * sol.z.get(comp, 0.0))
        assert n_out == pytest.approx(n_in, rel=1e-12)
    assert gas.P == pytest.approx(180e3 - d["dP_max_Pa"])


def test_book_baghouse_filtration_time_drag_model():
    """The book's sec. 12.3 question: the filtration time at which the dirty
    dP reaches 2 kPa. Drag model: t = (dP_max - S_E v) / (K2 c v^2). The
    book's HYSYS number (5.92e7 s) uses unpublished internal defaults and is
    NOT asserted; the model structure (linear in dP — the book's Fig. 12.18)
    and the hand-recomputed value with our cited defaults are."""
    fs = baghouse_fs({"solids": "sulfur"})
    fs.solve()
    d = fs.units["BH"].design
    t_hand = (2000.0 - d["dP_clean_Pa"]) / (5.0e4 * d["dust_loading_kg_m3"] * 0.01 ** 2)
    assert d["filtration_time_s"] == pytest.approx(t_hand, rel=1e-9)

    # linear in dP_max (book Fig. 12.18): equal dP steps give equal time steps
    times = []
    for dp_max in (2000.0, 3000.0, 4000.0):
        f = baghouse_fs({"solids": "sulfur", "dP_max": dp_max})
        f.solve()
        times.append(f.units["BH"].design["filtration_time_s"])
    assert times[2] - times[1] == pytest.approx(times[1] - times[0], rel=1e-9)
    assert times == sorted(times)


def test_baghouse_dust_free_gas_runs_forever():
    fs = baghouse_fs({"solids": "sulfur"},
                     z={"nitrogen": 0.79, "oxygen": 0.21, "sulfur": 0.0})
    fs.solve()
    assert fs.units["BH"].design["filtration_time_s"] == math.inf
    assert fs.streams["S"].molar_flow == 0.0


def test_baghouse_typed_errors():
    cases = [
        ({}, "solids"),
        ({"solids": "sulfur", "efficiency": 1.5}, "efficiency"),
        ({"solids": "sulfur", "efficiency": 0.0}, "efficiency"),
        ({"solids": "sulfur", "face_velocity": -0.01}, "face_velocity"),
        ({"solids": "sulfur", "K2": 0.0}, "K2"),
        ({"solids": "sulfur", "dP_max": 100.0}, "clean-bag"),
        ({"solids": "sulfur", "dP_max": 250e3}, "inlet pressure"),
        ({"solids": "sulfur", "bag_length_m": 0.0}, "bag_diameter_m"),
    ]
    for params, match in cases:
        with pytest.raises(SolidsOperationError, match=match):
            baghouse_fs(params).solve()


def test_baghouse_costing():
    fs = baghouse_fs({"solids": "sulfur"})
    rep = fs.solve()
    pp = make_package(fs.property_package, fs.component_ids)
    sizes = size_flowsheet(fs, rep, pp)
    assert len(sizes) == 1 and sizes[0].equipment_type == "baghouse"
    cost = cost_equipment(sizes[0])
    assert math.isfinite(cost.bare_module) and cost.bare_module > 0.0
    # EPA-linear correlation: 160 $/m^2 (2001 basis), clamped at the 10 m^2 floor
    a = max(sizes[0].attribute, 10.0)
    assert cost.purchased == pytest.approx(160.0 * a * (797.9 / 397.0), rel=1e-2)


# -- the book's cyclone + baghouse train (Exercise 12.3 / sec. 12.3.2) ----------
# HYSYS's log-probability PSD of kaolin (mean 8.343 um, sigma 5.227), as the
# 9 discrete bins HYSYS itself generated (book Fig. 12.22; mass-in-range %).
PSD_KAOLIN = [{"d_microns": d, "mass_frac": w / 100.0} for d, w in
              [(22.23e-3, 0.03), (0.1209, 0.80), (0.6577, 7.05),
               (3.577, 24.20), (19.46, 35.83), (105.8, 24.20),
               (575.6, 7.05), (3131.0, 0.80), (17030.0, 0.03)]]


def test_book_cyclone_baghouse_train():
    """Book p. 408-409 (Figs. 12.23-12.25): 10,000 kg/h air + 2000 kg/h
    kaolin at 30 C / 5 bar; a 75%-efficient cyclone leaves 500 kg/h of dust
    in the gas, and the baghouse then takes out ~100% of the rest. Kaolin
    has no critical constants in the databank, so carbon stands in as the
    dust *component* (the PSD and the 2600 kg/m^3 density are the book's
    kaolin values; mass rates are set in kg/h so the substitution only
    changes the dust molar flow, not the mass balance). The cyclone body
    (D = 0.6559 m, Stairmand HE) is sized so the Lapple model hits the
    book's 75% spec; the train then reproduces the 500 kg/h to <0.1%."""
    n_air = (10000.0 / 3600.0) / MW_AIR
    n_c = (2000.0 / 3600.0) / MW_C
    n = n_air + n_c
    zc = n_c / n
    fs = Flowsheet(
        components=[Component("nitrogen"), Component("oxygen"), Component("carbon")],
        property_package="thermo:PR")
    fs.add(Cyclone("CYC", {"solids": "carbon", "particle_density": 2600.0,
                           "geometry": "Stairmand_HE", "body_diameter": 0.6559,
                           "psd": PSD_KAOLIN}))
    # the PSD does not propagate (v1 model): the baghouse re-specifies its own
    # capture as a bulk efficiency (default 99.9%)
    fs.add(BaghouseFilter("BH", {"solids": "carbon"}))
    fs.feed("F", "CYC:gas_in", T=303.15, P=5e5, molar_flow=n,
            z={"nitrogen": (1 - zc) * 0.79, "oxygen": (1 - zc) * 0.21, "carbon": zc})
    fs.connect("G", "CYC:gas_out", "BH:gas_in")
    fs.connect("S1", "CYC:solids_out", None)
    fs.connect("G2", "BH:gas_out", None)
    fs.connect("S2", "BH:solids_out", None)
    rep = fs.solve()
    assert rep.converged

    eta = fs.units["CYC"].design["overall_efficiency"]
    assert eta == pytest.approx(0.75, abs=2e-3)              # the book's spec
    # book p. 409: 500 kg/h of dust remains in the cyclone's gas outlet
    m_gas_dust = solid_mass_rate(fs.streams["G"], "carbon", MW_C) * 3600.0
    assert m_gas_dust == pytest.approx(500.0, rel=3e-3)
    # ... and the baghouse removes ~100% of it (book: 0.0 kg/h; ours: 0.1%)
    m_final = solid_mass_rate(fs.streams["G2"], "carbon", MW_C) * 3600.0
    assert m_final == pytest.approx(0.5, rel=3e-3)
    assert m_final < 0.005 * 500.0
    # train-wide dust balance closes exactly
    m_caught = (solid_mass_rate(fs.streams["S1"], "carbon", MW_C)
                + solid_mass_rate(fs.streams["S2"], "carbon", MW_C)) * 3600.0
    assert m_caught + m_final == pytest.approx(2000.0, rel=1e-12)
