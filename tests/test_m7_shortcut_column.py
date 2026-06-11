"""M7 tests: the ShortcutColumn (Fenske-Underwood-Gilliland) unit op.

Validation reference: the classic equimolar benzene/toluene split at 1 atm.
The literature average relative volatility for benzene/toluene at 101.3 kPa is
alpha ~= 2.4 (Wankat, *Separation Process Engineering*, 3e, Ch. 7 McCabe-Thiele
examples; consistent with Perry's 8e Sec. 13 Antoine K-ratios over 80-110 C).
With alpha = 2.40, q = 1 (saturated-liquid feed), z = 0.50/0.50 and 95%/95% key
recoveries, the closed-form shortcut relations give:

  * Fenske:    N_min = ln[(0.95/0.05)^2] / ln 2.4            = 6.73
  * Underwood: theta from 2.4(0.5)/(2.4-theta) + 0.5/(1-theta) = 0
               -> theta = 24/17 = 1.4118
               R_min = 2.4(0.95)/(2.4-theta) + 0.05/(1-theta) - 1 = 1.186
  * Gilliland (Molokanov), R = 1.3 R_min = 1.542:
               X = 0.140, Y = 0.5146, N = (N_min + Y)/(1 - Y) = 14.9
  * Kirkbride: symmetric split -> N_R/N_S = 1, feed near mid-column.

The engine derives alpha from the PR EOS (geometric mean of top/bottom values)
rather than taking 2.40 as given, so the comparison tolerance is 10% — generous
for shortcut methods, tight enough to catch a wrong equation.
"""
import math

import pytest

from caldyr.core import Component, Flowsheet
from caldyr.economics import TEAConfig, analyze
from caldyr.io import from_dict, to_dict
from caldyr.unitops import ShortcutColumn, ShortcutColumnError

P_ATM = 101325.0

# Hand-calculated references (alpha = 2.40, q = 1; see module docstring).
REF_N_MIN = 6.73
REF_R_MIN = 1.186
REF_N = 14.9


def bt_column(**param_overrides) -> Flowsheet:
    """Equimolar benzene/toluene feed (100 mol/s, saturated liquid at 1 atm)
    into a shortcut column with 95%/95% key recoveries, R = 1.3 R_min."""
    params = {"light_key": "benzene", "heavy_key": "toluene",
              "recovery_light": 0.95, "recovery_heavy": 0.95,
              "rr_factor": 1.3, "P": P_ATM}
    params.update(param_overrides)
    fs = Flowsheet(components=[Component("benzene"), Component("toluene")],
                   property_package="thermo:PR")
    fs.add(ShortcutColumn("COL", params))
    # 365 K is just below the 50/50 bubble point (365.25 K) -> q = 1.
    fs.feed("FEED", "COL:in1", T=365.0, P=P_ATM, molar_flow=100.0,
            z={"benzene": 0.5, "toluene": 0.5})
    fs.connect("DIST", "COL:distillate", None)
    fs.connect("BOT", "COL:bottoms", None)
    fs.connect("QC", "COL:condenser_duty", None)
    fs.connect("QR", "COL:reboiler_duty", None)
    return fs


def btx_column() -> Flowsheet:
    """Benzene/toluene split with a heavy non-key (p-xylene) to distribute."""
    fs = Flowsheet(
        components=[Component("benzene"), Component("toluene"), Component("p-xylene")],
        property_package="thermo:PR")
    fs.add(ShortcutColumn("COL", {"light_key": "benzene", "heavy_key": "toluene",
                                  "recovery_light": 0.99, "recovery_heavy": 0.98,
                                  "rr_factor": 1.3, "P": P_ATM}))
    fs.feed("FEED", "COL:in1", T=370.0, P=P_ATM, molar_flow=100.0,
            z={"benzene": 0.45, "toluene": 0.40, "p-xylene": 0.15})
    fs.connect("DIST", "COL:distillate", None)
    fs.connect("BOT", "COL:bottoms", None)
    fs.connect("QC", "COL:condenser_duty", None)
    fs.connect("QR", "COL:reboiler_duty", None)
    return fs


# -- FUG vs the hand-calculated reference -----------------------------------
def test_fug_matches_benzene_toluene_reference():
    fs = bt_column()
    assert fs.solve().converged
    d = fs.units["COL"].design

    assert d["q"] == pytest.approx(1.0, abs=0.02)           # saturated-liquid feed
    # PR-derived alpha should land on the literature ~2.4 average.
    assert d["alpha"]["benzene"] == pytest.approx(2.40, rel=0.05)
    assert d["N_min"] == pytest.approx(REF_N_MIN, rel=0.10)
    assert d["R_min"] == pytest.approx(REF_R_MIN, rel=0.10)
    assert d["R"] == pytest.approx(1.3 * d["R_min"], rel=1e-12)
    assert d["N"] == pytest.approx(REF_N, rel=0.10)
    # Symmetric split -> Kirkbride puts the feed near mid-column.
    assert d["feed_stage"] == pytest.approx(d["N"] / 2, abs=1.5)
    # Products: 95% recoveries on a 50/50 feed -> 50 mol/s at 95% purity each.
    assert fs.streams["DIST"].molar_flow == pytest.approx(50.0, rel=1e-12)
    assert fs.streams["DIST"].z["benzene"] == pytest.approx(0.95, rel=1e-12)
    assert fs.streams["BOT"].z["toluene"] == pytest.approx(0.95, rel=1e-12)


def test_product_thermal_states_and_duty_signs():
    fs = bt_column()
    rep = fs.solve()
    dist, bot = fs.streams["DIST"], fs.streams["BOT"]
    # Total condenser: liquid distillate at its bubble point; liquid bottoms at
    # theirs; bottoms hotter than distillate; both between the pure-component
    # boiling points (benzene 353.2 K, toluene 383.8 K at 1 atm).
    assert dist.phase == "liquid" and bot.phase == "liquid"
    assert 352.0 < dist.T < 356.0
    assert 379.0 < bot.T < 384.0
    assert dist.T < bot.T
    # Heater sign convention: condenser removes heat (<0), reboiler adds (>0),
    # and the magnitudes are ~ V * latent heat (~30 kJ/mol for aromatics).
    assert rep.duties["QC"] < 0 < rep.duties["QR"]
    v_top = fs.units["COL"].design["V_top"]
    assert abs(rep.duties["QC"]) == pytest.approx(v_top * 31e3, rel=0.2)


def test_partial_condenser_gives_vapor_distillate_at_dew_point():
    fs_tot = bt_column()
    fs_tot.solve()
    fs_par = bt_column(partial_condenser=True)
    rep = fs_par.solve()
    dist = fs_par.streams["DIST"]
    assert dist.phase == "vapor"
    assert dist.T > fs_tot.streams["DIST"].T          # dew point > bubble point
    # Condensing only the reflux costs less than condensing all the overhead.
    assert abs(rep.duties["QC"]) < abs(fs_tot.solve().duties["QC"])


# -- conservation ------------------------------------------------------------
def test_component_mass_balance_closes_to_machine_precision():
    fs = btx_column()
    assert fs.solve().converged
    feed, dist, bot = fs.streams["FEED"], fs.streams["DIST"], fs.streams["BOT"]
    for c in feed.components:
        n_in = feed.molar_flow * feed.z[c]
        n_out = dist.molar_flow * dist.z[c] + bot.molar_flow * bot.z[c]
        assert n_out == pytest.approx(n_in, rel=1e-12, abs=1e-12)


def test_energy_balance_closes_including_both_duties():
    """F h_F + Q_cond + Q_reb = D h_D + B h_B — exact, because the engine
    enthalpy is formation-inclusive (absolute) and the reboiler duty is defined
    to close the balance."""
    fs = btx_column()
    rep = fs.solve()
    feed, dist, bot = fs.streams["FEED"], fs.streams["DIST"], fs.streams["BOT"]
    h_in = feed.molar_flow * feed.H + rep.duties["QC"] + rep.duties["QR"]
    h_out = dist.molar_flow * dist.H + bot.molar_flow * bot.H
    assert h_out == pytest.approx(h_in, rel=1e-12)


# -- non-key distribution -----------------------------------------------------
def test_heavy_nonkey_goes_overwhelmingly_to_bottoms():
    fs = btx_column()
    fs.solve()
    d = fs.units["COL"].design
    f_px = 15.0
    assert d["bottoms_flows"]["p-xylene"] / f_px > 0.999    # ~all to bottoms
    assert fs.streams["DIST"].z["p-xylene"] < 1e-4
    # ... and it still appears (in trace) in the distillate, per Fenske.
    assert d["distillate_flows"]["p-xylene"] > 0.0


# -- `.flow` round-trip --------------------------------------------------------
def test_flow_round_trip_is_exact():
    fs = btx_column()
    fs.solve()
    doc = to_dict(fs)
    again = to_dict(from_dict(doc))
    assert doc == again
    # The reloaded flowsheet re-solves to the same design.
    fs2 = from_dict(doc)
    fs2.solve()
    assert fs2.units["COL"].design["N"] == pytest.approx(fs.units["COL"].design["N"])


# -- both solver backends agree -------------------------------------------------
def test_sequential_and_equation_oriented_backends_agree():
    fa = btx_column()
    ra = fa.solve(backend="sequential", tol=1e-9)
    fb = btx_column()
    rb = fb.solve(backend="equation_oriented", tol=1e-10)
    assert ra.converged and rb.converged
    for sid in ("DIST", "BOT"):
        sa, sb = fa.streams[sid], fb.streams[sid]
        assert math.isclose(sa.molar_flow, sb.molar_flow, rel_tol=1e-6)
        assert math.isclose(sa.T, sb.T, rel_tol=1e-6)
        for c in sa.components:
            assert math.isclose(sa.z.get(c, 0.0), sb.z.get(c, 0.0),
                                rel_tol=1e-6, abs_tol=1e-9)
    for qid in ("QC", "QR"):
        assert math.isclose(ra.duties[qid], rb.duties[qid], rel_tol=1e-5, abs_tol=1.0)


# -- economics -------------------------------------------------------------------
def test_economics_analyze_costs_the_column():
    fs = btx_column()
    rep = fs.solve()
    res = analyze(fs, rep, TEAConfig(product_component="benzene",
                                     product_min_fraction=0.9))
    by_id = {s.unit_id: s for s in res.sizes}
    # Tower + trays + condenser + reboiler all sized.
    assert set(by_id) == {"COL", "COL.trays", "COL.condenser", "COL.reboiler"}
    assert by_id["COL"].equipment_type == "vessel_vertical"
    assert by_id["COL.trays"].quantity > 1
    # The duties map to a cooling and a heating utility (opex like a Heater).
    assert by_id["COL.condenser"].utility == "cooling_water"
    util_kind = {"cooling_water": "cool", "low_pressure_steam": "heat",
                 "fired_heat": "heat"}
    assert util_kind[by_id["COL.reboiler"].utility] == "heat"
    assert res.opex.breakdown["utilities"]["cooling_water"] > 0

    costs = {c.unit_id: c for c in res.costs}
    tower = costs["COL"]
    assert math.isfinite(tower.bare_module) and tower.bare_module > 0
    assert costs["COL.trays"].bare_module > 0
    assert math.isfinite(res.profitability.lcop) and res.profitability.lcop > 0
    assert res.capital.tci > res.capital.isbl > 0


# -- typed, actionable errors ------------------------------------------------------
def test_unknown_key_raises():
    fs = bt_column(light_key="hexane")
    with pytest.raises(ShortcutColumnError, match="not in the flowsheet component list"):
        fs.solve()


def test_inverted_keys_raise():
    fs = bt_column(light_key="toluene", heavy_key="benzene")
    with pytest.raises(ShortcutColumnError, match="not more volatile"):
        fs.solve()


def test_missing_keys_raise():
    fs = bt_column(light_key=None)
    with pytest.raises(ShortcutColumnError, match="required"):
        fs.solve()


def test_infeasible_recoveries_raise():
    fs = bt_column(recovery_light=0.10, recovery_heavy=0.10)
    with pytest.raises(ShortcutColumnError, match="separation factor"):
        fs.solve()


def test_recovery_of_one_raises():
    fs = bt_column(recovery_light=1.0)
    with pytest.raises(ShortcutColumnError, match="infinite stages"):
        fs.solve()


def test_rr_factor_at_or_below_one_raises():
    fs = bt_column(rr_factor=1.0)
    with pytest.raises(ShortcutColumnError, match="rr_factor"):
        fs.solve()
