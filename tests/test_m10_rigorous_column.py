"""M10 tests: the RigorousColumn (tray-by-tray MESH, bubble-point/Wang-Henke).

Validation strategy (sources cited per test):

1. **Consistency with FUG** (Fenske 1932 / Underwood 1948 / Gilliland 1940, as
   implemented and Wankat-validated in test_m7_shortcut_column): a rigorous
   column given the *same* reflux ratio, stage count, feed stage and
   distillate rate as a converged ShortcutColumn design must reproduce the
   FUG key recoveries within a few percent. FUG itself is a +/-10% class
   method, so agreement at the few-percent level on the easy benzene/toluene
   system is the expected behavior, not luck.
2. **Closed-form monotonicity** (any mass-transfer text, e.g. Wankat,
   *Separation Process Engineering* 3e ch. 7; Seader, Henley & Roper,
   *Separation Process Principles* 3e ch. 10): at fixed reflux more stages
   sharpen the split, and at fixed stages more reflux sharpens the split.
3. **Structural validation of a ternary case** (benzene/toluene/p-xylene).
   No published rigorous stage-by-stage table is faithfully reproducible
   here: textbook MESH examples (e.g. Seader 3e ch. 10 examples) are built on
   their own K-value charts (DePriester/ideal-Raoult) and enthalpy
   correlations, so matching their tables would test transcription of their
   thermo, not this column. Instead the ternary case is validated
   *structurally*: recoveries ordered by volatility, a strictly increasing
   temperature profile, monotone light-component composition profiles, and
   vapor/liquid traffic consistent with the specified reflux — each a known
   qualitative property of a simple column (Seader 3e ch. 10, Kister
   *Distillation Design* ch. 4).
4. Conservation, `.flow` round-trip, both solver backends, profile-array
   integrity, typed errors, and economics (sized exactly like a
   ShortcutColumn).
"""
import math

import pytest

from caldyr.core import Component, Flowsheet
from caldyr.economics import TEAConfig, analyze
from caldyr.io import from_dict, to_dict
from caldyr.unitops import RigorousColumn, RigorousColumnError, ShortcutColumn

P_ATM = 101325.0


def bt_shortcut() -> Flowsheet:
    """The Wankat-validated FUG benzene/toluene column from the M7 tests:
    equimolar saturated-liquid feed, 95%/95% key recoveries, R = 1.3 R_min."""
    fs = Flowsheet(components=[Component("benzene"), Component("toluene")],
                   property_package="thermo:PR")
    fs.add(ShortcutColumn("COL", {"light_key": "benzene", "heavy_key": "toluene",
                                  "recovery_light": 0.95, "recovery_heavy": 0.95,
                                  "rr_factor": 1.3, "P": P_ATM}))
    fs.feed("FEED", "COL:in1", T=365.0, P=P_ATM, molar_flow=100.0,
            z={"benzene": 0.5, "toluene": 0.5})
    for sid, port in (("DIST", "distillate"), ("BOT", "bottoms"),
                      ("QC", "condenser_duty"), ("QR", "reboiler_duty")):
        fs.connect(sid, f"COL:{port}", None)
    return fs


def bt_rigorous(**param_overrides) -> Flowsheet:
    """Benzene/toluene rigorous column with the FUG-equivalent default specs
    (16 stages incl. condenser+reboiler, feed stage 8, R=1.545, D=50 mol/s —
    see test_matches_fug for the mapping)."""
    params = {"n_stages": 16, "feed_stage": 8, "reflux_ratio": 1.545,
              "distillate_rate": 50.0, "P": P_ATM}
    params.update(param_overrides)
    fs = Flowsheet(components=[Component("benzene"), Component("toluene")],
                   property_package="thermo:PR")
    fs.add(RigorousColumn("COL", params))
    fs.feed("FEED", "COL:in1", T=365.0, P=P_ATM, molar_flow=100.0,
            z={"benzene": 0.5, "toluene": 0.5})
    for sid, port in (("DIST", "distillate"), ("BOT", "bottoms"),
                      ("QC", "condenser_duty"), ("QR", "reboiler_duty")):
        fs.connect(sid, f"COL:{port}", None)
    return fs


def btx_rigorous() -> Flowsheet:
    """Ternary benzene/toluene/p-xylene rigorous column (structural case)."""
    fs = Flowsheet(
        components=[Component("benzene"), Component("toluene"), Component("p-xylene")],
        property_package="thermo:PR")
    fs.add(RigorousColumn("COL", {"n_stages": 16, "feed_stage": 8,
                                  "reflux_ratio": 2.0, "distillate_to_feed": 0.45,
                                  "P": P_ATM}))
    fs.feed("FEED", "COL:in1", T=370.0, P=P_ATM, molar_flow=100.0,
            z={"benzene": 0.45, "toluene": 0.40, "p-xylene": 0.15})
    for sid, port in (("DIST", "distillate"), ("BOT", "bottoms"),
                      ("QC", "condenser_duty"), ("QR", "reboiler_duty")):
        fs.connect(sid, f"COL:{port}", None)
    return fs


def lk_recovery(fs) -> float:
    """Fraction of the feed benzene recovered to the distillate."""
    feed, dist = fs.streams["FEED"], fs.streams["DIST"]
    return (dist.molar_flow * dist.z["benzene"]) / (feed.molar_flow * feed.z["benzene"])


# -- 1. rigorous vs FUG -------------------------------------------------------
def test_matches_fug_benzene_toluene():
    """At the FUG design point (same R, stages, feed stage, D) the rigorous
    recoveries must land within a few percent of the FUG targets. Stage-count
    mapping: FUG's N counts the reboiler but not the total condenser, so
    n_stages = round(N) + 1 and feed_stage = FUG feed stage + 1 (module
    docstring of caldyr.unitops.rigorous_column)."""
    fs_fug = bt_shortcut()
    assert fs_fug.solve().converged
    d_fug = fs_fug.units["COL"].design

    fs_rig = bt_rigorous(n_stages=round(d_fug["N"]) + 1,
                         feed_stage=d_fug["feed_stage"] + 1,
                         reflux_ratio=d_fug["R"],
                         distillate_rate=d_fug["D"])
    assert fs_rig.solve().converged
    d_rig = fs_rig.units["COL"].design

    # FUG asked for 95%/95% key recoveries; the rigorous column at the same
    # design must agree to a few percent (FUG itself is a ~10% method).
    assert lk_recovery(fs_rig) == pytest.approx(0.95, abs=0.03)
    assert d_rig["x_B"]["toluene"] == pytest.approx(d_fug["x_B"]["toluene"], abs=0.03)
    # Duties agree to ~10% (FUG assumes constant molal overflow; MESH doesn't).
    assert d_rig["Q_condenser"] == pytest.approx(d_fug["Q_condenser"], rel=0.10)
    assert d_rig["Q_reboiler"] == pytest.approx(d_fug["Q_reboiler"], rel=0.10)
    # Product temperatures bracket the pure boiling points (benzene 353.2 K,
    # toluene 383.8 K at 1 atm; NIST) and the rigorous endpoints agree with FUG.
    assert d_rig["T_top"] == pytest.approx(d_fug["T_top"], abs=1.0)
    assert d_rig["T_bottom"] == pytest.approx(d_fug["T_bottom"], abs=1.5)


def test_more_stages_at_same_reflux_sharpen_the_split():
    recoveries = []
    for n in (10, 16, 24):
        fs = bt_rigorous(n_stages=n, feed_stage=n // 2)
        assert fs.solve().converged
        recoveries.append(lk_recovery(fs))
    assert recoveries[0] < recoveries[1] < recoveries[2]


def test_higher_reflux_at_same_stages_sharpens_the_split():
    recoveries = []
    for r in (1.0, 1.545, 3.0):
        fs = bt_rigorous(reflux_ratio=r)
        assert fs.solve().converged
        recoveries.append(lk_recovery(fs))
    assert recoveries[0] < recoveries[1] < recoveries[2]


# -- 2. ternary case, structural ----------------------------------------------
def test_ternary_structure_recoveries_ordered_by_volatility():
    """See the module docstring: no faithful textbook MESH table exists for
    our PR thermo, so the ternary case asserts the *structure* a simple column
    must have — recoveries to the distillate strictly ordered by volatility
    (benzene > toluene > p-xylene at 1 atm; their normal boiling points are
    353.2 / 383.8 / 411.5 K)."""
    fs = btx_rigorous()
    assert fs.solve().converged
    d = fs.units["COL"].design
    feed = fs.streams["FEED"]
    rec = {c: d["distillate_flows"][c] / (feed.molar_flow * feed.z[c])
           for c in ("benzene", "toluene", "p-xylene")}
    assert rec["benzene"] > 0.9 > rec["toluene"] > rec["p-xylene"]
    assert rec["p-xylene"] < 0.01


def test_ternary_profiles_are_monotone():
    """Temperature must increase strictly top->bottom and the lightest
    component's liquid fraction must fall monotonically down the column (no
    feed-composition pinch reversals for this wide-alpha system)."""
    fs = btx_rigorous()
    fs.solve()
    d = fs.units["COL"].design
    T = d["T_profile"]
    assert all(a < b for a, b in zip(T, T[1:]))
    xb = [row["benzene"] for row in d["x_profile"]]
    assert all(a >= b for a, b in zip(xb, xb[1:]))
    # The vapor is always enriched in the lightest component over the liquid.
    yb = [row["benzene"] for row in d["y_profile"]]
    assert all(y > x for y, x in zip(yb, xb))


# -- 3. conservation and result integrity ---------------------------------------
def test_component_mass_balance_closes_to_machine_precision():
    fs = btx_rigorous()
    assert fs.solve().converged
    feed, dist, bot = fs.streams["FEED"], fs.streams["DIST"], fs.streams["BOT"]
    for c in feed.components:
        n_in = feed.molar_flow * feed.z[c]
        n_out = dist.molar_flow * dist.z[c] + bot.molar_flow * bot.z[c]
        assert n_out == pytest.approx(n_in, rel=1e-12, abs=1e-12)


def test_energy_balance_closes_including_both_duties():
    """F h_F + Q_cond + Q_reb = D h_D + B h_B — exact, because the engine
    enthalpy is formation-inclusive (absolute) and the reboiler duty is
    defined to close the balance; the independent stage-N reboiler balance
    must agree to the MESH tolerance (design['energy_residual_rel'])."""
    fs = btx_rigorous()
    rep = fs.solve()
    feed, dist, bot = fs.streams["FEED"], fs.streams["DIST"], fs.streams["BOT"]
    h_in = feed.molar_flow * feed.H + rep.duties["QC"] + rep.duties["QR"]
    h_out = dist.molar_flow * dist.H + bot.molar_flow * bot.H
    assert h_out == pytest.approx(h_in, rel=1e-12)
    assert fs.units["COL"].design["energy_residual_rel"] < 1e-3
    assert rep.duties["QC"] < 0 < rep.duties["QR"]


def test_profile_arrays_present_and_consistent():
    fs = bt_rigorous()
    assert fs.solve().converged
    d = fs.units["COL"].design
    n = d["n_stages"]
    for key in ("T_profile", "P_profile", "L_profile", "V_profile",
                "x_profile", "y_profile"):
        assert isinstance(d[key], list) and len(d[key]) == n
    for row in d["x_profile"] + d["y_profile"]:
        assert sum(row.values()) == pytest.approx(1.0, abs=1e-9)
    # Traffic honors the specs: L_1 = R*D, V_2 = D(1+R) (total condenser:
    # V_1 = 0), L_N = B.
    assert d["V_profile"][0] == 0.0
    assert d["L_profile"][0] == pytest.approx(d["R"] * d["D"], rel=1e-9)
    assert d["V_profile"][1] == pytest.approx(d["D"] * (1 + d["R"]), rel=1e-9)
    assert d["L_profile"][-1] == pytest.approx(d["B"], rel=1e-12)
    # FUG-compatible sizing keys are published ("N" excludes the condenser).
    assert d["N"] == n - 1
    assert d["V_top"] == d["V_profile"][1]


def test_partial_condenser_gives_vapor_distillate():
    fs_tot = bt_rigorous()
    rep_tot = fs_tot.solve()
    fs_par = bt_rigorous(partial_condenser=True)
    rep_par = fs_par.solve()
    dist = fs_par.streams["DIST"]
    assert dist.phase == "vapor"
    assert dist.T > fs_tot.streams["DIST"].T          # dew point > bubble point
    # Condensing only the reflux costs less than condensing all the overhead.
    assert abs(rep_par.duties["QC"]) < abs(rep_tot.duties["QC"])


def test_distillate_to_feed_equals_distillate_rate():
    fa = bt_rigorous()
    fa.solve()
    fb = bt_rigorous(distillate_rate=None, distillate_to_feed=0.5)
    fb.solve()
    assert fb.streams["DIST"].molar_flow == pytest.approx(50.0, rel=1e-12)
    assert fb.streams["DIST"].z["benzene"] == pytest.approx(
        fa.streams["DIST"].z["benzene"], rel=1e-9)


def test_stage_pressure_profile_raises_bottom_temperature():
    fs0 = bt_rigorous()
    fs0.solve()
    fs1 = bt_rigorous(dP_stage=700.0)
    fs1.solve()
    d1 = fs1.units["COL"].design
    assert d1["P_profile"][-1] == pytest.approx(P_ATM + 15 * 700.0)
    assert d1["T_bottom"] > fs0.units["COL"].design["T_bottom"] + 1.0


# -- `.flow` round-trip ----------------------------------------------------------
def test_flow_round_trip_is_exact():
    fs = btx_rigorous()
    fs.solve()
    doc = to_dict(fs)
    again = to_dict(from_dict(doc))
    assert doc == again
    # The reloaded flowsheet re-solves to the same separation.
    fs2 = from_dict(doc)
    fs2.solve()
    assert fs2.streams["DIST"].z["benzene"] == pytest.approx(
        fs.streams["DIST"].z["benzene"], rel=1e-9)


# -- both solver backends agree ---------------------------------------------------
def test_sequential_and_equation_oriented_backends_agree():
    fa = bt_rigorous()
    ra = fa.solve(backend="sequential", tol=1e-9)
    fb = bt_rigorous()
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


# -- economics (sized/costed exactly like a ShortcutColumn) ------------------------
def test_economics_analyze_costs_the_column():
    fs = bt_rigorous()
    rep = fs.solve()
    res = analyze(fs, rep, TEAConfig(product_component="benzene",
                                     product_min_fraction=0.9))
    by_id = {s.unit_id: s for s in res.sizes}
    assert set(by_id) == {"COL", "COL.trays", "COL.condenser", "COL.reboiler"}
    assert by_id["COL"].equipment_type == "vessel_vertical"
    # n_stages=16 incl. condenser -> 15 equilibrium stages / 0.7 efficiency
    # = 22 real trays (same convention as ShortcutColumn).
    assert by_id["COL.trays"].quantity == math.ceil(15 / 0.7)
    assert by_id["COL.condenser"].utility == "cooling_water"
    costs = {c.unit_id: c for c in res.costs}
    assert costs["COL"].bare_module > 0
    assert costs["COL.trays"].bare_module > 0
    assert math.isfinite(res.profitability.lcop) and res.profitability.lcop > 0


# -- typed, actionable errors -------------------------------------------------------
def test_feed_stage_out_of_range_raises():
    for bad in (1, 16, 0, 99):
        fs = bt_rigorous(feed_stage=bad)
        with pytest.raises(RigorousColumnError, match="out of range"):
            fs.solve()


def test_nonpositive_reflux_raises():
    fs = bt_rigorous(reflux_ratio=0.0)
    with pytest.raises(RigorousColumnError, match="reflux_ratio"):
        fs.solve()
    fs = bt_rigorous(reflux_ratio=-2.0)
    with pytest.raises(RigorousColumnError, match="reflux_ratio"):
        fs.solve()


def test_distillate_rate_at_or_above_feed_raises():
    for bad in (100.0, 150.0, 0.0):
        fs = bt_rigorous(distillate_rate=bad)
        with pytest.raises(RigorousColumnError, match="strictly between"):
            fs.solve()


def test_both_or_neither_distillate_spec_raises():
    fs = bt_rigorous(distillate_to_feed=0.5)            # both given
    with pytest.raises(RigorousColumnError, match="exactly one"):
        fs.solve()
    fs = bt_rigorous(distillate_rate=None)              # neither given
    with pytest.raises(RigorousColumnError, match="exactly one"):
        fs.solve()


def test_too_few_stages_raises():
    fs = bt_rigorous(n_stages=2, feed_stage=2)
    with pytest.raises(RigorousColumnError, match="n_stages"):
        fs.solve()


def test_missing_stage_params_raise():
    fs = bt_rigorous(n_stages=None)
    with pytest.raises(RigorousColumnError, match="n_stages"):
        fs.solve()
