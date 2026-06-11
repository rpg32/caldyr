"""M8 tests: flowsheet-level logical ops (Set / Adjust), solver hints (the
Recycle-block equivalent: tear guesses + tolerance override), and the
post-solve balance report.

These are flowsheet-level constructs, not unit ops: they live in
``Flowsheet.logical`` / ``Flowsheet.solver_hints``, persist in `.flow` under
the ``"logical"`` / ``"solver_hints"`` keys, and are interpreted by
``fs.solve()`` in a documented order — Sets are applied before every inner
solve; Adjusts wrap the whole solve in an outer root find (Brent for one,
scipy hybr for several), so they are backend-agnostic.

The physics here is intentionally simple (a heater feeding an adiabatic flash;
the cited M1/M2 example flowsheets) — what is under test is the logical-op
machinery, which must hit specs within tolerance and keep balances closed.
"""
import pytest

from caldyr.core import Component, Flowsheet
from caldyr.io import from_dict, to_dict
from caldyr.solver import AdjustError, LogicalOpError, balance_report
from caldyr.unitops import EquilibriumReactor, FlashDrum, Heater, Mixer, Splitter

P_ATM = 101325.0


def heater_flash() -> Flowsheet:
    """FEED -> HEAT -> FL (adiabatic flash at 1 atm); heater T_out drives the
    flash vapor make."""
    fs = Flowsheet(components=[Component("n-pentane"), Component("n-octane")],
                   property_package="thermo:PR")
    fs.add(Heater("HEAT", {"T_out": 350.0}))
    fs.add(FlashDrum("FL", {"P": P_ATM}))
    fs.feed("FEED", "HEAT:in1", T=300.0, P=P_ATM, molar_flow=10.0,
            z={"n-pentane": 0.5, "n-octane": 0.5})
    fs.connect("S1", "HEAT:out", "FL:in1")
    fs.connect("VAP", "FL:vapor", None)
    fs.connect("LIQ", "FL:liquid", None)
    fs.connect("QH", "HEAT:duty", None)
    fs.connect("QF", "FL:duty", None)
    return fs


def flash_recycle() -> Flowsheet:
    """The M1 flash-recycle example (examples/02_flash_recycle.py)."""
    fs = Flowsheet(components=[Component("n-pentane"), Component("n-octane")],
                   property_package="thermo:PR")
    fs.add(Mixer("MIX", {"dP": 0.0}))
    fs.add(FlashDrum("FL", {"T": 360.0, "P": P_ATM}))
    fs.add(Splitter("SP", {"split": 0.6}))
    fs.feed("FEED", "MIX:in1", T=330.0, P=P_ATM, molar_flow=10.0,
            z={"n-pentane": 0.5, "n-octane": 0.5})
    fs.connect("MIXOUT", "MIX:out", "FL:in1")
    fs.connect("VAP", "FL:vapor", None)
    fs.connect("LIQ", "FL:liquid", "SP:in1")
    fs.connect("RECY", "SP:out1", "MIX:in2")
    fs.connect("BOTTOMS", "SP:out2", None)
    fs.connect("Q", "FL:duty", None)
    return fs


def ammonia_loop() -> Flowsheet:
    """The M2 ammonia-loop example (examples/04_ammonia_loop.py)."""
    rxn = {"stoich": {"nitrogen": -1, "hydrogen": -3, "ammonia": 2}, "key": "nitrogen"}
    fs = Flowsheet(
        components=[Component(c) for c in ("nitrogen", "hydrogen", "ammonia", "argon")],
        property_package="thermo:PR")
    fs.add(Mixer("MIX", {"dP": 0.0}))
    fs.add(Heater("PREHEAT", {"T_out": 673.15}))
    fs.add(EquilibriumReactor("RXN", {"reaction": rxn, "T": 673.15}))
    fs.add(Heater("COOL", {"T_out": 250.0}))
    fs.add(FlashDrum("SEP", {"T": 250.0, "P": 2.0e7}))
    fs.add(Splitter("SPLIT", {"split": 0.90}))
    fs.feed("MAKEUP", "MIX:in1", T=300.0, P=2.0e7, molar_flow=100.0,
            z={"nitrogen": 0.2475, "hydrogen": 0.7425, "ammonia": 0.0, "argon": 0.01})
    fs.connect("S1", "MIX:out", "PREHEAT:in1")
    fs.connect("S2", "PREHEAT:out", "RXN:in1")
    fs.connect("S3", "RXN:out", "COOL:in1")
    fs.connect("S4", "COOL:out", "SEP:in1")
    fs.connect("PRODUCT", "SEP:liquid", None)
    fs.connect("VAP", "SEP:vapor", "SPLIT:in1")
    fs.connect("RECYCLE", "SPLIT:out1", "MIX:in2")
    fs.connect("PURGE", "SPLIT:out2", None)
    for unit in ("PREHEAT", "RXN", "COOL", "SEP"):
        fs.connect(f"Q_{unit}", f"{unit}:duty", None)
    return fs


VAP_ADJUST = {"type": "adjust", "vary": ["HEAT", "T_out"], "bounds": [320.0, 420.0],
              "spec": {"type": "molar_flow", "stream": "VAP"},
              "value": 5.0, "tolerance": 1e-6}


# -- Adjust ---------------------------------------------------------------------
def test_adjust_hits_a_downstream_flash_vapor_target():
    fs = heater_flash()
    fs.logical.append(dict(VAP_ADJUST))
    rep = fs.solve()
    assert rep.converged
    assert fs.streams["VAP"].molar_flow == pytest.approx(5.0, abs=1e-6)
    assert 320.0 < fs.units["HEAT"].params["T_out"] < 420.0
    # Convergence is narrated in the report messages.
    assert any("adjust: HEAT.T_out" in m for m in rep.messages)
    assert any("flowsheet solves" in m for m in rep.messages)


def test_adjust_is_backend_agnostic():
    """Adjust wraps solve(), so it must give the same answer under both
    backends (which agree on the physics to ~1e-9)."""
    results = {}
    for backend in ("sequential", "equation_oriented"):
        fs = heater_flash()
        fs.logical.append(dict(VAP_ADJUST))
        rep = fs.solve(backend=backend)
        assert rep.converged
        results[backend] = fs.units["HEAT"].params["T_out"]
        assert fs.streams["VAP"].molar_flow == pytest.approx(5.0, abs=1e-5)
    assert results["sequential"] == pytest.approx(results["equation_oriented"], rel=1e-6)


def test_two_simultaneous_adjusts():
    """Two Adjusts -> a vector root find: heater T_out and flash P together hit
    a vapor-flow and a vapor-temperature spec."""
    fs = heater_flash()
    fs.logical.extend([
        dict(VAP_ADJUST),
        {"type": "adjust", "vary": ["FL", "P"], "bounds": [0.5e5, 2.0e5],
         "spec": {"type": "T", "stream": "VAP"}, "value": 345.0, "tolerance": 1e-3},
    ])
    rep = fs.solve()
    assert rep.converged
    assert fs.streams["VAP"].molar_flow == pytest.approx(5.0, abs=1e-5)
    assert fs.streams["VAP"].T == pytest.approx(345.0, abs=1e-3)


def test_adjust_unreachable_target_raises_with_diagnostics():
    fs = heater_flash()
    fs.logical.append({**VAP_ADJUST, "value": 11.0})   # > total feed: impossible
    with pytest.raises(AdjustError, match="does not change sign"):
        fs.solve()


def test_adjust_component_rate_metric():
    fs = heater_flash()
    fs.logical.append({"type": "adjust", "vary": ["HEAT", "T_out"],
                       "bounds": [320.0, 420.0],
                       "spec": {"type": "component_rate", "stream": "VAP",
                                "component": "n-pentane"},
                       "value": 3.0, "tolerance": 1e-6})
    rep = fs.solve()
    assert rep.converged
    vap = fs.streams["VAP"]
    assert vap.molar_flow * vap.z["n-pentane"] == pytest.approx(3.0, abs=1e-6)


# -- Set -----------------------------------------------------------------------
def test_set_locks_two_params_through_a_solve():
    """FL.T = 1.0*HEAT.T_out - 10 stays locked: the flash runs isothermal 10 K
    below whatever the heater delivers, even as the heater spec changes."""
    fs = heater_flash()
    fs.logical.append({"type": "set", "target": ["FL", "T"],
                       "source": ["HEAT", "T_out"], "multiplier": 1.0, "offset": -10.0})
    fs.solve()
    assert fs.units["FL"].params["T"] == pytest.approx(340.0)
    assert fs.streams["VAP"].T == pytest.approx(340.0)

    fs.units["HEAT"].params["T_out"] = 370.0        # re-solve: the Set re-applies
    fs.solve()
    assert fs.units["FL"].params["T"] == pytest.approx(360.0)
    assert fs.streams["VAP"].T == pytest.approx(360.0)


def test_set_inside_an_adjust_follows_every_iteration():
    """An Adjust varying the *source* of a Set drags the target along (the
    documented order: Sets re-apply before each inner solve)."""
    fs = heater_flash()
    fs.logical.extend([
        {"type": "set", "target": ["FL", "T"], "source": ["HEAT", "T_out"],
         "offset": -5.0},
        dict(VAP_ADJUST),
    ])
    rep = fs.solve()
    assert rep.converged
    assert fs.units["FL"].params["T"] == pytest.approx(
        fs.units["HEAT"].params["T_out"] - 5.0, rel=1e-12)
    assert fs.streams["VAP"].molar_flow == pytest.approx(5.0, abs=1e-6)


def test_set_with_unknown_unit_raises():
    fs = heater_flash()
    fs.logical.append({"type": "set", "target": ["NOPE", "T"],
                       "source": ["HEAT", "T_out"]})
    with pytest.raises(LogicalOpError, match="unknown unit"):
        fs.solve()


def test_unknown_logical_type_raises():
    fs = heater_flash()
    fs.logical.append({"type": "recycle?"})
    with pytest.raises(LogicalOpError, match="unknown logical op type"):
        fs.solve()


# -- `.flow` round-trip ------------------------------------------------------------
def test_logical_and_solver_hints_round_trip_exactly():
    fs = heater_flash()
    fs.logical = [
        {"type": "set", "target": ["FL", "T"], "source": ["HEAT", "T_out"],
         "multiplier": 1.0, "offset": -10.0},
        dict(VAP_ADJUST),
    ]
    fs.solver_hints = {
        "tear_guesses": {"RECY": {"T": 345.0, "P": P_ATM, "molar_flow": 4.2,
                                  "z": {"n-pentane": 0.3, "n-octane": 0.7}}},
        "tear_tolerance": 1e-7,
    }
    doc = to_dict(fs)
    assert doc["logical"] == fs.logical
    assert doc["solver_hints"] == fs.solver_hints
    fs2 = from_dict(doc)
    assert fs2.logical == fs.logical
    assert fs2.solver_hints == fs.solver_hints
    assert to_dict(fs2) == doc                       # byte-exact round-trip

    # Absent by default (no schema noise on plain flowsheets).
    plain = to_dict(heater_flash())
    assert "logical" not in plain and "solver_hints" not in plain


# -- solver hints (Recycle-block equivalent) -----------------------------------------
def test_tear_guess_reduces_iteration_count():
    fs = flash_recycle()
    base = fs.solve(tol=1e-8)
    assert base.converged and base.tear_streams == ["RECY"]
    recy = fs.streams["RECY"]
    guess = {"T": recy.T, "P": recy.P, "molar_flow": recy.molar_flow,
             "z": dict(recy.z)}

    hinted = flash_recycle()
    hinted.solver_hints = {"tear_guesses": {"RECY": guess}}
    rep = hinted.solve(tol=1e-8)
    assert rep.converged
    assert rep.iterations < base.iterations          # a good guess pays off
    assert any("seeded from solver_hints" in m for m in rep.messages)


def test_tear_tolerance_override():
    fs = flash_recycle()
    fs.solver_hints = {"tear_tolerance": 1e-3}
    rep = fs.solve(tol=1e-10)                        # hint wins over the call arg
    assert rep.converged
    assert rep.tol == pytest.approx(1e-3)

    loose, tight = rep.iterations, flash_recycle().solve(tol=1e-10).iterations
    assert loose < tight


# -- balance report -------------------------------------------------------------------
def test_balance_report_closes_on_the_ammonia_loop():
    """Per-unit and overall mass (kg/s) and energy (W, incl. duties) closure on
    the M2 ammonia loop. Units away from the tear close to machine precision;
    the tear-adjacent mixer closes to the solve tolerance."""
    fs = ammonia_loop()
    rep = fs.solve(tol=1e-9, max_iter=600)
    assert rep.converged
    br = balance_report(fs)                          # uses fs.last_report
    assert br["warnings"] == []
    assert br["overall"]["mass_rel"] < 1e-6
    assert br["overall"]["energy_rel"] < 1e-6
    for u in br["units"]:
        assert u["mass_rel"] < 1e-6, u
        assert u["energy_rel"] < 1e-6, u
    # Worst offenders first: the list is sorted by relative imbalance.
    rels = [max(u["mass_rel"], u["energy_rel"]) for u in br["units"]]
    assert rels == sorted(rels, reverse=True)


def test_balance_report_is_machine_precision_on_an_acyclic_flowsheet():
    fs = heater_flash()
    fs.solve()
    br = balance_report(fs)
    assert br["overall"]["mass_rel"] < 1e-12
    assert br["overall"]["energy_rel"] < 1e-12
    for u in br["units"]:
        assert max(u["mass_rel"], u["energy_rel"]) < 1e-12


def test_balance_report_requires_a_solve():
    with pytest.raises(ValueError, match="solve"):
        balance_report(heater_flash())
