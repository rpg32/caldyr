"""M4 acceptance tests: the equation-oriented backend and optimization.

The cross-check that matters: the sequential-modular and equation-oriented
backends share the same unit physics and property package, so they must converge
to the same flowsheet state — on acyclic graphs, on recycles (which the EO solver
handles with no tear stream), and with a reactor in the loop.
"""
import math

import pytest

from caldyr.core import Component, Flowsheet
from caldyr.solver import DesignVar, optimize
from caldyr.unitops import ConversionReactor, FlashDrum, Heater, Mixer, Splitter

AMMONIA = {"stoich": {"nitrogen": -1, "hydrogen": -3, "ammonia": 2}, "key": "nitrogen"}


def _streams_close(a: Flowsheet, b: Flowsheet, ids, rel=1e-6):
    for sid in ids:
        sa, sb = a.streams[sid], b.streams[sid]
        assert math.isclose(sa.molar_flow, sb.molar_flow, rel_tol=rel, abs_tol=1e-9)
        assert math.isclose(sa.T, sb.T, rel_tol=rel, abs_tol=1e-6)
        for c in sa.components:
            assert math.isclose(sa.z.get(c, 0.0), sb.z.get(c, 0.0), rel_tol=rel, abs_tol=1e-7)


def _solve_both(build, ids, duty_ids=()):
    fa = build()
    ra = fa.solve(backend="sequential", tol=1e-9)
    fb = build()
    rb = fb.solve(backend="equation_oriented", tol=1e-10)
    assert ra.converged and rb.converged
    _streams_close(fa, fb, ids)
    for q in duty_ids:
        assert math.isclose(ra.duties[q], rb.duties[q], rel_tol=1e-5, abs_tol=1.0)
    return fa, fb, ra, rb


def mixer_heater() -> Flowsheet:
    fs = Flowsheet(components=[Component("water"), Component("ethanol")],
                   property_package="thermo:PR")
    fs.add(Mixer("MIX", {"dP": 0.0}))
    fs.add(Heater("H", {"T_out": 350.0, "dP": 0.0}))
    fs.feed("S1", "MIX:in1", T=298.15, P=101325.0, molar_flow=10.0,
            z={"water": 0.6, "ethanol": 0.4})
    fs.feed("S2", "MIX:in2", T=320.0, P=101325.0, molar_flow=5.0, z={"water": 1.0})
    fs.connect("S3", "MIX:out", "H:in1")
    fs.connect("S4", "H:out", None)
    fs.connect("Q", "H:duty", None)
    return fs


def flash_recycle() -> Flowsheet:
    fs = Flowsheet(components=[Component("n-pentane"), Component("n-octane")],
                   property_package="thermo:PR")
    fs.add(Mixer("MIX", {"dP": 0.0}))
    fs.add(FlashDrum("FL", {"T": 360.0, "P": 101325.0}))
    fs.add(Splitter("SP", {"split": 0.6}))
    fs.feed("FEED", "MIX:in1", T=330.0, P=101325.0, molar_flow=10.0,
            z={"n-pentane": 0.5, "n-octane": 0.5})
    fs.connect("MIXOUT", "MIX:out", "FL:in1")
    fs.connect("VAP", "FL:vapor", None)
    fs.connect("LIQ", "FL:liquid", "SP:in1")
    fs.connect("RECY", "SP:out1", "MIX:in2")
    fs.connect("BOT", "SP:out2", None)
    fs.connect("Q", "FL:duty", None)
    return fs


def reactor_once_through() -> Flowsheet:
    fs = Flowsheet(
        components=[Component("nitrogen"), Component("hydrogen"), Component("ammonia")],
        property_package="thermo:PR")
    fs.add(ConversionReactor("R", {"reaction": AMMONIA, "conversion": 0.25, "T_out": 700.0}))
    fs.feed("F", "R:in1", T=673.15, P=2e7, molar_flow=4.0,
            z={"nitrogen": 0.25, "hydrogen": 0.75, "ammonia": 0.0})
    fs.connect("OUT", "R:out", None)
    fs.connect("Q", "R:duty", None)
    return fs


def test_eo_matches_sm_acyclic():
    _solve_both(mixer_heater, ["S3", "S4"], duty_ids=["Q"])


def test_eo_matches_sm_with_recycle_no_tear():
    fa, fb, ra, rb = _solve_both(flash_recycle, ["VAP", "BOT", "RECY"], duty_ids=["Q"])
    assert rb.tear_streams == []                 # EO needs no tear stream
    assert rb.method == "equation_oriented"


def test_eo_matches_sm_with_reactor():
    fa, fb, ra, rb = _solve_both(reactor_once_through, ["OUT"], duty_ids=["Q"])
    # Reaction physics carried through both backends identically.
    assert math.isclose(fa.streams["OUT"].z["ammonia"],
                        fb.streams["OUT"].z["ammonia"], rel_tol=1e-6)


def test_eo_closes_mass_balance_on_recycle():
    fs = flash_recycle()
    fs.solve(backend="equation_oriented", tol=1e-10)
    s = fs.streams
    n_in = s["FEED"].molar_flow
    n_out = s["VAP"].molar_flow + s["BOT"].molar_flow
    assert math.isclose(n_in, n_out, rel_tol=1e-6)


# -- optimization ----------------------------------------------------------
def test_optimization_meets_recovery_spec_at_min_duty():
    fs = flash_recycle()

    def overhead(fs):
        v = fs.streams["VAP"]
        return v.molar_flow * v.z["n-pentane"]

    res = optimize(
        fs,
        objective=lambda fs, rep: rep.duties["Q"] / 1e3,    # kW (well-scaled)
        design_vars=[DesignVar("FL", "T", 340.0, 370.0, initial=360.0)],
        constraints=[lambda fs, rep: overhead(fs) - 4.2],
        solve_kwargs={"tol": 1e-9},
    )
    assert res.success
    assert 340.0 <= res.design["FL.T"] <= 370.0
    # The recovery constraint is active at the optimum (minimizing duty pushes T
    # down, the spec pushes it up) -> overhead sits at the 4.2 mol/s bound.
    assert overhead(fs) == pytest.approx(4.2, abs=2e-3)


def test_optimization_respects_bounds_when_unconstrained():
    """Minimizing duty with no recovery constraint drives the flash to its lower
    temperature bound (least vaporization)."""
    fs = flash_recycle()
    res = optimize(
        fs,
        objective=lambda fs, rep: rep.duties["Q"] / 1e3,
        design_vars=[DesignVar("FL", "T", 345.0, 365.0, initial=360.0)],
        solve_kwargs={"tol": 1e-9},
    )
    assert res.success
    assert res.design["FL.T"] == pytest.approx(345.0, abs=0.5)
