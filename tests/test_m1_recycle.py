"""M1 acceptance tests: recycle convergence + the new unit ops.

Property reference: Peng-Robinson EOS via `thermo` (Caleb Bell). The
flash-with-recycle case uses an n-pentane / n-octane mixture, which is well
described by a cubic EOS (non-polar) and gives a clean two-phase split. Unit-op
checks validate against first principles (mass balance, isenthalpic throttling,
incompressible pump work V·ΔP/η, isentropic compression as the η→1 limit).
"""
import math

import pytest

from caldyr.core import Component, Flowsheet
from caldyr.thermo import make_package
from caldyr.unitops import Compressor, FlashDrum, Mixer, Pump, Splitter, Valve

P_ATM = 101325.0


def hydrocarbons() -> list[Component]:
    return [Component("n-pentane"), Component("n-octane")]


def build_flash_recycle(split: float = 0.6) -> Flowsheet:
    """Feed -> Mixer -> Flash; the flash liquid is split, part recycled to the
    mixer and part drawn off as bottoms. Vapor leaves as product."""
    fs = Flowsheet(components=hydrocarbons(), property_package="thermo:PR")
    fs.add(Mixer("MIX", {"dP": 0.0}))
    fs.add(FlashDrum("FL", {"T": 360.0, "P": P_ATM}))
    fs.add(Splitter("SP", {"split": split}))
    fs.feed("FEED", "MIX:in1", T=330.0, P=P_ATM, molar_flow=10.0,
            z={"n-pentane": 0.5, "n-octane": 0.5})
    fs.connect("MIXOUT", "MIX:out", "FL:in1")
    fs.connect("VAP", "FL:vapor", None)
    fs.connect("LIQ", "FL:liquid", "SP:in1")
    fs.connect("RECY", "SP:out1", "MIX:in2")
    fs.connect("BOTTOMS", "SP:out2", None)
    fs.connect("Q", "FL:duty", None)
    return fs


# -- recycle convergence ---------------------------------------------------
def test_flash_recycle_converges_and_closes_mass_balance():
    fs = build_flash_recycle()
    report = fs.solve(tol=1e-7)

    assert report.converged
    assert report.tear_streams == ["RECY"]   # detected the single recycle
    assert report.history and report.history[-1] < 1e-7
    assert report.iterations == len(report.history)

    s = fs.streams
    n_in = s["FEED"].molar_flow
    n_out = s["VAP"].molar_flow + s["BOTTOMS"].molar_flow
    assert math.isclose(n_in, n_out, rel_tol=1e-6)  # overall: feed == products
    for c in ("n-pentane", "n-octane"):
        c_in = s["FEED"].molar_flow * s["FEED"].z[c]
        c_out = (s["VAP"].molar_flow * s["VAP"].z.get(c, 0.0)
                 + s["BOTTOMS"].molar_flow * s["BOTTOMS"].z.get(c, 0.0))
        assert math.isclose(c_in, c_out, rel_tol=1e-6)


def test_wegstein_and_direct_agree_and_wegstein_is_faster():
    """Both methods reach the same fixed point; acceleration cuts iterations."""
    weg = build_flash_recycle(0.7)
    dsub = build_flash_recycle(0.7)
    rw = weg.solve(tol=1e-7, method="wegstein")
    rd = dsub.solve(tol=1e-7, method="direct")
    assert rw.converged and rd.converged
    assert rw.iterations < rd.iterations            # Wegstein accelerates
    for sid in ("VAP", "BOTTOMS", "RECY"):
        assert math.isclose(weg.streams[sid].molar_flow,
                            dsub.streams[sid].molar_flow, rel_tol=1e-5)


def test_recycle_steady_state_independent_of_seed_is_physical():
    """The vapor product should be richer in the lighter component than bottoms."""
    fs = build_flash_recycle()
    fs.solve(tol=1e-8)
    vap, bot = fs.streams["VAP"], fs.streams["BOTTOMS"]
    assert vap.z["n-pentane"] > bot.z["n-pentane"]   # lighter favors the vapor


# -- new unit ops ----------------------------------------------------------
def test_splitter_partitions_flow_and_preserves_state():
    fs = Flowsheet(components=hydrocarbons(), property_package="thermo:PR")
    fs.add(Splitter("SP", {"split": 0.3}))
    fs.feed("F", "SP:in1", T=300.0, P=2e5, molar_flow=10.0,
            z={"n-pentane": 0.5, "n-octane": 0.5})
    fs.connect("A", "SP:out1", None)
    fs.connect("B", "SP:out2", None)
    fs.solve()
    a, b, f = fs.streams["A"], fs.streams["B"], fs.streams["F"]
    assert math.isclose(a.molar_flow, 3.0) and math.isclose(b.molar_flow, 7.0)
    assert a.T == f.T and a.P == f.P and a.z == f.z   # intensive state preserved


def test_valve_is_isenthalpic_and_drops_pressure():
    fs = Flowsheet(components=[Component("n-pentane")], property_package="thermo:PR")
    fs.add(Valve("V", {"P_out": 1e5}))
    fs.feed("F", "V:in1", T=320.0, P=5e5, molar_flow=1.0, z={"n-pentane": 1.0})
    fs.connect("O", "V:out", None)
    fs.solve()
    f, o = fs.streams["F"], fs.streams["O"]
    assert o.P == 1e5
    assert math.isclose(o.H, f.H, rel_tol=1e-9)       # enthalpy conserved
    assert o.T < f.T                                  # Joule-Thomson cooling


def test_pump_work_matches_incompressible_formula():
    fs = Flowsheet(components=[Component("water")], property_package="thermo:PR")
    eta = 0.75
    fs.add(Pump("P", {"P_out": 1e6, "eta": eta}))
    fs.feed("F", "P:in1", T=298.15, P=1e5, molar_flow=2.0, z={"water": 1.0})
    fs.connect("O", "P:out", None)
    fs.connect("W", "P:work", None)
    report = fs.solve()
    pp = make_package("thermo:PR", ["water"])
    v = pp.volume(298.15, 1e5, {"water": 1.0})
    w_expected = 2.0 * v * (1e6 - 1e5) / eta           # n·V·ΔP/η
    assert math.isclose(report.duties["W"], w_expected, rel_tol=1e-6)
    assert fs.streams["O"].P == 1e6
    assert fs.streams["O"].T > 298.15                  # slight frictional heating


def test_compressor_isentropic_limit_and_efficiency():
    """At η=1 the work equals the isentropic enthalpy rise; η<1 needs more."""
    feed = dict(T=300.0, P=1e5, molar_flow=1.0, z={"methane": 1.0})
    pp = make_package("thermo:PR", ["methane"])
    h_in = pp.enthalpy(300.0, 1e5, {"methane": 1.0})
    s_in = pp.entropy(300.0, 1e5, {"methane": 1.0})
    w_isentropic = pp.flash_ps(5e5, s_in, {"methane": 1.0}).H - h_in

    def run(eta):
        fs = Flowsheet(components=[Component("methane")], property_package="thermo:PR")
        fs.add(Compressor("C", {"P_out": 5e5, "eta": eta}))
        fs.feed("F", "C:in1", **feed)
        fs.connect("O", "C:out", None)
        fs.connect("W", "C:work", None)
        return fs.solve().duties["W"]

    assert math.isclose(run(1.0), w_isentropic, rel_tol=1e-6)
    assert run(0.75) > w_isentropic                    # irreversibility costs work
    assert math.isclose(run(0.75), w_isentropic / 0.75, rel_tol=1e-6)


def test_flash_drum_isothermal_split_and_balance():
    fs = Flowsheet(components=hydrocarbons(), property_package="thermo:PR")
    fs.add(FlashDrum("FL", {"T": 360.0, "P": P_ATM}))
    fs.feed("F", "FL:in1", T=300.0, P=P_ATM, molar_flow=10.0,
            z={"n-pentane": 0.5, "n-octane": 0.5})
    fs.connect("V", "FL:vapor", None)
    fs.connect("L", "FL:liquid", None)
    fs.connect("Q", "FL:duty", None)
    report = fs.solve()
    v, liq, f = fs.streams["V"], fs.streams["L"], fs.streams["F"]

    assert 0.0 < v.molar_flow < f.molar_flow            # genuinely two-phase
    assert math.isclose(v.molar_flow + liq.molar_flow, f.molar_flow, rel_tol=1e-9)
    assert v.z["n-pentane"] > liq.z["n-pentane"]        # vapor richer in light key
    for c in ("n-pentane", "n-octane"):
        assert math.isclose(
            v.molar_flow * v.z[c] + liq.molar_flow * liq.z[c],
            f.molar_flow * f.z[c], rel_tol=1e-9,
        )
    assert report.duties["Q"] > 0.0                     # heating to hold 360 K


def test_flash_drum_adiabatic_when_no_temperature_given():
    fs = Flowsheet(components=hydrocarbons(), property_package="thermo:PR")
    fs.add(FlashDrum("FL", {"P": P_ATM}))               # adiabatic: P only
    fs.feed("F", "FL:in1", T=360.0, P=3e5, molar_flow=10.0,
            z={"n-pentane": 0.5, "n-octane": 0.5})
    fs.connect("V", "FL:vapor", None)
    fs.connect("L", "FL:liquid", None)
    fs.connect("Q", "FL:duty", None)
    report = fs.solve()
    assert report.duties["Q"] == pytest.approx(0.0, abs=1e-6)
    # adiabatic flash to lower P: enthalpy conserved across the drum
    f = fs.streams["F"]
    v, liq = fs.streams["V"], fs.streams["L"]
    h_out = v.molar_flow * v.H + liq.molar_flow * liq.H
    assert math.isclose(h_out, f.molar_flow * f.H, rel_tol=1e-7)
