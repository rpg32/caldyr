"""M13 tests: specialty distillation — extractive distillation (two-feed
RigorousColumn) and the azeotropic pinch that motivates it.

Validation strategy (sources cited per test):

1. **Extractive distillation** (Hameed, *Chemical Process Simulations using
   Aspen Hysys*, Wiley 2025, sec. 9.5.5): the book separates close-boiling
   n-heptane/toluene by feeding phenol *above* the main feed of a 50-stage
   NRTL column (solvent on stage 4, feed on stage 37), cutting the total
   heating duty from 1.909e7 kJ/h (one 80-stage column, no solvent) to
   7.549e6 kJ/h. **Solvent substitution, documented here**: the bundled
   ChemSep NRTL table has *no* parameters for any of the
   heptane/toluene/phenol pairs (all three would silently fall back to an
   ideal liquid, in which phenol changes no relative volatility and the
   test would prove nothing), so the same flowsheet *structure* is
   demonstrated on acetone/methanol with water as the heavy solvent — a
   classic extractive-distillation system (Seader, Henley & Roper 3e ch. 11;
   Luyben's standard benchmark) for which ChemSep parameterizes all three
   binaries. The book's structural claims are asserted exactly: the same
   column (same stages, reflux, distillate rate) cannot pass the
   acetone/methanol azeotrope without solvent, and feeding the heavy solvent
   a few stages below the condenser lifts the distillate purity well past
   the azeotrope while the solvent leaves in the bottoms (the book's
   Rich-Solvent stream, sent to solvent recovery).
2. **The azeotropic pinch** (Hameed 2025 sec. 9.5.6 problem statement: the
   ethanol/water azeotrope "is a barrier to separation" that forces the
   entrainer flowsheet): a rigorous NRTL column on ethanol/water cannot
   produce a distillate past the azeotrope no matter the stage count or
   reflux — asserted at increasing N and R, including a feed already at
   80 mol% ethanol whose mass balance would happily allow a pure-ethanol
   distillate. Our NRTL puts the azeotrope at ~87.7 mol% ethanol
   (experimental: 89.4 mol% at 1 atm, Gmehling/DECHEMA; the ActivityPackage
   is validated against that point in test_thermo_activity), and the
   distillates pin at it from below to < 0.2 mol%.
"""
import pytest

from caldyr.core import Component, Flowsheet
from caldyr.thermo import make_package
from caldyr.unitops import RigorousColumn

P_ATM = 101325.0


# -- builders -------------------------------------------------------------------
def acetone_methanol_column(solvent: float | None, n=16, R=3.0, D=5.0,
                            feed_stage=9, solvent_stage=4) -> Flowsheet:
    """An equimolar acetone/methanol splitter at 1 atm: 10 mol/s of feed,
    distillate rate = the acetone in the feed. ``solvent`` mol/s of water
    enters above the feed (book sec. 9.5.5 structure: solvent near the top,
    main feed mid-column) when given; None builds the same column without
    the solvent feed."""
    fs = Flowsheet(components=[Component("acetone"), Component("methanol"),
                               Component("water")],
                   property_package="thermo:NRTL")
    params: dict = {"n_stages": n, "reflux_ratio": R, "distillate_rate": D,
                    "P": P_ATM, "max_iter": 600}
    if solvent is None:
        params["feed_stage"] = feed_stage
    else:
        params["feeds"] = [{"stage": feed_stage}, {"stage": solvent_stage}]
    fs.add(RigorousColumn("COL", params))
    fs.feed("FEED1", "COL:in1", T=330.0, P=P_ATM, molar_flow=10.0,
            z={"acetone": 0.5, "methanol": 0.5})
    if solvent is not None:
        fs.feed("SOLV", "COL:in2", T=330.0, P=P_ATM, molar_flow=solvent,
                z={"water": 1.0})
    fs.connect("DIST", "COL:distillate", None)
    fs.connect("BOT", "COL:bottoms", None)
    fs.connect("QC", "COL:condenser_duty", None)
    fs.connect("QR", "COL:reboiler_duty", None)
    return fs


def ethanol_column(n, R, D=4.0, z_etoh=0.5) -> Flowsheet:
    fs = Flowsheet(components=[Component("ethanol"), Component("water")],
                   property_package="thermo:NRTL")
    fs.add(RigorousColumn("COL", {"n_stages": n, "feed_stage": n // 2,
                                  "reflux_ratio": R, "distillate_rate": D,
                                  "P": P_ATM, "max_iter": 600}))
    fs.feed("FEED1", "COL:in1", T=358.0, P=P_ATM, molar_flow=10.0,
            z={"ethanol": z_etoh, "water": 1.0 - z_etoh})
    fs.connect("DIST", "COL:distillate", None)
    fs.connect("BOT", "COL:bottoms", None)
    fs.connect("QC", "COL:condenser_duty", None)
    fs.connect("QR", "COL:reboiler_duty", None)
    return fs


def azeotrope_x(pp, light: str, heavy: str, lo=0.5, hi=0.99) -> float:
    """The minimum-boiling azeotrope composition under the given package:
    the root of y(x) - x by bisection on the bubble-point curve (y > x on
    the dilute side of a positive-deviation azeotrope, y < x above it)."""
    for _ in range(40):
        mid = 0.5 * (lo + hi)
        res = pp.bubble_point(P_ATM, {light: mid, heavy: 1.0 - mid})
        if res.y[light] > mid:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


@pytest.fixture(scope="module")
def solved_extractive() -> Flowsheet:
    """The two-feed extractive column, solved once for the module (the
    azeotropic MESH solve is the expensive part of these tests)."""
    fs = acetone_methanol_column(solvent=10.0)
    assert fs.solve().converged
    return fs


# -- 1. extractive distillation (book sec. 9.5.5 structure) -------------------------
def test_extractive_solvent_above_feed_beats_the_azeotrope(solved_extractive):
    """The book's two-feed structure: without solvent the distillate pins
    below the minimum-boiling azeotrope; the *same* column with the heavy
    solvent fed above the main feed passes it decisively, and the solvent
    leaves in the bottoms (the book's Rich-Solvent stream)."""
    pp = make_package("thermo:NRTL", ["acetone", "methanol"])
    azeo = azeotrope_x(pp, "acetone", "methanol", lo=0.6, hi=0.95)
    # Our ChemSep NRTL puts the acetone/methanol azeotrope at ~79 mol%
    # acetone (experimental ~78-80 mol% at 1 atm; Seader 3e ch. 11).
    assert 0.75 < azeo < 0.85

    base = acetone_methanol_column(solvent=None)
    assert base.solve().converged
    x_base = base.streams["DIST"].z["acetone"]
    assert x_base < azeo + 0.01            # pinned at/below the azeotrope

    extr = solved_extractive
    d, b = extr.streams["DIST"], extr.streams["BOT"]
    # The solvent breaks the pinch: well past the azeotrope and far above
    # the solvent-free column (achieved: 0.71 -> 0.88).
    assert d.z["acetone"] > azeo + 0.05
    assert d.z["acetone"] > x_base + 0.10
    # Methanol is pushed down, and the heavy solvent leaves in the bottoms.
    assert d.z["methanol"] < base.streams["DIST"].z["methanol"] - 0.10
    assert b.z["water"] > 0.5
    assert d.z["water"] < 0.10             # only a little solvent slips up


def test_extractive_two_feed_balances_machine_exact(solved_extractive):
    fs = solved_extractive
    rep = fs.last_report
    d, b = fs.streams["DIST"], fs.streams["BOT"]
    for c in fs.component_ids:
        n_in = sum(fs.streams[sid].molar_flow * fs.streams[sid].z.get(c, 0.0)
                   for sid in ("FEED1", "SOLV"))
        n_out = d.molar_flow * d.z[c] + b.molar_flow * b.z[c]
        assert n_out == pytest.approx(n_in, rel=1e-12, abs=1e-12)
    # Energy closes exactly through the two duty ports.
    e_in = sum(fs.streams[sid].molar_flow * fs.streams[sid].H
               for sid in ("FEED1", "SOLV"))
    e_in += rep.duties["QC"] + rep.duties["QR"]
    e_out = d.molar_flow * d.H + b.molar_flow * b.H
    assert e_out == pytest.approx(e_in, rel=1e-12)


# -- 2. the ethanol/water azeotropic pinch (book sec. 9.5.6 motivation) -------------
def test_azeotrope_pinch_caps_the_distillate():
    """More stages and more reflux push the distillate *toward* the
    azeotrope but never past it — the pinch the book's entrainer flowsheet
    exists to break. D = 4 mol/s with 5 mol/s of ethanol fed, so the mass
    balance would allow a pure-ethanol distillate; equilibrium forbids it."""
    pp = make_package("thermo:NRTL", ["ethanol", "water"])
    azeo = azeotrope_x(pp, "ethanol", "water", lo=0.7, hi=0.99)
    # ~87.7 mol% ethanol under ChemSep NRTL (experimental 89.4 mol%).
    assert 0.85 < azeo < 0.92

    x_d = []
    for n, r in ((16, 2.0), (22, 4.0)):
        fs = ethanol_column(n, r)
        assert fs.solve().converged
        x_d.append(fs.streams["DIST"].z["ethanol"])
    assert x_d[0] < x_d[1]                 # more stages/reflux: closer...
    for x in x_d:
        assert x < azeo + 0.002            # ...but never past the azeotrope


def test_azeotrope_pinch_holds_for_near_azeotropic_feed():
    """Even feeding 80 mol% ethanol (8 mol/s ethanol vs D = 4 mol/s, so the
    balance is no constraint at all), the distillate stays below the
    azeotrope."""
    pp = make_package("thermo:NRTL", ["ethanol", "water"])
    azeo = azeotrope_x(pp, "ethanol", "water", lo=0.7, hi=0.99)
    fs = ethanol_column(20, 3.0, z_etoh=0.8)
    assert fs.solve().converged
    x = fs.streams["DIST"].z["ethanol"]
    assert 0.8 < x < azeo + 0.002
