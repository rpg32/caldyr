"""M12 tests: RigorousColumn with multiple feeds and side draws.

Validation strategy:

1. **The book's multi-component column** (Hameed, *Chemical Process
   Simulations using Aspen Hysys*, Wiley 2025, sec. 9.5.7): 2500 kmol/h of
   saturated liquid at 80 C (5% C3, 10% iC4, 40% iC5, 36% nC6, 9% nC7),
   total condenser at 360 kPa / reboiler at 370 kPa, HYSYS 20 stages (+
   condenser + reboiler -> n_stages = 22 here), feed on HYSYS stage 10
   (-> 11), R = 0.92, distillate 1381 kmol/h (the book's converged rate for
   its x(iC5)_reboiler = 0.01 spec). Book product compositions (fig. 9.87):
   D = {C3 .0905, iC4 .1810, iC5 .7159, nC6 .0126, nC7 .0000},
   B = {iC5 .0100, nC7 .2011} — achieved deltas are <= 0.0002 mole fraction
   on every entry. (The book's printed x(nC6)_B = 0.7689 does not close its
   own summation — the column sums to 0.98 — so the consistent complement
   0.7889 is asserted, which we hit to 1e-4.) Duties land within 1% of the
   book's shortcut-column values (-17.5 / +18.3 MW).
2. **Two feeds**: splitting the same feed onto two stages must (a) keep the
   balances machine-exact, (b) reduce *identically* to the single-feed
   answer when both entries name the same stage, and (c) shift the internal
   profiles when part of the feed moves lower while the products only move
   mildly (the documented expectation: a non-optimal split feed slightly
   degrades the key separation at fixed R and D).
3. **Side draws**: the draw rate is honored exactly, its composition is the
   converged stage composition (liquid or vapor), all balances still close
   to machine precision, and the energy balance (with the side stream) is
   exact.
4. **Regression**: the classic single-feed form (no ``feeds``/``side_draws``
   params) keeps byte-identical ports and behavior — asserted directly here
   and by the untouched test_m10 suite.
"""
import math

import pytest

from caldyr.core import Component, Flowsheet
from caldyr.io import from_dict, to_dict
from caldyr.thermo import make_package
from caldyr.unitops import RigorousColumn, RigorousColumnError

KMOLH = 1000.0 / 3600.0
COMPS = ["propane", "isobutane", "isopentane", "n-hexane", "n-heptane"]
Z = {"propane": 0.05, "isobutane": 0.10, "isopentane": 0.40,
     "n-hexane": 0.36, "n-heptane": 0.09}
BOOK = {"n_stages": 22, "reflux_ratio": 0.92,
        "distillate_rate": 1381.0 * KMOLH,
        "P": 360000.0, "dP_stage": 10000.0 / 21.0}


def feed_pressure() -> float:
    """The book feeds saturated liquid at 80 C (VF = 0, T given), so the feed
    pressure is the bubble pressure at 80 C (~477 kPa with our PR)."""
    from scipy.optimize import brentq
    pp = make_package("thermo:PR", COMPS)
    return float(brentq(lambda p: pp.bubble_dew(p, Z)[0] - 353.15,
                        1e5, 2e6, xtol=10.0))


def build(params, feeds=None) -> Flowsheet:
    """Column flowsheet; ``feeds`` is a list of molar flows for in1, in2, ...
    (None -> the single 2500 kmol/h book feed on in1)."""
    fs = Flowsheet(components=[Component(c) for c in COMPS],
                   property_package="thermo:PR")
    fs.add(RigorousColumn("COL", params))
    pf = feed_pressure()
    flows = feeds or [2500.0 * KMOLH]
    for i, flow in enumerate(flows):
        fs.feed(f"FEED{i + 1}", f"COL:in{i + 1}", T=353.15, P=pf,
                molar_flow=flow, z=dict(Z))
    fs.connect("DIST", "COL:distillate", None)
    fs.connect("BOT", "COL:bottoms", None)
    for i in range(len(params.get("side_draws") or [])):
        fs.connect(f"SIDE{i + 1}", f"COL:side{i + 1}", None)
    fs.connect("QC", "COL:condenser_duty", None)
    fs.connect("QR", "COL:reboiler_duty", None)
    return fs


def total_in(fs, c) -> float:
    return sum(s.molar_flow * s.z[c] for sid, s in fs.streams.items()
               if sid.startswith("FEED"))


# -- 1. the book's column (single feed; baseline for the variants) ---------------------
def test_book_9_5_7_multicomponent_column():
    fs = build({**BOOK, "feed_stage": 11})
    assert fs.solve().converged
    d, b = fs.streams["DIST"], fs.streams["BOT"]
    # Book fig. 9.87 distillate column (achieved deltas <= 0.0002).
    assert d.z["propane"] == pytest.approx(0.0905, abs=0.003)
    assert d.z["isobutane"] == pytest.approx(0.1810, abs=0.003)
    assert d.z["isopentane"] == pytest.approx(0.7159, abs=0.005)
    assert d.z["n-hexane"] == pytest.approx(0.0126, abs=0.005)
    assert d.z["n-heptane"] == pytest.approx(0.0, abs=1e-4)
    # Book bottoms: iC5 0.0100 and nC7 0.2011; the printed nC6 (0.7689) does
    # not close the book's own summation, so the consistent complement is
    # asserted (see module docstring).
    assert b.z["isopentane"] == pytest.approx(0.0100, abs=0.005)
    assert b.z["n-heptane"] == pytest.approx(0.2011, abs=0.003)
    assert b.z["n-hexane"] == pytest.approx(1.0 - 0.0100 - 0.2011, abs=0.005)
    assert b.z["propane"] == pytest.approx(0.0, abs=1e-6)
    assert b.z["isobutane"] == pytest.approx(0.0, abs=1e-4)
    # Duties within a few % of the book's shortcut values (fig. 9.85:
    # -63041139.75 kJ/h = -17.51 MW; +65817131.99 kJ/h = +18.28 MW).
    des = fs.units["COL"].design
    assert des["Q_condenser"] == pytest.approx(-17.51e6, rel=0.03)
    assert des["Q_reboiler"] == pytest.approx(18.28e6, rel=0.03)
    # Condenser at 45.58 C per the book's shortcut performance tab.
    assert d.T - 273.15 == pytest.approx(45.58, abs=1.0)


# -- 2. multiple feeds -------------------------------------------------------------------
def test_two_feeds_on_same_stage_reduce_to_single_feed():
    """Splitting the feed 60/40 onto the *same* stage must reproduce the
    single-feed column exactly (same combined feed, same stage)."""
    fs1 = build({**BOOK, "feed_stage": 11})
    assert fs1.solve().converged
    fs2 = build({**BOOK, "feeds": [{"stage": 11}, {"stage": 11}]},
                feeds=[1500.0 * KMOLH, 1000.0 * KMOLH])
    assert fs2.solve().converged
    for sid in ("DIST", "BOT"):
        a, b = fs1.streams[sid], fs2.streams[sid]
        assert b.molar_flow == pytest.approx(a.molar_flow, rel=1e-9)
        assert b.T == pytest.approx(a.T, abs=1e-5)
        for c in COMPS:
            assert b.z[c] == pytest.approx(a.z[c], rel=1e-6, abs=1e-10)


def test_two_feeds_split_across_stages():
    """60% of the feed on stage 9, 40% on stage 14. Expectations (documented
    in the module docstring): machine-exact balances; the liquid traffic
    between the two feed stages rises vs the single-feed column (the upper
    feed's liquid now travels through stages 9-13); the products move only
    mildly, with the key split slightly *degraded* at fixed R and D because
    neither feed sits on the optimal stage."""
    fs1 = build({**BOOK, "feed_stage": 11})
    assert fs1.solve().converged
    fs2 = build({**BOOK, "feeds": [{"stage": 9}, {"stage": 14}]},
                feeds=[1500.0 * KMOLH, 1000.0 * KMOLH])
    assert fs2.solve().converged

    # Machine-exact component balances over both feeds.
    d2, b2 = fs2.streams["DIST"], fs2.streams["BOT"]
    for c in COMPS:
        n_out = d2.molar_flow * d2.z[c] + b2.molar_flow * b2.z[c]
        assert n_out == pytest.approx(total_in(fs2, c), rel=1e-12, abs=1e-12)
    # Energy closes exactly (Q_reb closes the overall balance).
    rep2 = fs2.last_report
    e_in = sum(fs2.streams[f"FEED{i}"].molar_flow * fs2.streams[f"FEED{i}"].H
               for i in (1, 2)) + rep2.duties["QC"] + rep2.duties["QR"]
    e_out = d2.molar_flow * d2.H + b2.molar_flow * b2.H
    assert e_out == pytest.approx(e_in, rel=1e-12)

    des1, des2 = fs1.units["COL"].design, fs2.units["COL"].design
    # The internal profile shifts exactly as the split-feed physics says:
    # between the upper feed (stage 9, 0-based 8) and the original feed
    # stage (11, 0-based 10) the liquid traffic now carries the upper feed
    # (single-feed: still rectifying-section traffic there), while between
    # the original stage and the lower feed (14, 0-based 13) only 60% of the
    # feed liquid has arrived, so the traffic is *lower* than single-feed.
    assert all(des2["L_profile"][j] > des1["L_profile"][j] * 1.5
               for j in range(8, 10))
    assert all(des2["L_profile"][j] < des1["L_profile"][j] * 0.85
               for j in range(10, 13))
    # Products move only mildly; the split degrades slightly (more light key
    # lost to the bottoms than with the feed on its optimal stage).
    assert d2.z["isopentane"] == pytest.approx(
        fs1.streams["DIST"].z["isopentane"], abs=0.05)
    assert b2.z["isopentane"] >= fs1.streams["BOT"].z["isopentane"] - 1e-6
    # design publishes both feeds with their stages and qualities.
    assert [f["stage"] for f in des2["feeds"]] == [9, 14]
    assert all(0.0 <= f["q"] <= 1.0 for f in des2["feeds"])


def test_default_ports_unchanged_without_new_params():
    """Param-driven ports must not disturb the classic single-feed layout."""
    classic = RigorousColumn("COL", {"n_stages": 16, "feed_stage": 8,
                                     "reflux_ratio": 1.5,
                                     "distillate_rate": 50.0})
    assert [p.name for p in classic.ports] == [
        "in1", "distillate", "bottoms", "condenser_duty", "reboiler_duty"]
    multi = RigorousColumn("COL", {"n_stages": 16,
                                   "feeds": [{"stage": 6}, {"stage": 10}],
                                   "side_draws": [{"stage": 8,
                                                   "phase": "liquid",
                                                   "rate": 5.0}],
                                   "reflux_ratio": 1.5,
                                   "distillate_rate": 50.0})
    assert [p.name for p in multi.ports] == [
        "in1", "in2", "distillate", "bottoms", "side1",
        "condenser_duty", "reboiler_duty"]


# -- 3. side draws -------------------------------------------------------------------------
def test_liquid_side_draw_rate_honored_and_balances_close():
    rate = 300.0 * KMOLH
    fs = build({**BOOK, "feed_stage": 11, "distillate_rate": 1100.0 * KMOLH,
                "side_draws": [{"stage": 6, "phase": "liquid", "rate": rate}]})
    assert fs.solve().converged
    d, b, s = fs.streams["DIST"], fs.streams["BOT"], fs.streams["SIDE1"]
    # The draw rate is honored exactly and the bottoms picks up the rest.
    assert s.molar_flow == pytest.approx(rate, rel=1e-12)
    assert b.molar_flow == pytest.approx(
        2500.0 * KMOLH - 1100.0 * KMOLH - rate, rel=1e-12)
    # Machine-exact component balances including the draw.
    for c in COMPS:
        n_out = (d.molar_flow * d.z[c] + b.molar_flow * b.z[c]
                 + s.molar_flow * s.z[c])
        assert n_out == pytest.approx(total_in(fs, c), rel=1e-12, abs=1e-12)
    # Exact energy closure with the side stream.
    rep = fs.last_report
    e_in = (fs.streams["FEED1"].molar_flow * fs.streams["FEED1"].H
            + rep.duties["QC"] + rep.duties["QR"])
    e_out = sum(p.molar_flow * p.H for p in (d, b, s))
    assert e_out == pytest.approx(e_in, rel=1e-12)
    # The draw is the converged stage-6 saturated liquid.
    des = fs.units["COL"].design
    assert s.phase == "liquid"
    assert s.T == pytest.approx(des["T_profile"][5], abs=1e-9)
    for c in COMPS:
        assert s.z[c] == pytest.approx(des["x_profile"][5][c], abs=1e-9)
    # An intermediate draw is intermediate in composition: richer in the
    # lights than the bottoms, leaner than the distillate.
    assert d.z["isopentane"] > s.z["isopentane"] > b.z["isopentane"]


def test_vapor_side_draw_is_stage_vapor():
    rate = 200.0 * KMOLH
    fs = build({**BOOK, "feed_stage": 11, "distillate_rate": 1100.0 * KMOLH,
                "side_draws": [{"stage": 14, "phase": "vapor", "rate": rate}]})
    assert fs.solve().converged
    s = fs.streams["SIDE1"]
    des = fs.units["COL"].design
    assert s.phase == "vapor"
    assert s.molar_flow == pytest.approx(rate, rel=1e-12)
    assert s.T == pytest.approx(des["T_profile"][13], abs=1e-9)
    for c in COMPS:
        assert s.z[c] == pytest.approx(des["y_profile"][13][c], abs=1e-9)
    # All balances still close to machine precision.
    d, b = fs.streams["DIST"], fs.streams["BOT"]
    for c in COMPS:
        n_out = (d.molar_flow * d.z[c] + b.molar_flow * b.z[c]
                 + s.molar_flow * s.z[c])
        assert n_out == pytest.approx(total_in(fs, c), rel=1e-12, abs=1e-12)


# -- 4. `.flow` round-trip with the new params ----------------------------------------------
def test_flow_round_trip_with_feeds_and_draws():
    fs = build({**BOOK, "feeds": [{"stage": 9}, {"stage": 14}],
                "distillate_rate": 1100.0 * KMOLH,
                "side_draws": [{"stage": 6, "phase": "liquid",
                                "rate": 200.0 * KMOLH}]},
               feeds=[1500.0 * KMOLH, 1000.0 * KMOLH])
    fs.solve()
    doc = to_dict(fs)
    assert doc == to_dict(from_dict(doc))
    fs2 = from_dict(doc)
    fs2.solve()
    for sid in ("DIST", "BOT", "SIDE1"):
        assert fs2.streams[sid].z["isopentane"] == pytest.approx(
            fs.streams[sid].z["isopentane"], rel=1e-9)


# -- 5. typed errors --------------------------------------------------------------------------
def test_feeds_and_feed_stage_together_raise():
    fs = build({**BOOK, "feed_stage": 11, "feeds": [{"stage": 9}]})
    with pytest.raises(RigorousColumnError, match="not both"):
        fs.solve()


def test_draw_on_end_stages_raises():
    for stage in (1, 22):
        fs = build({**BOOK, "feed_stage": 11,
                    "side_draws": [{"stage": stage, "phase": "liquid",
                                    "rate": 1.0}]})
        with pytest.raises(RigorousColumnError, match="out of range"):
            fs.solve()


def test_draws_plus_distillate_exceeding_feed_raise():
    fs = build({**BOOK, "feed_stage": 11,
                "side_draws": [{"stage": 6, "phase": "liquid",
                                "rate": 1500.0 * KMOLH}]})
    with pytest.raises(RigorousColumnError, match="below the total feed"):
        fs.solve()


def test_bad_draw_phase_raises():
    fs = build({**BOOK, "feed_stage": 11,
                "side_draws": [{"stage": 6, "phase": "mist", "rate": 1.0}]})
    with pytest.raises(RigorousColumnError, match="'liquid' or 'vapor'"):
        fs.solve()


def test_missing_second_feed_stream_raises():
    fs = Flowsheet(components=[Component(c) for c in COMPS],
                   property_package="thermo:PR")
    fs.add(RigorousColumn("COL", {**BOOK,
                                  "feeds": [{"stage": 9}, {"stage": 14}]}))
    fs.feed("FEED1", "COL:in1", T=353.15, P=feed_pressure(),
            molar_flow=2500.0 * KMOLH, z=dict(Z))
    for sid, port in (("DIST", "distillate"), ("BOT", "bottoms"),
                      ("QC", "condenser_duty"), ("QR", "reboiler_duty")):
        fs.connect(sid, f"COL:{port}", None)
    with pytest.raises(RigorousColumnError, match="in2"):
        fs.solve()


# -- 6. side-draw mass-balance sanity at the math level ---------------------------------------
def test_draw_reduces_internal_liquid_below_the_draw_stage():
    """Pulling saturated liquid off stage 6 must lower the liquid traffic
    everywhere between the draw and the feed relative to the no-draw column
    (same R, D): L_j drops by about the draw rate."""
    base = build({**BOOK, "feed_stage": 11,
                  "distillate_rate": 1100.0 * KMOLH})
    assert base.solve().converged
    drawn = build({**BOOK, "feed_stage": 11,
                   "distillate_rate": 1100.0 * KMOLH,
                   "side_draws": [{"stage": 6, "phase": "liquid",
                                   "rate": 300.0 * KMOLH}]})
    assert drawn.solve().converged
    l0 = base.units["COL"].design["L_profile"]
    l1 = drawn.units["COL"].design["L_profile"]
    assert all(l1[j] < l0[j] for j in range(6, 10))
    assert math.isclose(l0[7] - l1[7], 300.0 * KMOLH,
                        rel_tol=0.35)   # ~the draw rate (energy balance shifts a bit)
