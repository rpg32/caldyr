"""M18 tests: reactive distillation on RigorousColumn (book §9.5.3).

A reactive distillation column hosts a chemical reaction on a band of trays while
the column simultaneously separates the products — continuously pulling a product
off shifts a reversible equilibrium toward more conversion (Le Chatelier), the
whole point of the technique. The engine adds kinetic reactions (forward and
reversible, see :class:`caldyr.unitops.reaction.KineticReaction`) on named stages
of the bubble-point MESH: the per-stage generation enters the component balances
as an extra source, and the heat of reaction is carried automatically by the
formation-inclusive stage enthalpies.

Validation references:

* **Methyl-acetate synthesis** (Hameed §9.5.3 — methanol + acetic acid ⇌ methyl
  acetate + water, the book's reactive-distillation worked example). The book
  uses HYSYS with the Wilson fluid package and a kinetic reaction; caldyr's
  activity-model analogue is ``thermo:NRTL`` (caldyr has no Wilson). The book
  prints no stream numbers, so the checks are first-principles: the reaction
  consumes methanol/acetic acid and makes methyl acetate, the light ester
  enriches in the distillate, element balances close exactly, and conversion
  rises with the tray holdup (vanishing as holdup → 0). **Documented thermo
  limit:** pushing to high conversion drives the trays into the acetic-acid
  vapour-association regime where NRTL returns h_V < h_L (and PR into the
  ester/water immiscibility region) — so the validated regime is the moderate
  conversion the explicit reactive-MESH reaches; the equilibrium-limited
  fast-kinetics regime needs a simultaneous reactive solver (follow-up).

* **Toluene disproportionation** (2 toluene ⇌ benzene + xylene) — a real reactive
  separation among fully-miscible aromatics that ``thermo:PR`` handles cleanly, a
  robustness cross-check that reaches ~60 % conversion with benzene taken
  overhead and xylene as bottoms. The mole-conserving stoichiometry keeps the
  carbon/hydrogen balances machine-exact.
"""
import pytest

from caldyr.core import Component, Flowsheet
from caldyr.io import from_dict, to_dict
from caldyr.unitops import RigorousColumn, RigorousColumnError
from caldyr.unitops.reaction import KineticReaction

# ---- methyl acetate (book §9.5.3) -------------------------------------------
MA_COMPS = ["methanol", "acetic acid", "methyl acetate", "water"]
MA_RXN = {
    "stoich": {"methanol": -1, "acetic acid": -1,
               "methyl acetate": 1, "water": 1},
    "key": "methanol", "k0": 3e-8, "Ea": 0.0,
    "orders": {"methanol": 1, "acetic acid": 1},
    "k0_rev": 6e-9, "Ea_rev": 0.0,
    "orders_rev": {"methyl acetate": 1, "water": 1},
    "stages": [5, 10],
}
MA_ATOMS = {
    "methanol": {"C": 1, "H": 4, "O": 1},
    "acetic acid": {"C": 2, "H": 4, "O": 2},
    "methyl acetate": {"C": 3, "H": 6, "O": 2},
    "water": {"H": 2, "O": 1},
}

# ---- toluene disproportionation ---------------------------------------------
TL_COMPS = ["benzene", "toluene", "p-xylene"]
TL_RXN = {
    "stoich": {"toluene": -2, "benzene": 1, "p-xylene": 1},
    "key": "toluene", "k0": 3e-8, "Ea": 0.0, "orders": {"toluene": 2},
    "k0_rev": 7.5e-8, "Ea_rev": 0.0,
    "orders_rev": {"benzene": 1, "p-xylene": 1},
    "stages": [4, 11],
}
TL_C = {"benzene": 6, "toluene": 7, "p-xylene": 8}
TL_H = {"benzene": 6, "toluene": 8, "p-xylene": 10}


def _atom_in_out(streams, atoms, comps, element):
    return sum(s.molar_flow * sum(s.z.get(c, 0.0) * atoms[c].get(element, 0.0)
                                  for c in comps) for s in streams)


def _methyl_acetate_column(holdup: float) -> Flowsheet:
    params = dict(n_stages=15, feed_stage=10, reflux_ratio=5.0,
                  distillate_rate=5.556, P=90000.0, dP_stage=500.0,
                  max_iter=200)
    if holdup > 0:
        params["reactions"] = [MA_RXN]
        params["tray_holdup"] = holdup
    fs = Flowsheet(components=[Component(c) for c in MA_COMPS],
                   property_package="thermo:NRTL")
    fs.add(RigorousColumn("COL", params))
    fs.feed("F", "COL:in1", T=348.15, P=101325.0, molar_flow=12.5,
            z={"methanol": 0.4, "acetic acid": 0.4,
               "methyl acetate": 0.1, "water": 0.1})
    fs.connect("D", "COL:distillate", None)
    fs.connect("B", "COL:bottoms", None)
    fs.connect("qc", "COL:condenser_duty", None)
    fs.connect("qr", "COL:reboiler_duty", None)
    return fs


def _toluene_column(holdup: float) -> Flowsheet:
    params = dict(n_stages=15, feed_stage=8, reflux_ratio=3.0,
                  distillate_to_feed=0.5, max_iter=200)
    if holdup > 0:
        params["reactions"] = [TL_RXN]
        params["tray_holdup"] = holdup
    fs = Flowsheet(components=[Component(c) for c in TL_COMPS],
                   property_package="thermo:PR")
    fs.add(RigorousColumn("COL", params))
    fs.feed("F", "COL:in1", T=400.0, P=120000.0, molar_flow=10.0,
            z={"benzene": 0.05, "toluene": 0.90, "p-xylene": 0.05})
    fs.connect("D", "COL:distillate", None)
    fs.connect("B", "COL:bottoms", None)
    fs.connect("qc", "COL:condenser_duty", None)
    fs.connect("qr", "COL:reboiler_duty", None)
    return fs


def _meoh_conversion(fs) -> float:
    d, b = fs.streams["D"], fs.streams["B"]
    fed = 12.5 * 0.4
    out = d.molar_flow * d.z["methanol"] + b.molar_flow * b.z["methanol"]
    return (fed - out) / fed


# ============================ reversible kinetics ============================
def test_kinetic_reaction_reverse_term():
    """A reversible KineticReaction's net rate is forward minus reverse, and it
    vanishes at the equilibrium concentration ratio k_f/k_r."""
    rxn = KineticReaction.from_param(MA_RXN)
    # forward-only at zero product concentration
    conc0 = {"methanol": 1000.0, "acetic acid": 1000.0,
             "methyl acetate": 0.0, "water": 0.0}
    assert rxn.rate(conc0, 350.0) > 0
    # equilibrium: k_f C_M C_A = k_r C_E C_W -> net rate ~ 0.
    # K = k_f/k_r = 3e-8/6e-9 = 5 = (C_E C_W)/(C_M C_A).
    conc_eq = {"methanol": 100.0, "acetic acid": 100.0,
               "methyl acetate": 500.0, "water": 100.0}  # (500*100)/(100*100)=5
    assert rxn.rate(conc_eq, 350.0) == pytest.approx(0.0, abs=1e-12)
    # past equilibrium (excess products) -> net reverse (negative)
    conc_rev = {"methanol": 10.0, "acetic acid": 10.0,
                "methyl acetate": 900.0, "water": 900.0}
    assert rxn.rate(conc_rev, 350.0) < 0


def test_forward_only_reaction_unchanged():
    """Without reverse parameters the rate is forward-only (back-compat)."""
    fwd = KineticReaction.from_param(
        {"stoich": {"methanol": -1, "methyl acetate": 1}, "key": "methanol",
         "k0": 2.0, "Ea": 0.0})
    assert fwd.k0_rev == 0.0
    assert fwd.rate({"methanol": 3.0}, 300.0) == pytest.approx(6.0)


# ============================ methyl acetate (book) ==========================
def test_methyl_acetate_reaction_proceeds_and_enriches_distillate():
    fs = _methyl_acetate_column(holdup=0.3)
    rep = fs.solve()
    assert rep.converged
    d, b = fs.streams["D"], fs.streams["B"]
    # the reaction makes methyl acetate: more leaves than the feed carries
    meac_out = (d.molar_flow * d.z["methyl acetate"]
                + b.molar_flow * b.z["methyl acetate"])
    assert meac_out > 12.5 * 0.1 + 1e-3
    # the light ester concentrates in the distillate vs the bottoms
    assert d.z["methyl acetate"] > b.z["methyl acetate"]
    assert _meoh_conversion(fs) > 0.05


def test_methyl_acetate_element_balances_exact():
    fs = _methyl_acetate_column(holdup=0.3)
    fs.solve()
    feed = [fs.streams["F"]]
    prod = [fs.streams["D"], fs.streams["B"]]
    for el in ("C", "H", "O"):
        a_in = _atom_in_out(feed, MA_ATOMS, MA_COMPS, el)
        a_out = _atom_in_out(prod, MA_ATOMS, MA_COMPS, el)
        assert a_out == pytest.approx(a_in, rel=1e-9), el


def test_reactive_distillation_vanishes_as_holdup_goes_to_zero():
    """At negligible holdup the reactive column reproduces the non-reactive
    one — the reaction source is continuous in the holdup."""
    base = _methyl_acetate_column(holdup=0.0)
    base.solve()
    tiny = _methyl_acetate_column(holdup=1e-6)
    tiny.solve()
    assert _meoh_conversion(tiny) < 1e-3
    assert tiny.streams["D"].z["methyl acetate"] == pytest.approx(
        base.streams["D"].z["methyl acetate"], abs=2e-3)


def test_conversion_rises_with_holdup():
    lo = _methyl_acetate_column(holdup=0.1)
    hi = _methyl_acetate_column(holdup=0.3)
    lo.solve()
    hi.solve()
    assert _meoh_conversion(hi) > _meoh_conversion(lo) > 0.0


# ============================ toluene disproportionation =====================
def test_toluene_disproportionation_separates_products():
    fs = _toluene_column(holdup=0.5)
    rep = fs.solve()
    assert rep.converged
    d, b = fs.streams["D"], fs.streams["B"]
    # benzene (light) overhead, xylene (heavy) in the bottoms
    assert d.z["benzene"] > 0.4
    assert b.z["p-xylene"] > 0.4
    tol_out = d.molar_flow * d.z["toluene"] + b.molar_flow * b.z["toluene"]
    conv = (10.0 * 0.9 - tol_out) / (10.0 * 0.9)
    assert 0.3 < conv < 0.9


def test_toluene_carbon_hydrogen_balance_exact():
    fs = _toluene_column(holdup=0.5)
    fs.solve()
    d, b = fs.streams["D"], fs.streams["B"]
    f = fs.streams["F"]
    for atoms in (TL_C, TL_H):
        a_in = f.molar_flow * sum(f.z.get(c, 0.0) * atoms[c] for c in TL_COMPS)
        a_out = sum(s.molar_flow * sum(s.z.get(c, 0.0) * atoms[c]
                    for c in TL_COMPS) for s in (d, b))
        assert a_out == pytest.approx(a_in, rel=1e-9)


# ============================ error paths + round-trip =======================
def test_reactions_require_reboiled_bubble_point():
    fs = Flowsheet(components=[Component(c) for c in TL_COMPS],
                   property_package="thermo:PR")
    fs.add(RigorousColumn("COL", dict(
        n_stages=10, feed_stage=5, reflux_ratio=2.0, distillate_to_feed=0.4,
        reboiled=False, method="naphtali_sandholm",
        reactions=[TL_RXN], tray_holdup=0.5)))
    fs.feed("F", "COL:in1", T=400.0, P=120000.0, molar_flow=10.0,
            z={"benzene": 0.05, "toluene": 0.90, "p-xylene": 0.05})
    fs.connect("D", "COL:distillate", None)
    fs.connect("B", "COL:bottoms", None)
    fs.connect("qc", "COL:condenser_duty", None)
    fs.connect("qr", "COL:reboiler_duty", None)
    with pytest.raises(RigorousColumnError, match="reactive distillation"):
        fs.solve()


def test_reaction_stage_range_validated():
    bad = dict(MA_RXN)
    bad["stages"] = [5, 20]                        # 20 > n_stages
    fs = Flowsheet(components=[Component(c) for c in MA_COMPS],
                   property_package="thermo:NRTL")
    fs.add(RigorousColumn("COL", dict(
        n_stages=15, feed_stage=10, reflux_ratio=5.0, distillate_rate=5.556,
        reactions=[bad], tray_holdup=0.3)))
    fs.feed("F", "COL:in1", T=348.15, P=101325.0, molar_flow=12.5,
            z={"methanol": 0.4, "acetic acid": 0.4,
               "methyl acetate": 0.1, "water": 0.1})
    fs.connect("D", "COL:distillate", None)
    fs.connect("B", "COL:bottoms", None)
    fs.connect("qc", "COL:condenser_duty", None)
    fs.connect("qr", "COL:reboiler_duty", None)
    with pytest.raises(RigorousColumnError, match="out of range"):
        fs.solve()


def test_flow_round_trip_with_reactions():
    fs = _toluene_column(holdup=0.5)
    fs.solve()
    doc = to_dict(fs)
    assert doc == to_dict(from_dict(doc))
