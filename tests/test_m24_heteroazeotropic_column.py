"""Integrated heteroazeotropic (3-phase) distillation column.

This is the column-solver half of the P5/P6 thermo wall. With the predictive
UNIFAC VLE+LLE package (test_thermo_unifac.py) the FLASH is honest about a
two-liquid tray, but the single-liquid bubble-point MESH still cannot host one:
a tray whose liquid lies inside the miscibility gap (water + cyclohexane
coexisting in the entrainer column) has no single-liquid bubble point, and the
formation-enthalpy spread (water -242 vs cyclohexane -123 kJ/mol) makes the bare
energy-balance pivot sign-indefinite.

The fix has two admissible parts, both validated here:
  * a VLLE-aware stage flash — the tray holdup is ONE liquid stream of the
    combined composition, but its equilibrium vapour is the one in equilibrium
    with BOTH settled layers (isoactivity), and its enthalpy is the molar-
    combined two-liquid enthalpy (caldyr.thermo.activity_pkg);
  * a sensible-basis energy balance in the reboiled MESH (subtract each stream's
    formation-enthalpy offset — admissible, a per-component constant), so the
    recurrence pivot stays a true latent heat.

The sequential bubble-point method still *oscillates* on the real ternary, so
the column is solved by the simultaneous-correction (Naphtali-Sandholm) method,
which converges it to a machine-precision energy balance. Reference system:
ethanol/water/cyclohexane, the entrainer of anhydrous-ethanol dehydration
(Hameed 2025 sec. 9.5.6).
"""
import pytest

from caldyr.core import Component, Flowsheet
from caldyr.unitops import RigorousColumn

P_ATM = 101325.0


def _ternary_flowsheet(method: str, n_stages: int = 16, max_iter: int = 80):
    fs = Flowsheet(
        components=[Component("ethanol"), Component("water"), Component("cyclohexane")],
        property_package="thermo:UNIFAC")
    fs.add(RigorousColumn("T100", {
        "n_stages": n_stages,
        "feeds": [{"stage": 3}, {"stage": 10}],     # cyclohexane near top, feed mid
        "reflux_ratio": 2.0,
        "distillate_rate": 18.0,
        "partial_condenser": True,
        "method": method,
        "max_iter": max_iter,
    }))
    fs.feed("SOLV", "T100:in1", T=298.15, P=P_ATM, molar_flow=15.0,
            z={"cyclohexane": 1.0})
    fs.feed("FEED", "T100:in2", T=343.0, P=P_ATM, molar_flow=27.8,
            z={"ethanol": 0.87, "water": 0.13})
    fs.connect("DIST", "T100:distillate", None)
    fs.connect("BOT", "T100:bottoms", None)
    fs.connect("QC", "T100:condenser_duty", None)
    fs.connect("QR", "T100:reboiler_duty", None)
    return fs


@pytest.mark.slow
def test_heteroazeotropic_column_converges_via_ns():
    """The two-liquid-tray entrainer column converges (Naphtali-Sandholm) to a
    machine-precision energy balance and an exact mass balance, and its overhead
    is a genuine heteroazeotrope (it decants into two liquids)."""
    fs = _ternary_flowsheet("naphtali_sandholm")
    fs.solve()
    col = fs.streams
    bot, dist = col["BOT"], col["DIST"]

    # Energy balance closed to machine precision (the headline diagnostic).
    col_unit = fs.units["T100"]
    assert abs(col_unit.design["energy_residual_rel"]) < 1e-6

    # Overall mass balance: feeds = products, component by component.
    for c in ("ethanol", "water", "cyclohexane"):
        fed = 15.0 * (1.0 if c == "cyclohexane" else 0.0) + 27.8 * (
            0.87 if c == "ethanol" else 0.13 if c == "water" else 0.0)
        out = bot.molar_flow * bot.z[c] + dist.molar_flow * dist.z[c]
        assert out == pytest.approx(fed, rel=1e-6, abs=1e-6)

    # The overhead is the cyclohexane/water/ethanol heteroazeotrope — it carries
    # the entrainer and the water overhead, and it splits into two liquids.
    assert dist.z["cyclohexane"] > 0.4
    assert dist.z["water"] > 0.05
    from caldyr.thermo import make_package
    pp = make_package("thermo:UNIFAC", ["ethanol", "water", "cyclohexane"])
    split = pp.flash_pt_3p(315.0, P_ATM, dict(dist.z))
    assert split.beta_light > 0.05 and split.beta_heavy > 0.05

    # Ethanol concentrates in the bottoms (the heavy product of the column).
    assert bot.z["ethanol"] > dist.z["ethanol"]


@pytest.mark.slow
def test_heteroazeotropic_bubble_point_oscillates():
    """Documents WHY the NS method is required: the sequential bubble-point MESH
    does not converge on the same two-liquid-tray column (it oscillates)."""
    fs = _ternary_flowsheet("bubble_point", max_iter=40)
    from caldyr.unitops.rigorous_column import RigorousColumnError
    with pytest.raises(RigorousColumnError):
        fs.solve()


# -- integrated decanting condenser (Hameed §9.5.6 anhydrous ethanol) ----------
def _decant_column(D: float, solv: float, n_stages: int = 20,
                   feeds=((2, ), (12, )), condenser_T: float = 313.15,
                   max_iter: int = 140) -> Flowsheet:
    """A reboiled entrainer column with an INTEGRATED decanting condenser: stage
    0's overhead is settled at ``condenser_T``, the cyclohexane-rich organic
    layer is refluxed in full (internally), and the ethanol/water aqueous layer
    is the distillate. Anhydrous-ethanol dehydration, Hameed 2025 §9.5.6."""
    fs = Flowsheet(
        components=[Component("ethanol"), Component("water"),
                    Component("cyclohexane")],
        property_package="thermo:UNIFAC")
    fs.add(RigorousColumn("T100", {
        "n_stages": n_stages,
        "feeds": [{"stage": feeds[0][0]}, {"stage": feeds[1][0]}],
        "reflux_ratio": 3.0,
        "distillate_rate": D,
        "method": "naphtali_sandholm",
        "reboiled": True,
        "decant_condenser": True,
        "condenser_T": condenser_T,
        "reflux_layer": "organic",
        "max_iter": max_iter,
    }))
    fs.feed("SOLV", "T100:in1", T=298.15, P=P_ATM, molar_flow=solv,
            z={"cyclohexane": 1.0})
    fs.feed("FEED", "T100:in2", T=343.0, P=P_ATM, molar_flow=27.8,
            z={"ethanol": 0.87, "water": 0.13})
    fs.connect("DIST", "T100:distillate", None)
    fs.connect("BOT", "T100:bottoms", None)
    fs.connect("QC", "T100:condenser_duty", None)
    fs.connect("QR", "T100:reboiler_duty", None)
    return fs


def test_decant_condenser_param_validation():
    """The decant mode validates its specs up front (no solve required)."""
    from caldyr.unitops.rigorous_column import RigorousColumnError

    def _solve(extra):
        params = {"n_stages": 12, "feed_stage": 6, "reflux_ratio": 2.0,
                  "distillate_rate": 5.0, "method": "naphtali_sandholm",
                  "reboiled": True, "decant_condenser": True}
        params.update(extra)
        fs = Flowsheet(
            components=[Component("ethanol"), Component("water"),
                        Component("cyclohexane")],
            property_package="thermo:UNIFAC")
        fs.add(RigorousColumn("C", params))
        fs.feed("F", "C:in1", T=343.0, P=P_ATM, molar_flow=27.8,
                z={"ethanol": 0.8, "water": 0.1, "cyclohexane": 0.1})
        fs.connect("D", "C:distillate", None)
        fs.connect("B", "C:bottoms", None)
        fs.connect("QC", "C:condenser_duty", None)
        fs.connect("QR", "C:reboiler_duty", None)
        fs.solve()

    # condenser_T is required in decant mode
    with pytest.raises(RigorousColumnError, match="condenser_T"):
        _solve({})
    # decant needs the simultaneous-correction reboiled column
    with pytest.raises(RigorousColumnError, match="naphtali_sandholm"):
        _solve({"condenser_T": 313.15, "method": "bubble_point"})
    # reflux_layer must name a layer
    with pytest.raises(RigorousColumnError, match="reflux_layer"):
        _solve({"condenser_T": 313.15, "reflux_layer": "middle"})
    # decant + partial condenser is contradictory
    with pytest.raises(RigorousColumnError, match="partial_condenser"):
        _solve({"condenser_T": 313.15, "partial_condenser": True})


@pytest.mark.slow
def test_decant_condenser_converges_and_recirculates():
    """The integrated decanting condenser converges to a machine-precision
    energy balance and an exact mass balance, its overhead is a genuine
    heteroazeotrope (decants into two liquids), and the cyclohexane entrainer
    RECIRCULATES internally (a cyclohexane-rich organic reflux that never leaves
    the column)."""
    fs = _decant_column(D=4.0, solv=8.0)
    fs.solve()
    col = fs.units["T100"]
    d = col.design
    bot, dist = fs.streams["BOT"], fs.streams["DIST"]

    assert d["decant_condenser"] is True
    # machine-precision energy balance + exact overall mass balance.
    assert abs(d["energy_residual_rel"]) < 1e-6
    for c in ("ethanol", "water", "cyclohexane"):
        fed = 8.0 * (1.0 if c == "cyclohexane" else 0.0) + 27.8 * (
            0.87 if c == "ethanol" else 0.13 if c == "water" else 0.0)
        out = bot.molar_flow * bot.z[c] + dist.molar_flow * dist.z[c]
        assert out == pytest.approx(fed, rel=1e-6, abs=1e-6)

    # the overhead vapour is the ternary heteroazeotrope: it decants in two.
    from caldyr.thermo import make_package
    pp = make_package("thermo:UNIFAC", ["ethanol", "water", "cyclohexane"])
    y0 = d["y_profile"][0]
    split = pp.flash_pt_3p(313.15, P_ATM, dict(y0))
    assert split.beta_light > 0.05 and split.beta_heavy > 0.05

    # the entrainer recirculates internally: a real organic reflux that is
    # cyclohexane-rich (far richer than the net aqueous distillate).
    assert d["organic_reflux"] > 0.0
    assert d["R"] > 0.0
    assert d["x_organic"]["cyclohexane"] > 0.6
    assert d["x_organic"]["cyclohexane"] > dist.z["cyclohexane"]


@pytest.mark.slow
def test_decant_condenser_reaches_anhydrous_regime():
    """The integrated decant makes the LARGE-distillate regime tractable — the
    regime an external decanter cannot converge (P6). Warm-started distillate-
    rate continuation drives the bottoms toward anhydrous ethanol: water and
    cyclohexane both fall monotonically and the column stays converged with an
    exact mass balance at every step."""
    col = RigorousColumn("T100", {
        "n_stages": 30, "feeds": [{"stage": 3}, {"stage": 14}],
        "reflux_ratio": 3.0, "distillate_rate": 4.2,
        "method": "naphtali_sandholm", "reboiled": True,
        "decant_condenser": True, "condenser_T": 305.0,
        "reflux_layer": "organic", "max_iter": 120,
    })
    last_water = 1.0
    bot = None
    for D, solv in [(4.2, 8.0), (8.0, 6.0), (12.0, 4.0), (16.0, 3.2),
                    (18.0, 3.0)]:
        col.params["distillate_rate"] = D
        fs = Flowsheet(
            components=[Component("ethanol"), Component("water"),
                        Component("cyclohexane")],
            property_package="thermo:UNIFAC")
        fs.add(col)
        fs.feed("SOLV", "T100:in1", T=298.15, P=P_ATM, molar_flow=solv,
                z={"cyclohexane": 1.0})
        fs.feed("FEED", "T100:in2", T=343.0, P=P_ATM, molar_flow=27.8,
                z={"ethanol": 0.87, "water": 0.13})
        fs.connect("DIST", "T100:distillate", None)
        fs.connect("BOT", "T100:bottoms", None)
        fs.connect("QC", "T100:condenser_duty", None)
        fs.connect("QR", "T100:reboiler_duty", None)
        fs.solve()
        d = col.design
        bot = fs.streams["BOT"]
        assert abs(d["energy_residual_rel"]) < 1e-6
        # water in the bottoms falls monotonically as D rises (entrainer drags
        # the water overhead in the ethanol-rich aqueous layer).
        assert bot.z["water"] < last_water + 1e-9
        last_water = bot.z["water"]

    # deep in the large-D regime the bottoms is near-anhydrous, near-pure ethanol
    assert bot is not None
    assert bot.z["ethanol"] > 0.98
    assert bot.z["water"] < 0.015
    assert bot.z["cyclohexane"] < 0.015
