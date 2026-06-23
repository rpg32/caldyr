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
