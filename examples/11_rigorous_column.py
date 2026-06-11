"""Rigorous (MESH) vs shortcut (FUG) distillation, side by side.

The same benzene/toluene separation solved two ways:

1. ShortcutColumn designs the column with Fenske-Underwood-Gilliland for
   95%/95% key recoveries at R = 1.3 R_min.
2. RigorousColumn re-rates that exact design (same stages, feed stage, reflux
   ratio and distillate rate) with full tray-by-tray MESH balances — the
   bubble-point (Wang-Henke) method, real PR K-values and enthalpies on every
   stage — and reports the converged stage profiles.

The two should (and do) agree on the split to within a few percent, which is
both a validation of the rigorous column and a demonstration of what the
shortcut method glosses over (stage-by-stage traffic and temperature). The
rigorous column is then sized and costed exactly like the shortcut one
(tower + trays + condenser + reboiler, Turton bare-module).

    python examples/11_rigorous_column.py
"""
import sys
from pathlib import Path

_ENGINE = Path(__file__).resolve().parent.parent / "engine"
if _ENGINE.is_dir() and str(_ENGINE) not in sys.path:
    sys.path.insert(0, str(_ENGINE))

from caldyr.core import Component, Flowsheet  # noqa: E402
from caldyr.economics import TEAConfig, analyze  # noqa: E402
from caldyr.unitops import RigorousColumn, ShortcutColumn  # noqa: E402

P_ATM = 101325.0


def build(column) -> Flowsheet:
    fs = Flowsheet(components=[Component("benzene"), Component("toluene")],
                   property_package="thermo:PR")
    fs.add(column)
    fs.feed("FEED", "COL:in1", T=365.0, P=P_ATM, molar_flow=100.0,
            z={"benzene": 0.5, "toluene": 0.5})
    for sid, port in (("DIST", "distillate"), ("BOT", "bottoms"),
                      ("QC", "condenser_duty"), ("QR", "reboiler_duty")):
        fs.connect(sid, f"COL:{port}", None)
    return fs


def recovery(fs, comp, product) -> float:
    feed, prod = fs.streams["FEED"], fs.streams[product]
    return prod.molar_flow * prod.z[comp] / (feed.molar_flow * feed.z[comp])


def money(x: float) -> str:
    return f"${x:,.0f}"


def main() -> None:
    # -- 1. shortcut design (FUG) -------------------------------------------
    fs_fug = build(ShortcutColumn("COL", {
        "light_key": "benzene", "heavy_key": "toluene",
        "recovery_light": 0.95, "recovery_heavy": 0.95,
        "rr_factor": 1.3, "P": P_ATM,
    }))
    rep_fug = fs_fug.solve()
    assert rep_fug.converged
    d_fug = fs_fug.units["COL"].design

    # -- 2. rigorous re-rate of the same design (MESH) ------------------------
    # FUG's N counts the reboiler but not the total condenser; RigorousColumn
    # numbers stages 1..n_stages including both, so n_stages = N + 1 and the
    # feed moves down one slot.
    n_stages = round(d_fug["N"]) + 1
    fs_rig = build(RigorousColumn("COL", {
        "n_stages": n_stages, "feed_stage": d_fug["feed_stage"] + 1,
        "reflux_ratio": d_fug["R"], "distillate_rate": d_fug["D"], "P": P_ATM,
    }))
    rep_rig = fs_rig.solve()
    assert rep_rig.converged
    d_rig = fs_rig.units["COL"].design

    print("=" * 70)
    print("  BENZENE/TOLUENE: SHORTCUT (FUG) vs RIGOROUS (MESH, bubble-point)")
    print("=" * 70)

    rows = [
        ("theoretical stages", f"{d_fug['N']:.1f} (+ total cond.)",
         f"{n_stages} (incl. cond. + reb.)"),
        ("feed stage (from top)", f"{d_fug['feed_stage']}",
         f"{d_fug['feed_stage'] + 1}"),
        ("reflux ratio R", f"{d_fug['R']:.3f}", f"{d_rig['R']:.3f}  (specified)"),
        ("distillate D, mol/s", f"{d_fug['D']:.1f}", f"{d_rig['D']:.1f}  (specified)"),
        ("benzene -> distillate", "95.0%  (specified)",
         f"{recovery(fs_rig, 'benzene', 'DIST') * 100:.1f}%  (computed)"),
        ("toluene -> bottoms", "95.0%  (specified)",
         f"{recovery(fs_rig, 'toluene', 'BOT') * 100:.1f}%  (computed)"),
        ("condenser duty, MW", f"{d_fug['Q_condenser'] / 1e6:.3f}",
         f"{d_rig['Q_condenser'] / 1e6:.3f}"),
        ("reboiler duty, MW", f"{d_fug['Q_reboiler'] / 1e6:.3f}",
         f"{d_rig['Q_reboiler'] / 1e6:.3f}"),
        ("T top, K", f"{d_fug['T_top']:.2f}", f"{d_rig['T_top']:.2f}"),
        ("T bottom, K", f"{d_fug['T_bottom']:.2f}", f"{d_rig['T_bottom']:.2f}"),
    ]
    print(f"\n  {'':<24} {'shortcut (FUG)':>22} {'rigorous (MESH)':>26}")
    for name, a, b in rows:
        print(f"  {name:<24} {a:>22} {b:>26}")

    print(f"\nMESH convergence: {d_rig['iterations']} bubble-point iterations, "
          f"{d_rig['flash_calls']} stage flashes, "
          f"max |dT| = {d_rig['max_dT']:.2e} K, "
          f"energy residual = {d_rig['energy_residual_rel']:.1e}")

    print("\nConverged stage profiles (stage 1 = condenser, "
          f"{n_stages} = reboiler):")
    print(f"  {'stage':>5} {'T (K)':>8} {'L (mol/s)':>10} {'V (mol/s)':>10} "
          f"{'x_benzene':>10} {'y_benzene':>10}")
    for j in range(n_stages):
        marker = " <- feed" if j + 1 == d_fug["feed_stage"] + 1 else ""
        print(f"  {j + 1:>5} {d_rig['T_profile'][j]:>8.2f} "
              f"{d_rig['L_profile'][j]:>10.2f} {d_rig['V_profile'][j]:>10.2f} "
              f"{d_rig['x_profile'][j]['benzene']:>10.4f} "
              f"{d_rig['y_profile'][j]['benzene']:>10.4f}{marker}")

    # -- 3. economics (identical sizing/costing path to the shortcut column) --
    res = analyze(fs_rig, rep_rig, TEAConfig(product_component="benzene",
                                             product_min_fraction=0.9))
    print("\nEquipment (Turton bare-module, installed):")
    print(f"  {'item':<16} {'type':<18} {'size':>16} {'Cbm':>12}")
    for size, cost in zip(res.sizes, res.costs):
        attr = f"{size.attribute:.1f} {size.attribute_name.split('_')[0]}"
        if size.quantity > 1:
            attr = f"{size.quantity} x {attr}"
        print(f"  {size.unit_id:<16} {size.equipment_type:<18} {attr:>16} "
              f"{money(cost.bare_module):>12}")
    print(f"\nCapital (TCI) ............. {money(res.capital.tci)}")
    print(f"Opex (COM, $/yr) .......... {money(res.opex.total)} "
          f"(utilities {money(res.opex.utilities)})")
    print(f"LCOP ...................... ${res.profitability.lcop:.3f}/kg benzene")


if __name__ == "__main__":
    main()
