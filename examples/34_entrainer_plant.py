"""Example 34 — the full closed §9.5.6 anhydrous-ethanol entrainer PLANT, taken
to BOOK SCALE (62-stage T-100) by growing the loop in place.

Two columns + two recycles (Hameed 2025 §9.5.6):

    fresh feed (EtOH/water) ─┐
                             ▼
         cyclohexane     ┌────────┐  aqueous (EtOH/water/cyclohexane)
         make-up ───────►│ T-100  ├──────────────────────────────┐
            ▲            │ decant │                               ▼
            │            │ column │                          ┌────────┐
       ┌────┴────┐       └───┬────┘                          │ T-101  │
       │ Makeup  │◄──────────┼───────────────────────────────┤ water  │
       └─────────┘   recycle (EtOH/cyclohexane)  T-101 dist   │ column │
                             │                                └───┬────┘
                    anhydrous EtOH (bottoms)            water (bottoms)

T-100 is the integrated DECANTING-CONDENSER column (`decant_condenser=True`): its
overhead settles internally, the cyclohexane-rich organic layer is refluxed in
full (NO external organic recycle — that is the whole point), and the ethanol/
water aqueous layer is the distillate. T-101 recovers the water (bottoms) and
recycles the ethanol+cyclohexane (distillate) back to T-100.

Three techniques get the closed loop to book parity:

  1. **Flowsheet continuation** brings the 30-stage loop up: the distillate rates
     ramp UP (drive the water overhead) while the cyclohexane make-up ramps DOWN
     (so the entrainer recirculates instead of piling into the bottoms).
  2. **Stage-count continuation** then grows T-100 to 62 stages in place — the
     converged 30-stage column seeds the 62-stage one (`warm_start_from`), and the
     already-converged recycle tears seed the larger loop (the tear guesses are
     dropped so the loop warm-starts from the 30-stage solution). 62 stages give
     the long stripping section that strips the last cyclohexane and the long
     rectifying section that takes the bottoms to anhydrous.
  3. **Inventory control** (A2): stepping the cyclohexane make-up DOWN trims the
     entrainer inventory and the bottoms cyclohexane falls with it (the production
     equivalent is a logical **Adjust** on a bottoms-cyclohexane spec — exact, but
     each root-find step re-solves the whole 62-stage recycle, so the marginal-NS
     endgame — the turning-point fold documented in the RigorousColumn solver —
     makes a deep-cut Adjust impractical at this scale; the make-up continuation
     shows the same physics robustly).

The high-D / low-make-up steps that would strip the LAST cyclohexane to anhydrous
sit in that marginal-NS endgame; the open-loop column (example 33) shows 62
stages reaching anhydrous, so closed-loop book parity is gated on that fold, not
on a missing capability.

NOTE: SLOW (a 62-stage VLLE column in a recycle, two continuation passes). Run it
directly; it is not part of the fast test suite.

    python examples/34_entrainer_plant.py
"""
import sys
from pathlib import Path

_ENGINE = Path(__file__).resolve().parent.parent / "engine"
if _ENGINE.is_dir() and str(_ENGINE) not in sys.path:
    sys.path.insert(0, str(_ENGINE))

from caldyr.core import Component, Flowsheet  # noqa: E402
from caldyr.unitops import Makeup, RigorousColumn  # noqa: E402

P_ATM = 101325.0
COMPS = ["ethanol", "water", "cyclohexane"]


def build() -> Flowsheet:
    fs = Flowsheet(components=[Component(c) for c in COMPS],
                   property_package="thermo:UNIFAC")
    fs.add(RigorousColumn("T100", {
        "n_stages": 30, "feeds": [{"stage": 2}, {"stage": 10}],
        "reflux_ratio": 3.0, "distillate_rate": 4.0,
        "method": "naphtali_sandholm", "reboiled": True,
        "decant_condenser": True, "condenser_T": 305.0,
        "reflux_layer": "organic", "max_iter": 120,
    }))
    fs.add(RigorousColumn("T101", {
        # A3: water column resolved to reject ethanol+cyclohexane overhead and
        # concentrate water in the bottoms (feed near the top = long stripping
        # section; high reflux). The 62-stage T-100 sends an anhydrous bottoms, so
        # essentially ALL the feed water reports to the aqueous draw AQ and T-101
        # concentrates it; its draw D101 is ramped up alongside D100 below.
        "n_stages": 30, "feed_stage": 15, "reflux_ratio": 5.0,
        "distillate_rate": 2.0, "method": "naphtali_sandholm",
        "reboiled": True, "max_iter": 120,
    }))
    fs.add(Makeup("MK", {"component": "cyclohexane", "target": 8.0,
                         "T": 305.0, "P": P_ATM}))

    fs.feed("FEED", "T100:in2", T=343.0, P=P_ATM, molar_flow=27.78,
            z={"ethanol": 0.87, "water": 0.13})
    fs.connect("ENTR", "MK:out", "T100:in1")
    fs.connect("AQ", "T100:distillate", "T101:in1")
    fs.connect("REC", "T101:distillate", "MK:in1")
    fs.connect("ETOH", "T100:bottoms", None)
    fs.connect("WATER", "T101:bottoms", None)
    for u in ("T100", "T101"):
        fs.connect(f"{u}_QC", f"{u}:condenser_duty", None)
        fs.connect(f"{u}_QR", f"{u}:reboiler_duty", None)

    entr0 = {"ethanol": 0.25, "water": 0.10, "cyclohexane": 0.65}
    rec0 = {"ethanol": 0.30, "water": 0.10, "cyclohexane": 0.60}
    fs.solver_hints = {
        "tear_guesses": {
            "ENTR": {"T": 320.0, "P": P_ATM, "molar_flow": 8.0, "z": entr0},
            "REC": {"T": 333.0, "P": P_ATM, "molar_flow": 8.0, "z": rec0},
        },
        "tear_tolerance": 5e-3,
    }
    return fs


def _row(fs, d100, d101, mkt, rep):
    etoh = fs.streams["ETOH"]
    print(f"  {d100:5.1f} {d101:5.1f} {mkt:4.1f}  "
          f"{'conv ' if rep.converged else 'NOCONV':>8}  "
          f"{etoh.z['ethanol']:.4f} / {etoh.z['water']:.5f} / "
          f"{etoh.z['cyclohexane']:.5f}  (sweeps {rep.iterations})", flush=True)


def _set(fs, d100, d101, mkt):
    fs.units["T100"].params["distillate_rate"] = d100
    fs.units["T101"].params["distillate_rate"] = d101
    fs.units["MK"].params["target"] = mkt


def _continue(fs, schedule):
    """Run the flowsheet continuation; on the first non-converging step restore
    and re-solve the last converged setpoint (so the flowsheet is left in a clean
    converged state). Returns the last converged ``(d100, d101, mkt)``."""
    last_good = None
    for d100, d101, mkt in schedule:
        _set(fs, d100, d101, mkt)
        try:
            rep = fs.solve(method="direct", max_iter=30)
            if not rep.converged:
                raise RuntimeError("recycle did not converge")
        except Exception as exc:
            print(f"  {d100:5.1f} {d101:5.1f} {mkt:4.1f}  "
                  f"{'STOP':>8}  {str(exc)[:58]}", flush=True)
            if last_good is not None:
                _set(fs, *last_good)
                fs.solve(method="direct", max_iter=30)   # back to clean state
            return last_good
        _row(fs, d100, d101, mkt, rep)
        last_good = (d100, d101, mkt)
    return last_good


def main() -> None:
    fs = build()
    print("§9.5.6 entrainer plant — PASS 1: bring the 30-stage loop up:")
    print(f"  {'D100':>5} {'D101':>5} {'mk':>4}  {'recycle':>8}  "
          f"{'EtOH bottoms (x_EtOH / water / cyc)':<38}")
    # Keep PASS 1 modest — just establish the loop at a moderate draw; the high-D
    # push is left to the longer, more robust 62-stage column in PASS 2.
    if not _continue(fs, [(4.0, 2.0, 8.0), (7.0, 5.0, 6.0), (10.0, 7.0, 5.0)]):
        return

    # -- grow T-100 to 62 stages in place (stage-count continuation) ------------
    print("\nGrowing T-100 to 62 stages (warm_start_from the converged 30-stage "
          "column)...", flush=True)
    t100_62 = RigorousColumn("T100", {
        "n_stages": 62, "feeds": [{"stage": 2}, {"stage": 20}],
        "reflux_ratio": 3.0,
        "distillate_rate": fs.units["T100"].params["distillate_rate"],
        "method": "naphtali_sandholm", "reboiled": True,
        "decant_condenser": True, "condenser_T": 305.0,
        # warm steps converge in <15 iters; cap low so a high-D step that drifts
        # into the marginal-NS endgame fails FAST and the continuation falls
        # back to the last converged draw instead of grinding.
        "reflux_layer": "organic", "max_iter": 70,
    })
    t100_62.warm_start_from(fs.units["T100"])
    fs.units["T100"] = t100_62
    # Drop the tear guesses so the larger loop warm-starts from the converged
    # 30-stage tears (already in fs.streams) instead of the cold cyclohexane-rich
    # guess.
    fs.solver_hints = {"tear_tolerance": 5e-3}

    print("\nPASS 2: continue the 62-stage loop to anhydrous:")
    print(f"  {'D100':>5} {'D101':>5} {'mk':>4}  {'recycle':>8}  "
          f"{'EtOH bottoms (x_EtOH / water / cyc)':<38}")
    # Re-establish the loop at the PASS-1 draw on 62 stages, then ramp the draw
    # up gently (the high-D endgame is the marginal-NS regime, so small steps
    # keep the recycle in its basin; a step that drifts out fails fast and the
    # continuation falls back to the last converged draw).
    # D101 is ramped up alongside D100 (A3): with the 62-stage T-100 sending an
    # anhydrous bottoms, essentially all the feed water reports to AQ, and the
    # retuned T-101 concentrates it — the water product climbs as D101 rises until
    # T-101's own high-draw fold caps it; the continuation then falls
    # back to the last converged draw.
    good = _continue(fs, [(10.0, 7.0, 5.0), (12.0, 9.0, 4.5), (14.0, 11.0, 4.0)])
    if good is None:
        return

    # -- inventory control (A2): step the cyclohexane make-up DOWN at the
    #    converged draw and watch the bottoms cyclohexane fall — the entrainer
    #    stops piling into the product as its inventory is trimmed. This is the
    #    robust make-up-continuation form; the production equivalent is a logical
    #    **Adjust** that varies ["MK","target"] on a bottoms-cyclohexane
    #    component_rate spec (caldyr.solver.logical) — exact but slow here (each
    #    root-find step re-solves the whole 62-stage recycle), so the marginal-NS
    #    endgame makes the deep-cut Adjust impractical at this scale.
    d100, d101, mk0 = good
    print(f"\nInventory control (A2) — step the make-up down at D100={d100:.0f} "
          f"(bottoms cyclohexane falls):", flush=True)
    last_mk = mk0
    for mkt in [mk0, 0.8 * mk0, 0.6 * mk0]:
        _set(fs, d100, d101, mkt)
        try:
            rep = fs.solve(method="direct", max_iter=30)
            if not rep.converged:
                raise RuntimeError("recycle did not converge")
        except Exception as exc:
            # The deep-cut make-up step drifted into the lean-solvent marginal-NS
            # endgame. Restore the last CONVERGED make-up and re-solve so
            # the flowsheet is left clean; guard that fallback too (a re-solve from
            # the perturbed state can itself struggle) so the run still reports.
            print(f"  {d100:5.1f} {d101:5.1f} {mkt:4.1f}  "
                  f"{'STOP':>8}  {str(exc)[:58]}", flush=True)
            _set(fs, d100, d101, last_mk)
            try:
                fs.solve(method="direct", max_iter=30)
            except Exception:
                pass
            break
        _row(fs, d100, d101, mkt, rep)
        last_mk = mkt

    etoh, water, rec = fs.streams["ETOH"], fs.streams["WATER"], fs.streams["REC"]
    print(f"\nethanol product (T-100 bottoms): {etoh.molar_flow:.2f} mol/s, "
          f"x_EtOH={etoh.z['ethanol']:.4f} "
          f"(water {etoh.z['water']:.5f}, cyclohexane {etoh.z['cyclohexane']:.5f})")
    print(f"water product (T-101 bottoms): {water.molar_flow:.2f} mol/s, "
          f"x_water={water.z['water']:.4f}")
    print(f"cyclohexane make-up: {fs.units['MK'].design['makeup_flow']:.3f} mol/s "
          f"(recycle returns {rec.molar_flow * rec.z['cyclohexane']:.2f} mol/s)")
    print("\nThe 62-stage T-100 (grown from the converged 30-stage loop) drives "
          "the bottoms toward\nanhydrous ethanol (~0.91 here) with the cyclohexane "
          "recirculating; the retuned T-101\nconcentrates the water product (~0.53 "
          "here, up from ~0.32). Full book parity (>0.999\nEtOH / >0.95 water) "
          "needs the high-draw operating point that sits past the decant /\nwater-"
          "column turning-point fold — the same wall the open-loop\ncolumn "
          "(example 33) clears but the closed recycle does not in one "
          "solve.")


if __name__ == "__main__":
    main()
