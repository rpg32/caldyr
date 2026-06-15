"""M15 demo: a steam-stripped crude atmospheric distillation tower.

The tower follows the configuration of Hameed, "Chemical Process Simulations
using Aspen HYSYS" (Wiley 2025), sec. 10.2.1 "Refinery of Crude Oil-1" (book
pages 330-348, figs. 10.16 / 10.33-10.45):

    crude -> PreFlash drum -> furnace (the liquid only) -> mix back with the
    flash vapor -> the combined feed into a partial-condenser column that has
    NO reboiler. The bottom vapor is open stripping steam; kerosene and diesel
    leave as liquid side draws; pumparounds remove heat at intermediate stages
    (modeled as `stage_duties`, Kister's intermediate-condenser reduction).

The non-reboiled RigorousColumn runs the bubble-point MESH method: K-value
bubble temperatures per stage, vapor traffic from formation-offset-conditioned
envelope energy balances, with the open-steam bottom stage damped out of its
limit cycle (see `caldyr/unitops/rigorous_column.py`). The single column
specification is the distillate rate; the reflux ratio is an OUTPUT.

SCOPE / DEVIATIONS FROM THE BOOK (stated honestly):

* **Light crude.** The book's crude runs to a 1410 F endpoint (resid up to
  NBP ~736 C). With the cubic-EOS `thermo` backend, a bubble-point flash of a
  near-pure heavy petroleum pseudo-component returns an UNPHYSICAL saturated
  state (vapor enthalpy below the liquid's) once the NBP passes ~light gas oil,
  which no energy-balance MESH method can resolve. This example therefore uses
  a *light* synthetic crude (TBP 95-435 F) whose every cut stays in the
  property package's physically valid range. A full resid-bearing crude needs
  a property backend with corrected pseudo-component enthalpies (out of scope
  for the column solver, which sits above the PropertyPackage boundary).
* **Side draws, not side strippers; pumparounds as stage duties** (the book's
  kerosene/diesel are reboiled/steam-stripped side columns).
* **No book product table.** Sec. 10.2.1 ends at "Click Run and see the
  results" with no numeric slate, so the printed cut yields are compared to the
  assay's own TBP-implied cuts, not a book table.
* **method='sum_rates' does NOT converge this tower** (it limit-cycles); the
  converging method is the default 'bubble_point'.

    python examples/20_crude_tower.py
"""
import sys
import time
from pathlib import Path

_ENGINE = Path(__file__).resolve().parent.parent / "engine"
if _ENGINE.is_dir() and str(_ENGINE) not in sys.path:
    sys.path.insert(0, str(_ENGINE))

from caldyr.assay import (  # noqa: E402
    LIGHT_END_SG,
    api_to_sg,
    characterize_assay,
)
from caldyr.core import Flowsheet  # noqa: E402
from caldyr.core.components_db import molar_mass, resolve_component  # noqa: E402
from caldyr.unitops import (  # noqa: E402
    FlashDrum,
    Heater,
    Mixer,
    RigorousColumn,
)

PSIA = 6894.757293
BTU_PER_H = 0.29307107          # W per Btu/h


def degF(f: float) -> float:
    return (f - 32.0) / 1.8 + 273.15


def toF(k: float) -> float:
    return k * 1.8 - 459.67


# -- a light synthetic crude (naphtha .. light gas oil) ------------------------
TBP = [(0, 95.0), (10, 150.0), (20, 195.0), (30, 235.0), (40, 275.0),
       (50, 310.0), (60, 345.0), (70, 375.0), (80, 400.0), (90, 420.0),
       (98, 435.0)]
MW = [(0, 76.0), (10, 92.0), (20, 106.0), (30, 120.0), (40, 133.0),
      (50, 145.0), (60, 157.0), (70, 168.0), (80, 178.0), (90, 186.0),
      (98, 193.0)]
API = [(13, 72.0), (33, 64.0), (57, 56.0), (74, 50.0), (91, 44.0)]
LIGHT_ENDS = {"isobutane": 0.8, "n-butane": 1.5,
              "isopentane": 2.0, "n-pentane": 2.5}
BULK_API = 58.0
N_CUTS = 8
N_STAGES = 20
NAPHTHA_HI, KERO_HI, DIESEL_HI = 0.30, 0.52, 0.70


def main() -> None:
    t_all = time.perf_counter()
    res = characterize_assay(
        [(v, degF(f)) for v, f in TBP], kind="TBP", api_gravity=BULK_API,
        sg_curve=[(v, api_to_sg(a)) for v, a in API], mw_curve=MW,
        n_cuts=N_CUTS, light_ends=LIGHT_ENDS)

    ladder = [(cid, fr["vol_frac"], fr["mole_frac"])
              for cid, fr in res.light_end_fractions.items()]
    ladder += [(c.id, c.vol_frac, c.mole_frac) for c in res.cuts]

    def slice_moles(lo: float, hi: float) -> dict[str, float]:
        out: dict[str, float] = {}
        cum = 0.0
        for cid, dv, dn in ladder:
            a, b = cum, cum + dv
            cum = b
            ov = max(0.0, min(b, hi) - max(a, lo))
            if ov > 0.0:
                out[cid] = out.get(cid, 0.0) + dn * ov / dv
        return out

    vdot = 100000 * 0.158987294928 / 86400.0
    rho_w = 999.016
    sg_of = {c.id: c.SG for c in res.cuts}
    sg_of.update({cid: LIGHT_END_SG[cid] for cid in LIGHT_ENDS})
    mw_of = {c.id: c.MW for c in res.cuts}
    mw_of.update({cid: molar_mass(cid) for cid in LIGHT_ENDS})
    z = res.mole_fractions()
    vol_per_mol = sum(zi * mw_of[c] / (sg_of[c] * rho_w) for c, zi in z.items())
    F = vdot / vol_per_mol
    steam = 7500 * 0.45359237 / 3600.0 / 0.01801528

    naph = sum(slice_moles(0.0, NAPHTHA_HI).values()) * F
    kero = sum(slice_moles(NAPHTHA_HI, KERO_HI).values()) * F
    diesel = sum(slice_moles(KERO_HI, DIESEL_HI).values()) * F

    print("=" * 78)
    print("  CRUDE ATMOSPHERIC TOWER  (after Hameed 2025 sec. 10.2.1)")
    print("=" * 78)
    print(f"\n  light synthetic crude: {len(res.cuts)} pseudo-cuts + "
          f"{len(LIGHT_ENDS)} light ends + water (steam)")
    print(f"  crude feed   {F:8.1f} mol/s   "
          f"(100,000 bbl/day, {BULK_API:.0f} API)")
    print(f"  strip steam  {steam:8.1f} mol/s   (7500 lb/h on the bottom stage)")
    print("  TBP-implied cut yields (LV%):"
          f" naphtha {naph:.0f}, kerosene {kero:.0f}, diesel {diesel:.0f} mol/s")

    # -- flowsheet: preflash -> furnace -> mixer -> tower ----------------------
    comps = res.components() + [resolve_component("water")]
    fs = Flowsheet(components=comps, property_package="thermo:PR")
    fs.add(FlashDrum("PREFLASH", {"P": 75.0 * PSIA}))
    fs.add(Heater("FURNACE", {"T_out": degF(410.0), "dP": 10.0 * PSIA}))
    fs.add(Mixer("MIX"))
    n = N_STAGES
    fs.add(RigorousColumn("TOWER", {
        "n_stages": n, "feeds": [{"stage": n - 2}, {"stage": n}],
        "distillate_rate": naph + steam,
        "side_draws": [{"stage": 5, "phase": "liquid", "rate": kero},
                       {"stage": 10, "phase": "liquid", "rate": diesel}],
        "stage_duties": [{"stage": 2, "duty": -25e6 * BTU_PER_H},
                         {"stage": 10, "duty": -15e6 * BTU_PER_H}],
        "P": 19.7 * PSIA, "dP_stage": 13.0 * PSIA / (n - 1),
        "reboiled": False, "method": "bubble_point",
        "partial_condenser": True, "max_iter": 500,
    }))
    fs.connect("crude", None, "PREFLASH:in1")
    fs.connect("preflash_vap", "PREFLASH:vapor", "MIX:in1")
    fs.connect("preflash_liq", "PREFLASH:liquid", "FURNACE:in1")
    fs.connect("hot_crude", "FURNACE:out", "MIX:in2")
    fs.connect("tower_feed", "MIX:out", "TOWER:in1")
    fs.connect("naphtha", "TOWER:distillate", None)
    fs.connect("residue", "TOWER:bottoms", None)
    fs.connect("kerosene", "TOWER:side1", None)
    fs.connect("diesel", "TOWER:side2", None)
    fs.connect("qc", "TOWER:condenser_duty", None)
    fs.connect("qr", "TOWER:reboiler_duty", None)
    fs.feed("crude", "PREFLASH:in1", T=degF(450.0), P=75.0 * PSIA,
            molar_flow=F, z=z)
    fs.feed("steam", "TOWER:in2", T=degF(375.0), P=150.0 * PSIA,
            molar_flow=steam, z={"water": 1.0})

    t0 = time.perf_counter()
    report = fs.solve()
    t_solve = time.perf_counter() - t0
    assert report.converged

    des = fs.units["TOWER"].design
    print("\n" + "=" * 78)
    print("  CONVERGED  (method='bubble_point'; sum_rates limit-cycles here)")
    print("=" * 78)
    print(f"\n  iterations {des['iterations']}   "
          f"reflux ratio R = {des['R']:.2f} (an OUTPUT — no reboiler)   "
          f"in {t_solve:.1f} s")
    print(f"  condenser duty   {des['Q_condenser'] / 1e6:8.1f} MW")
    print(f"  energy residual  {des['energy_residual_rel']:.1e}  "
          f"(condenser-stage vs overall balance)")
    print("\n  stage temperature profile (F), top -> bottom:")
    print("   ", " ".join(f"{toF(t):.0f}" for t in des["T_profile"]))
    print("    (rises down the main column; the last two stages are cooler — "
          "the cold")
    print("     stripping steam quenches the open-steam bottom)")

    # -- product slate ---------------------------------------------------------
    print("\n  PRODUCT SLATE (ordered by boiling range):")
    print(f"  {'product':<10}{'mol/s':>10}{'T (F)':>9}{'mean MW':>9}"
          f"{'water':>8}")
    crude = fs.streams["crude"]
    for sid in ("naphtha", "kerosene", "diesel", "residue"):
        s = fs.streams[sid]
        mean_mw = sum(s.z.get(c, 0.0) * molar_mass(c)
                      for c in fs.component_ids) * 1000.0
        print(f"  {sid:<10}{s.molar_flow:>10.1f}{toF(s.T):>9.1f}"
              f"{mean_mw:>9.0f}{s.z.get('water', 0.0):>8.3f}")

    # the book reports NO numeric slate for this example (sec. 10.2.1 ends at
    # "Click Run and see the results"); the oracle is the TBP-implied cut.
    dry_naphtha = (fs.streams["naphtha"].molar_flow
                   * (1.0 - fs.streams["naphtha"].z.get("water", 0.0)))
    print(f"\n  naphtha (dry) yield  {dry_naphtha:6.1f} mol/s  vs "
          f"TBP-implied {naph:6.1f} mol/s   "
          f"(delta {(dry_naphtha - naph) / naph * 100:+.1f}%)")

    total_in = F + steam
    total_out = sum(fs.streams[s].molar_flow
                    for s in ("naphtha", "kerosene", "diesel", "residue"))
    print(f"  mass balance   in {total_in:.4f} = out {total_out:.4f} mol/s  "
          f"(closes machine-exact)")
    # heaviest cut recovery to residue (the rest goes to the deepest draw)
    heavy = sorted((c for c in fs.component_ids if c.startswith("NBP_")),
                   key=molar_mass)[-1]
    heavy_in = F * crude.z.get(heavy, 0.0)
    heavy_resid = (fs.streams["residue"].molar_flow
                   * fs.streams["residue"].z.get(heavy, 0.0))
    print(f"  heaviest cut ({heavy}) to residue: "
          f"{heavy_resid / heavy_in * 100:.0f}%")
    print(f"\n  total wall time {time.perf_counter() - t_all:.1f} s")


if __name__ == "__main__":
    main()
