"""Liquid-liquid extraction and extractive distillation.

Two worked examples from Hameed, *Chemical Process Simulations using Aspen
Hysys* (Wiley, 2025), ch. 9:

1. **The acetone extractor of sec. 9.4** — 100 kg/h of 50/50 wt%
   acetone/water (the heavy liquid, fed at the top) washed counter-currently
   by 150 kg/h of 3-methylhexane rising from the bottom; 8 trays, 1 bar,
   25 C.

   Book (HYSYS-PRSV) answer: the bottoms heavy-liquid Water stream leaves
   with a water mole fraction of 1.0 (sec. 9.4.2 step 9) — the acetone goes
   overhead with the solvent as the Rich-Sol stream. Reproduced here with
   stock-`thermo` PR (which behaves like the book's PRSV on this point: a
   cubic EOS expels essentially all non-water species from the aqueous
   phase). Honest caveat: that makes the acetone transfer *total* under
   either EOS, so this system exercises the structure, not NRTL-grade LLE
   partitioning — the ExtractionColumn's quantitative validation is the
   Kremser cross-check in tests/test_m13_extraction.py (Seader 3e ch. 8).

2. **Extractive distillation, the sec. 9.5.5 structure** — the book feeds
   phenol *above* the main feed of a 50-stage NRTL column to split
   close-boiling n-heptane/toluene, cutting total heating duty from
   1.909e7 kJ/h (one 80-stage column, no solvent) to 7.549e6 kJ/h. The
   bundled ChemSep NRTL table carries no parameters for any
   heptane/toluene/phenol pair (an ideal-liquid fallback would show no
   extractive effect at all), so the same two-feed structure is demonstrated
   on the classic acetone/methanol azeotrope with water as the heavy solvent
   (Seader 3e ch. 11) — all three binaries parameterized: the same column
   that pins below the ~79 mol% azeotrope without solvent passes 88 mol%
   acetone with it.

    python examples/17_extraction.py
"""
import sys
from pathlib import Path

_ENGINE = Path(__file__).resolve().parent.parent / "engine"
if _ENGINE.is_dir() and str(_ENGINE) not in sys.path:
    sys.path.insert(0, str(_ENGINE))

from caldyr.core import Component, Flowsheet  # noqa: E402
from caldyr.economics.sizing import SizingOptions, size_flowsheet  # noqa: E402
from caldyr.thermo import make_package  # noqa: E402
from caldyr.unitops import ExtractionColumn, RigorousColumn  # noqa: E402

P_ATM = 101325.0
KMOLH = 1000.0 / 3600.0          # kmol/h -> mol/s
MW = {"acetone": 58.0791, "water": 18.01528, "3-methylhexane": 100.2019}


def acetone_extractor() -> None:
    print("=" * 76)
    print("1. Acetone extractor -- Hameed 2025 sec. 9.4 (book: bottoms Water")
    print("   stream x(H2O) = 1.0; acetone leaves overhead in Rich-Sol)")
    print("=" * 76)
    f_acetone = 50.0 / MW["acetone"] * KMOLH
    f_water = 50.0 / MW["water"] * KMOLH

    fs = Flowsheet(components=[Component("acetone"), Component("water"),
                               Component("3-methylhexane")],
                   property_package="thermo:PR")
    fs.add(ExtractionColumn("EXT", {"n_stages": 8, "T": 298.15, "P": 1e5}))
    # The book's orientation rule (sec. 9.4.2 step 4): the heavier liquid
    # (the aqueous feed) enters at the top, the lighter (the alkane solvent)
    # at the bottom.
    fs.feed("FEED", "EXT:feed_in", T=298.15, P=1e5,
            molar_flow=f_acetone + f_water,
            z={"acetone": f_acetone / (f_acetone + f_water),
               "water": f_water / (f_acetone + f_water)})
    fs.feed("SOLV", "EXT:solvent_in", T=298.15, P=1e5,
            molar_flow=150.0 / MW["3-methylhexane"] * KMOLH,
            z={"3-methylhexane": 1.0})
    fs.connect("EXTRACT", "EXT:extract_out", None)
    fs.connect("RAFF", "EXT:raffinate_out", None)
    fs.connect("Q", "EXT:duty", None)
    rep = fs.solve()
    assert rep.converged

    e, r = fs.streams["EXTRACT"], fs.streams["RAFF"]
    d = fs.units["EXT"].design
    print(f"\nExtract (top, Rich-Sol): {e.molar_flow / KMOLH:6.3f} kmol/h   "
          f"x = {{acetone {e.z['acetone']:.4f}, water {e.z['water']:.4f}, "
          f"3MH {e.z['3-methylhexane']:.4f}}}")
    print(f"Raffinate (bottom, Water): {r.molar_flow / KMOLH:5.3f} kmol/h   "
          f"x(water) = {r.z['water']:.4f}   (book: 1.0)")
    print(f"Acetone recovery to the extract: "
          f"{d['recovery']['acetone']:.4f}")
    print(f"Converged in {d['iterations']} sweeps "
          f"({d['flash_calls']} LLE flashes); isothermal duty "
          f"{d['duty'] / 1e3:.2f} kW")

    print("\nStage profiles (ascending E / descending R, top -> bottom):")
    print(f"{'stage':>5} {'E':>8} {'R':>8} {'x_E(acetone)':>13} "
          f"{'x_R(acetone)':>13}")
    for j in range(d["n_stages"]):
        print(f"{j + 1:>5} {d['E_profile'][j]:>8.3f} "
              f"{d['R_profile'][j]:>8.3f} "
              f"{d['x_extract_profile'][j]['acetone']:>13.5f} "
              f"{d['x_raffinate_profile'][j]['acetone']:>13.5f}")

    pp = make_package("thermo:PR", fs.component_ids)
    sizes = [s for s in size_flowsheet(fs, rep, pp, SizingOptions())
             if s.unit_id.startswith("EXT")]
    print("\nTower sizing (no condenser/reboiler -- liquid throughput sets D):")
    for s in sizes:
        for note in s.notes:
            print(f"  {note}")


def extractive_distillation() -> None:
    print()
    print("=" * 76)
    print("2. Extractive distillation -- the Hameed 2025 sec. 9.5.5 structure")
    print("   (book: phenol above the feed cuts heating from 1.909e7 to")
    print("   7.549e6 kJ/h; demonstrated on acetone/methanol + water because")
    print("   ChemSep NRTL has no heptane/toluene/phenol parameters)")
    print("=" * 76)

    def column(solvent: float | None) -> Flowsheet:
        fs = Flowsheet(components=[Component("acetone"),
                                   Component("methanol"),
                                   Component("water")],
                       property_package="thermo:NRTL")
        params: dict = {"n_stages": 16, "reflux_ratio": 3.0,
                        "distillate_rate": 5.0, "P": P_ATM, "max_iter": 600}
        if solvent is None:
            params["feed_stage"] = 9
        else:
            # The book's two-feed layout: heavy solvent a few stages below
            # the condenser, the main feed mid-column.
            params["feeds"] = [{"stage": 9}, {"stage": 4}]
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

    # Where is the azeotrope under our NRTL?
    pp = make_package("thermo:NRTL", ["acetone", "methanol"])
    lo, hi = 0.6, 0.95
    for _ in range(40):
        mid = 0.5 * (lo + hi)
        res = pp.bubble_point(P_ATM, {"acetone": mid, "methanol": 1.0 - mid})
        lo, hi = (mid, hi) if res.y["acetone"] > mid else (lo, mid)
    azeo = 0.5 * (lo + hi)
    print(f"\nacetone/methanol azeotrope under ChemSep NRTL: "
          f"{azeo:.4f} mol acetone (exp. ~0.78-0.80 at 1 atm)")

    base = column(solvent=None)
    assert base.solve().converged
    xb = base.streams["DIST"].z
    print("\nWithout solvent (16 stages, R=3, D=5 mol/s):")
    print(f"  distillate x(acetone) = {xb['acetone']:.4f}  "
          f"<- pinned below the azeotrope")

    extr = column(solvent=10.0)
    assert extr.solve().converged
    xd, xbot = extr.streams["DIST"].z, extr.streams["BOT"].z
    qc = extr.last_report.duties["QC"]
    qr = extr.last_report.duties["QR"]
    print("\nSame column + 10 mol/s water fed on stage 4 (above the feed):")
    print(f"  distillate x(acetone) = {xd['acetone']:.4f}  "
          f"(methanol {xd['methanol']:.4f}, water {xd['water']:.4f})")
    print(f"  bottoms   x(water)    = {xbot['water']:.4f}  "
          f"<- the book's Rich-Solvent, off to solvent recovery")
    print(f"  duties: condenser {qc / 1e6:.2f} MW, reboiler {qr / 1e6:.2f} MW")
    print("\nThe solvent lifts the distillate past the azeotrope -- the")
    print("book's structural result (its phenol column does the same for")
    print("heptane/toluene at a third of the single-column energy).")


if __name__ == "__main__":
    acetone_extractor()
    extractive_distillation()
