"""M14 demo: petroleum assay characterization -> pseudo-component flowsheet.

The worked crude assay of Hameed, "Chemical Process Simulations using Aspen
HYSYS" (Wiley 2025), sec. 10.2.1 "Refinery of Crude Oil-1" (figs. 10.22-10.31;
the classic HYSYS oil-characterization tutorial data):

  * TBP curve 80-1410 F over 0-98 LV%, molecular-weight curve 68-838 g/mol,
    API-gravity curve 63.28-26.01 API
  * light ends (LV%): i-C4 0.19, n-C4 0.11, i-C5 0.37, n-C5 0.46  (1.13 total)
  * bulk inputs: MW 300, standard density 48.75 API
  * cut scheme: start 100 F; to 800 F in 28 cuts, 800-1200 F in 8 cuts,
    1200-1400 F in 2 cuts -> 38 hypocomponents (book step 20-2b)

`caldyr.assay.characterize_assay` reproduces the 38-cut blend; each cut gets
NBP (mid-volume), SG (book API curve), MW (book MW curve here; Riazi-Daubert
correlation when no curve is given), Tc/Pc (Kesler-Lee 1976), omega
(Lee-Kesler) and ideal-gas Cp (Kesler-Lee) — citations in caldyr/assay.py.
The pseudo-components then run as ordinary species in a PR flowsheet: the
book's column front end (crude at 450 F / 75 psia, fired to 650 F with a
10 psi drop, flashed at ~65 psia tower-feed conditions).

    python examples/19_crude_assay.py
"""
import sys
from pathlib import Path

_ENGINE = Path(__file__).resolve().parent.parent / "engine"
if _ENGINE.is_dir() and str(_ENGINE) not in sys.path:
    sys.path.insert(0, str(_ENGINE))

from caldyr.assay import api_to_sg, characterize_assay  # noqa: E402
from caldyr.core import Flowsheet  # noqa: E402
from caldyr.unitops import FlashDrum, Heater  # noqa: E402


def degF(f: float) -> float:
    return (f - 32.0) / 1.8 + 273.15


def toF(k: float) -> float:
    return k * 1.8 - 459.67


# -- the book's assay data (Hameed 2025 figs. 10.22-10.27) ---------------------
TBP_CURVE = [(0, 80.0), (10, 255.0), (20, 349.0), (30, 430.0), (40, 527.0),
             (50, 635.0), (60, 751.0), (70, 915.0), (80, 1095.0),
             (90, 1277.0), (98, 1410.0)]                       # (LV%, F)
MW_CURVE = [(0, 68.0), (10, 119.0), (20, 150.0), (30, 182.0), (40, 225.0),
            (50, 282.0), (60, 350.0), (70, 456.0), (80, 585.0),
            (90, 713.0), (98, 838.0)]                          # (LV%, g/mol)
API_CURVE = [(13, 63.28), (33, 54.86), (57, 45.91), (74, 38.21), (91, 26.01)]
LIGHT_ENDS = {"isobutane": 0.19, "n-butane": 0.11,
              "isopentane": 0.37, "n-pentane": 0.46}           # LV%
BULK_API = 48.75                                               # book fig. 10.23
BULK_MW = 300.0                                                # book fig. 10.23
CUT_POINTS_F = ([100.0 + 25.0 * i for i in range(29)]          # book fig. 10.31
                + [800.0 + 50.0 * i for i in range(1, 9)]
                + [1200.0 + 100.0 * i for i in range(1, 3)])


def main() -> None:
    result = characterize_assay(
        [(v, degF(f)) for v, f in TBP_CURVE],
        kind="TBP",
        api_gravity=BULK_API,
        sg_curve=[(v, api_to_sg(a)) for v, a in API_CURVE],
        mw_curve=MW_CURVE,
        cut_points=[degF(f) for f in CUT_POINTS_F],
        light_ends=LIGHT_ENDS,
    )

    print("=" * 78)
    print("  CRUDE ASSAY CHARACTERIZATION - Hameed (2025) sec. 10.2.1")
    print("=" * 78)
    print(f"\nPseudo-components generated: {len(result.cuts)}   "
          f"(book: 28 + 8 + 2 = 38 hypocomponents)")

    print(f"\n{'cut':<10}{'range F':>16}{'NBP F':>9}{'SG':>8}{'API':>7}"
          f"{'MW':>7}{'Tc K':>8}{'Pc bar':>8}{'omega':>7}{'LV%':>7}")
    for c in result.cuts:
        print(f"{c.id:<10}{toF(c.T_lo):>8.0f}-{toF(c.T_hi):<7.0f}"
              f"{toF(c.Tb):>9.1f}{c.SG:>8.4f}{141.5 / c.SG - 131.5:>7.2f}"
              f"{c.MW * 1000:>7.1f}{c.Tc:>8.1f}{c.Pc / 1e5:>8.2f}"
              f"{c.omega:>7.3f}{c.vol_frac * 100:>7.2f}")

    print("\nLight ends (book fig. 10.22, LV% of assay):")
    for cid, fr in result.light_end_fractions.items():
        print(f"  {cid:<14} LV {fr['vol_frac'] * 100:5.2f}%   "
              f"mol {fr['mole_frac'] * 100:5.2f}%")

    print("\nBulk properties of the characterized blend vs the book's inputs:")
    print(f"  API gravity ....... {result.bulk['API']:6.2f}   (book bulk input "
          f"{BULK_API}; Delta {result.bulk['API'] - BULK_API:+.2f} - the book's "
          f"curves and bulk")
    print("                       value are slightly inconsistent; HYSYS forces "
          "the bulk to win)")
    print(f"  molar mass ........ {result.bulk['MW'] * 1000:6.1f}   (book bulk "
          f"input {BULK_MW:.0f} g/mol; the book's own MW + density")
    print("                       curves imply ~234 - reported honestly, not "
          "rescaled)")
    print(f"  Watson K .......... {result.bulk['watson_k']:6.2f}   "
          f"(intermediate/paraffinic crude)")

    # -- the book's column front end (steps 27-29) ----------------------------
    print("\n" + "=" * 78)
    print("  FRONT END: 450 F crude -> furnace 650 F (-10 psi) -> flash at 65 psia")
    print("=" * 78)

    fs = Flowsheet(components=result.components(), property_package="thermo:PR")
    fs.add(Heater("FURNACE", {"T_out": degF(650.0), "dP": 68947.6}))
    fs.add(FlashDrum("FLASH", {"P": 448159.0}))
    fs.connect("crude", None, "FURNACE:in1")
    fs.connect("hot", "FURNACE:out", "FLASH:in1")
    fs.connect("vap", "FLASH:vapor", None)
    fs.connect("liq", "FLASH:liquid", None)
    fs.feed("crude", "FURNACE:in1", T=degF(450.0), P=517107.0,
            molar_flow=100.0, z=result.mole_fractions())
    report = fs.solve()
    assert report.converged

    vap, liq = fs.streams["vap"], fs.streams["liq"]
    print(f"\n  feed     100.00 mol/s at {toF(fs.streams['crude'].T):>6.1f} F")
    print(f"  vapor   {vap.molar_flow:>7.2f} mol/s ({vap.molar_flow:.1f}% of "
          f"feed) - tower feed vapor")
    print(f"  liquid  {liq.molar_flow:>7.2f} mol/s - atmospheric-resid-rich")

    ids = [c.id for c in result.cuts]
    show = [ids[0], ids[len(ids) // 2], ids[-1]]
    print("\n  split sanity (K = y/x falls monotonically with NBP):")
    for cid in show:
        k = (vap.z.get(cid, 0.0) + 1e-300) / (liq.z.get(cid, 0.0) + 1e-300)
        print(f"    {cid:<10} K = {k:9.3g}")
    heavy = ids[-1]
    rec = liq.molar_flow * liq.z[heavy] / (100.0 * fs.streams["crude"].z[heavy])
    print(f"  heaviest cut to liquid: {rec * 100:.2f}%   "
          f"(the resid stays in the furnace liquid)")


if __name__ == "__main__":
    main()
