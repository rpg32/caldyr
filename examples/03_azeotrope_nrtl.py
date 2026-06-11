"""Why the property package matters: ethanol + water VLE.

This mixture forms a minimum-boiling azeotrope (~89 mol% ethanol, 78.2 C at
1 atm) — the reason you cannot distil past ~95% ethanol. A cubic EOS (PR) cannot
represent it; the NRTL activity-coefficient package can. This script prints the
bubble-temperature curve for both so the difference is visible.

    python examples/03_azeotrope_nrtl.py
"""
import sys
from pathlib import Path

_ENGINE = Path(__file__).resolve().parent.parent / "engine"
if _ENGINE.is_dir() and str(_ENGINE) not in sys.path:
    sys.path.insert(0, str(_ENGINE))

from caldyr.thermo import make_package  # noqa: E402

P = 101325.0
COMPS = ["ethanol", "water"]


def bubble_curve(spec: str):
    # Scan the interior only: a binary VF flash at an exactly-pure endpoint trips
    # a thermo edge case, and the azeotrope (the point of interest) is interior.
    pp = make_package(spec, COMPS)
    rows = []
    for i in range(1, 20):
        x = i / 20.0
        z = {"ethanol": x, "water": 1.0 - x}
        rows.append((x, pp.bubble_dew(P, z)[0]))
    x_az, t_az = min(rows, key=lambda r: r[1])
    return rows, x_az, t_az


def main() -> None:
    print("Ethanol/water bubble temperature at 1 atm\n")
    print(f"{'x_ethanol':>9} | {'PR  T/K':>9} | {'NRTL T/K':>9}")
    print("-" * 33)
    pr_rows, pr_xaz, pr_taz = bubble_curve("thermo:PR")
    nr_rows, nr_xaz, nr_taz = bubble_curve("thermo:NRTL")
    for (x, tpr), (_, tnr) in zip(pr_rows, nr_rows):
        if abs(x * 10 - round(x * 10)) < 1e-9:   # print every 0.1
            print(f"{x:>9.2f} | {tpr:>9.2f} | {tnr:>9.2f}")

    print()
    print(f"PR   min bubble-T: {pr_taz:.2f} K at x={pr_xaz:.2f} "
          f"(monotone toward pure - no azeotrope)")
    print(f"NRTL min bubble-T: {nr_taz:.2f} K at x={nr_xaz:.2f}  "
          f"(azeotrope; textbook ~351.3 K at x~0.89)")


if __name__ == "__main__":
    main()
