"""M1 demo: a flash drum with a liquid recycle.

    FEED ─▶ MIX ─▶ FLASH ─▶ VAP (product)
            ▲        │
            │        ▼
            └── SP ◀ LIQ
                 │
                 ▼
              BOTTOMS (product)

The flash liquid is split: part recycles to the mixer, part leaves as bottoms.
The recycle stream is torn automatically and converged with Wegstein-accelerated
direct substitution. Run from the repo root:

    python examples/02_flash_recycle.py
"""
import sys
from pathlib import Path

_ENGINE = Path(__file__).resolve().parent.parent / "engine"
if _ENGINE.is_dir() and str(_ENGINE) not in sys.path:
    sys.path.insert(0, str(_ENGINE))

from caldyr.core import Component, Flowsheet  # noqa: E402
from caldyr.unitops import FlashDrum, Mixer, Splitter  # noqa: E402


def build() -> Flowsheet:
    fs = Flowsheet(
        components=[Component("n-pentane"), Component("n-octane")],
        property_package="thermo:PR",
    )
    fs.add(Mixer("MIX", {"dP": 0.0}))
    fs.add(FlashDrum("FL", {"T": 360.0, "P": 101325.0}))
    fs.add(Splitter("SP", {"split": 0.6}))     # 60% of flash liquid recycled

    fs.feed("FEED", "MIX:in1", T=330.0, P=101325.0, molar_flow=10.0,
            z={"n-pentane": 0.5, "n-octane": 0.5})
    fs.connect("MIXOUT", "MIX:out", "FL:in1")
    fs.connect("VAP", "FL:vapor", None)        # vapor product
    fs.connect("LIQ", "FL:liquid", "SP:in1")
    fs.connect("RECY", "SP:out1", "MIX:in2")   # recycle to mixer
    fs.connect("BOTTOMS", "SP:out2", None)     # liquid product
    fs.connect("Q", "FL:duty", None)
    return fs


def stream_table(fs: Flowsheet, ids) -> str:
    head = f"{'stream':>8} {'T/K':>7} {'n/(mol/s)':>10} {'x_C5':>7} {'phase':>7}"
    rows = [head, "-" * len(head)]
    for sid in ids:
        s = fs.streams[sid]
        rows.append(f"{sid:>8} {s.T:>7.2f} {s.molar_flow:>10.4f} "
                    f"{s.z.get('n-pentane', 0.0):>7.4f} {s.phase or '':>7}")
    return "\n".join(rows)


def main() -> None:
    fs = build()
    report = fs.solve(tol=1e-8)

    print(f"Torn stream(s): {report.tear_streams}")
    print(f"Method: {report.method}   converged: {report.converged}   "
          f"iterations: {report.iterations}")
    print("Residual per iteration: "
          + ", ".join(f"{r:.1e}" for r in report.history))
    print()
    print(stream_table(fs, ["FEED", "MIXOUT", "VAP", "LIQ", "RECY", "BOTTOMS"]))

    s = fs.streams
    n_in = s["FEED"].molar_flow
    n_out = s["VAP"].molar_flow + s["BOTTOMS"].molar_flow
    print(f"\nOverall mass balance: in {n_in:.5f} == out {n_out:.5f} mol/s "
          f"(residual {abs(n_in - n_out):.2e})")
    print(f"Flash duty: {report.duties['Q'] / 1e3:.2f} kW")


if __name__ == "__main__":
    main()
