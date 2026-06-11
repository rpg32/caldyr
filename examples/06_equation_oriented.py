"""M4 demo: equation-oriented solve + optimization.

Part 1 solves the flash-with-recycle flowsheet two ways — sequential-modular
(tear + Wegstein) and equation-oriented (all equations at once, no tearing) — and
shows they agree. Part 2 optimizes a design variable: choose the flash
temperature that *minimizes heating duty* while still recovering a target amount
of pentane overhead.

    python examples/06_equation_oriented.py
"""
import sys
from pathlib import Path

_ENGINE = Path(__file__).resolve().parent.parent / "engine"
if _ENGINE.is_dir() and str(_ENGINE) not in sys.path:
    sys.path.insert(0, str(_ENGINE))

from caldyr.core import Component, Flowsheet  # noqa: E402
from caldyr.solver import DesignVar, optimize  # noqa: E402
from caldyr.unitops import FlashDrum, Mixer, Splitter  # noqa: E402


def build() -> Flowsheet:
    fs = Flowsheet(components=[Component("n-pentane"), Component("n-octane")],
                   property_package="thermo:PR")
    fs.add(Mixer("MIX", {"dP": 0.0}))
    fs.add(FlashDrum("FL", {"T": 360.0, "P": 101325.0}))
    fs.add(Splitter("SP", {"split": 0.6}))
    fs.feed("FEED", "MIX:in1", T=330.0, P=101325.0, molar_flow=10.0,
            z={"n-pentane": 0.5, "n-octane": 0.5})
    fs.connect("MIXOUT", "MIX:out", "FL:in1")
    fs.connect("VAP", "FL:vapor", None)
    fs.connect("LIQ", "FL:liquid", "SP:in1")
    fs.connect("RECY", "SP:out1", "MIX:in2")
    fs.connect("BOT", "SP:out2", None)
    fs.connect("Q", "FL:duty", None)
    return fs


def c5_overhead(fs) -> float:
    v = fs.streams["VAP"]
    return v.molar_flow * v.z["n-pentane"]


def main() -> None:
    print("Part 1 - same flowsheet, two solvers")
    print(f"  {'backend':<20}{'iters/nfev':>11}{'VAP n':>10}{'RECY n':>10}{'residual':>12}")
    results = {}
    for backend in ("sequential", "equation_oriented"):
        fs = build()
        rep = fs.solve(backend=backend, tol=1e-9)
        results[backend] = fs
        print(f"  {backend:<20}{rep.iterations:>11}{fs.streams['VAP'].molar_flow:>10.5f}"
              f"{fs.streams['RECY'].molar_flow:>10.5f}{rep.residual:>12.1e}")
    a, b = results["sequential"], results["equation_oriented"]
    diff = max(abs(a.streams[s].molar_flow - b.streams[s].molar_flow)
               for s in ("VAP", "BOT", "RECY"))
    print(f"  -> max stream-flow difference: {diff:.2e} mol/s (the recycle was solved\n"
          f"     simultaneously, with no tear stream)")

    print("\nPart 2 - optimize flash temperature: min duty s.t. C5 overhead >= 4.2 mol/s")
    fs = build()
    base = fs.solve(tol=1e-9)
    print(f"  baseline T = 360.0 K:  duty = {base.duties['Q']/1e3:6.1f} kW, "
          f"C5 overhead = {c5_overhead(fs):.3f} mol/s  (infeasible)")

    res = optimize(
        fs,
        objective=lambda fs, rep: rep.duties["Q"] / 1e3,            # minimize kW
        design_vars=[DesignVar("FL", "T", 340.0, 370.0, initial=360.0)],
        constraints=[lambda fs, rep: c5_overhead(fs) - 4.2],        # >= 0
        solve_kwargs={"tol": 1e-9},
    )
    t_opt = res.design["FL.T"]
    print(f"  optimum   T = {t_opt:5.1f} K:  duty = {res.objective:6.1f} kW, "
          f"C5 overhead = {c5_overhead(fs):.3f} mol/s  (binding)")
    print(f"  success={res.success}  flowsheet solves={res.n_solves}")


if __name__ == "__main__":
    main()
