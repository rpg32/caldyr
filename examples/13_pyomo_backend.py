"""Pyomo grey-box EO backend demo.

Solves the flash-with-recycle flowsheet (the M4 example) with
``backend="pyomo"`` — the flowsheet residual system wrapped as a PyNumero
ExternalGreyBoxBlock and solved with IPOPT through cyipopt — and compares it
against the sequential-modular solution. Where no grey-box-capable NLP solver
is installed (cyipopt; pip-only Windows usually cannot build it), the demo
prints the backend's install guidance and still shows that the Pyomo model
*constructs*: variables, constraints, and an evaluable residual/Jacobian.

    python examples/13_pyomo_backend.py
"""
import sys
from pathlib import Path

_ENGINE = Path(__file__).resolve().parent.parent / "engine"
if _ENGINE.is_dir() and str(_ENGINE) not in sys.path:
    sys.path.insert(0, str(_ENGINE))

from caldyr.core import Component, Flowsheet  # noqa: E402
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


def main() -> None:
    fs_seq = build()
    rep_seq = fs_seq.solve(backend="sequential", tol=1e-9)
    print("sequential-modular reference:")
    print(f"  converged={rep_seq.converged}  iterations={rep_seq.iterations}  "
          f"VAP={fs_seq.streams['VAP'].molar_flow:.5f} mol/s  "
          f"RECY={fs_seq.streams['RECY'].molar_flow:.5f} mol/s")

    from caldyr.solver.pyomo_backend import solver_unavailable_reason
    reason = solver_unavailable_reason()
    if reason is None:
        fs = build()
        rep = fs.solve(backend="pyomo", tol=1e-8)
        print("\npyomo grey-box backend (IPOPT via cyipopt):")
        print(f"  converged={rep.converged}  residual={rep.residual:.2e}  "
              f"method={rep.method}")
        for msg in rep.messages:
            print(f"  {msg}")
        diff = max(abs(fs.streams[s].molar_flow - fs_seq.streams[s].molar_flow)
                   for s in ("VAP", "BOT", "RECY"))
        print(f"  -> max stream-flow difference vs sequential: {diff:.2e} mol/s")
        return

    print("\npyomo backend: live solve unavailable on this machine -")
    print(f"  {reason}")
    try:
        import pyomo  # noqa: F401
    except ImportError:
        return

    # pyomo itself is here: show that the translation layer works regardless.
    from caldyr.solver import PyomoEOSolver
    from caldyr.thermo import make_package

    fs = build()
    pp = make_package(fs.property_package, fs.component_ids)
    model, system = PyomoEOSolver().build_model(fs, pp)
    egb = model.residuals.get_external_model()
    egb.set_input_values(system.x0)
    r = egb.evaluate_equality_constraints()
    print("\nbut the Pyomo model constructs and evaluates fine:")
    print(f"  {egb.n_inputs()} variables, {egb.n_equality_constraints()} grey-box "
          f"equality constraints (5 streams x [n_C5, n_C8, T, P])")
    print(f"  max |residual| at the warm start: {max(abs(v) for v in r):.3e} "
          f"(nonzero -> the recycle is not yet converged; IPOPT would drive it to 0)")


if __name__ == "__main__":
    main()
