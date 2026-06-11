"""M9 demo: rigorous reactors — a Gibbs reformer and a kinetic PFR.

A miniature hydrogen plant front end:

  1. **Steam-methane reformer** as a GibbsReactor (Cantera Gibbs minimization
     at 1100 K / 1 bar): CH4 + H2O -> CO + 3H2 and the water-gas shift run
     simultaneously, so the full syngas slate (H2/CO/CO2/CH4/H2O) comes out of
     one equilibrium without enumerating reactions.
  2. **Cooler** to low-temperature-shift conditions (473 K).
  3. **Water-gas shift PFR** with power-law kinetics, first order in CO
     (activation energy 47.4 kJ/mol from Choi & Stenger, J. Power Sources 124
     (2003) 432 for Cu/ZnO LTS; the pre-exponential is illustrative, sized for
     ~90% CO conversion here). Treated as irreversible — defensible at 473 K
     where the WGS equilibrium constant is ~230, far from equilibrium at this
     conversion.

    python examples/10_reactors.py
"""
import sys
from pathlib import Path

_ENGINE = Path(__file__).resolve().parent.parent / "engine"
if _ENGINE.is_dir() and str(_ENGINE) not in sys.path:
    sys.path.insert(0, str(_ENGINE))

from caldyr.core import Component, Flowsheet, resolve_component  # noqa: E402
from caldyr.unitops import PFR, GibbsReactor, Heater  # noqa: E402

P = 1e5                  # Pa
WGS = {
    "stoich": {"carbon monoxide": -1, "water": -1,
               "carbon dioxide": 1, "hydrogen": 1},
    "key": "carbon monoxide",
    "k0": 2e4,           # 1/s (illustrative magnitude)
    "Ea": 47.4e3,        # J/mol — Choi & Stenger (2003), Cu/ZnO LTS
}


def build() -> Flowsheet:
    comps = ["methane", "water", "carbon monoxide", "carbon dioxide", "hydrogen"]
    fs = Flowsheet(components=[Component(c) for c in comps],
                   property_package="thermo:PR")
    fs.add(GibbsReactor("SMR", {"T": 1100.0}))
    fs.add(Heater("COOL", {"T_out": 473.15}))
    fs.add(PFR("WGS", {"V": 5.0, "T": 473.15, "reactions": [WGS]}))

    # steam:carbon = 3, natural-gas basis 1 mol/s CH4
    fs.feed("FEED", "SMR:in1", T=800.0, P=P, molar_flow=4.0,
            z={"methane": 0.25, "water": 0.75})
    fs.connect("SYNGAS", "SMR:out", "COOL:in1")
    fs.connect("COLD", "COOL:out", "WGS:in1")
    fs.connect("SHIFTED", "WGS:out", None)
    fs.connect("Q_SMR", "SMR:duty", None)
    fs.connect("Q_COOL", "COOL:duty", None)
    fs.connect("Q_WGS", "WGS:duty", None)
    return fs


def show(fs: Flowsheet, sid: str) -> None:
    s = fs.streams[sid]
    comp = "  ".join(f"{resolve_component(c).formula}={s.z.get(c, 0.0):.4f}"
                     for c in s.components)
    print(f"  {sid:8s} T={s.T:7.1f} K  n={s.molar_flow:6.3f} mol/s   {comp}")


def main() -> None:
    fs = build()
    report = fs.solve()
    assert report.converged

    print("=" * 72)
    print("  STEAM-METHANE REFORMING (Gibbs) + LT WATER-GAS SHIFT (kinetic PFR)")
    print("=" * 72)
    print("\nStreams:")
    for sid in ("FEED", "SYNGAS", "COLD", "SHIFTED"):
        show(fs, sid)

    feed, syn, shift = fs.streams["FEED"], fs.streams["SYNGAS"], fs.streams["SHIFTED"]
    ch4_conv = 1.0 - (syn.molar_flow * syn.z["methane"]) / (feed.molar_flow * feed.z["methane"])
    co_in = syn.molar_flow * syn.z["carbon monoxide"]
    co_out = shift.molar_flow * shift.z["carbon monoxide"]
    h2_out = shift.molar_flow * shift.z["hydrogen"]
    ch4_in = feed.molar_flow * feed.z["methane"]

    print("\nPerformance:")
    print(f"  reformer CH4 conversion ......... {ch4_conv:6.1%}")
    print(f"  shift CO conversion ............. {1 - co_out / co_in:6.1%}")
    print(f"  H2 yield ........................ {h2_out / ch4_in:5.2f} mol H2 / mol CH4")

    print("\nDuties (caldyr enthalpy basis, formation-inclusive):")
    for q in ("Q_SMR", "Q_COOL", "Q_WGS"):
        print(f"  {q:8s} {report.duties[q] / 1e3:10.1f} kW")
    print("\n  (reformer endothermic: Q_SMR > 0; shift exothermic at fixed T: "
          "Q_WGS < 0)")


if __name__ == "__main__":
    main()
