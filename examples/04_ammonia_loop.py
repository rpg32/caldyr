"""M2 demo: an ammonia synthesis (Haber-Bosch) loop.

    MAKEUP ─▶ MIX ─▶ PREHEAT ─▶ REACTOR ─▶ COOL ─▶ SEP ─┬─▶ PRODUCT (liquid NH3)
               ▲                                         │
               │                                  (vapor)│
               └──────────── RECYCLE ◀── SPLIT ◀─────────┘
                                          │
                                          ▼
                                       PURGE (bleeds inert argon)

N2 + 3H2 ⇌ 2NH3 is exothermic and mole-reducing, so equilibrium per pass is low;
unreacted gas is recycled. A small argon inert enters with the makeup and would
accumulate without the purge — so the splitter bleeds a little vapor. The recycle
is torn and converged automatically. Run from the repo root:

    python examples/04_ammonia_loop.py
"""
import sys
from pathlib import Path

_ENGINE = Path(__file__).resolve().parent.parent / "engine"
if _ENGINE.is_dir() and str(_ENGINE) not in sys.path:
    sys.path.insert(0, str(_ENGINE))

from caldyr.core import Component, Flowsheet  # noqa: E402
from caldyr.unitops import (  # noqa: E402
    EquilibriumReactor, FlashDrum, Heater, Mixer, Splitter,
)

P_LOOP = 2.0e7      # 200 bar
T_RXN = 673.15      # 400 C
T_SEP = 250.0       # chilled separator
AMMONIA = {"stoich": {"nitrogen": -1, "hydrogen": -3, "ammonia": 2}, "key": "nitrogen"}


def build() -> Flowsheet:
    fs = Flowsheet(
        components=[Component(c) for c in ("nitrogen", "hydrogen", "ammonia", "argon")],
        property_package="thermo:PR",
    )
    fs.add(Mixer("MIX", {"dP": 0.0}))
    fs.add(Heater("PREHEAT", {"T_out": T_RXN}))
    fs.add(EquilibriumReactor("RXN", {"reaction": AMMONIA, "T": T_RXN}))
    fs.add(Heater("COOL", {"T_out": T_SEP}))
    fs.add(FlashDrum("SEP", {"T": T_SEP, "P": P_LOOP}))
    fs.add(Splitter("SPLIT", {"split": 0.90}))     # 10% of the vapor is purged

    # Makeup: 3:1 H2:N2 with 1% argon inert, at loop pressure.
    fs.feed("MAKEUP", "MIX:in1", T=300.0, P=P_LOOP, molar_flow=100.0,
            z={"nitrogen": 0.2475, "hydrogen": 0.7425, "ammonia": 0.0, "argon": 0.01})
    fs.connect("S1", "MIX:out", "PREHEAT:in1")
    fs.connect("S2", "PREHEAT:out", "RXN:in1")
    fs.connect("S3", "RXN:out", "COOL:in1")
    fs.connect("S4", "COOL:out", "SEP:in1")
    fs.connect("PRODUCT", "SEP:liquid", None)
    fs.connect("VAP", "SEP:vapor", "SPLIT:in1")
    fs.connect("RECYCLE", "SPLIT:out1", "MIX:in2")
    fs.connect("PURGE", "SPLIT:out2", None)
    for unit in ("PREHEAT", "RXN", "COOL", "SEP"):
        fs.connect(f"Q_{unit}", f"{unit}:duty", None)
    return fs


def stream_table(fs: Flowsheet, ids) -> str:
    cols = ("nitrogen", "hydrogen", "ammonia", "argon")
    head = f"{'stream':>8} {'T/K':>7} {'n/(mol/s)':>10} " + " ".join(f"{c[:4]:>6}" for c in cols)
    rows = [head, "-" * len(head)]
    for sid in ids:
        s = fs.streams[sid]
        comps = " ".join(f"{s.z.get(c, 0.0):>6.3f}" for c in cols)
        rows.append(f"{sid:>8} {s.T:>7.1f} {s.molar_flow:>10.3f} {comps}")
    return "\n".join(rows)


def main() -> None:
    fs = build()
    report = fs.solve(tol=1e-7, max_iter=400)

    print(f"Torn stream(s): {report.tear_streams}   method: {report.method}")
    print(f"Converged: {report.converged} in {report.iterations} iterations "
          f"(residual {report.residual:.1e})\n")
    print(stream_table(fs, ["MAKEUP", "S2", "S3", "PRODUCT", "RECYCLE", "PURGE"]))

    s = fs.streams
    nh3 = s["PRODUCT"].molar_flow * s["PRODUCT"].z["ammonia"]
    per_pass = s["S3"].z["ammonia"]
    print(f"\nPer-pass reactor NH3 mole fraction: {per_pass:.3f}")
    print(f"Recycle ratio (recycle / makeup): {s['RECYCLE'].molar_flow / 100.0:.2f}")
    print(f"NH3 product rate: {nh3:.2f} mol/s at {s['PRODUCT'].z['ammonia']*100:.1f}% purity")

    # Atom balances close around the whole loop.
    def atoms(element):
        coeff = {"N": {"nitrogen": 2, "ammonia": 1}, "H": {"hydrogen": 2, "ammonia": 3},
                 "Ar": {"argon": 1}}[element]
        def total(sid):
            return sum(s[sid].molar_flow * s[sid].z.get(c, 0.0) * k for c, k in coeff.items())
        return total("MAKEUP"), total("PRODUCT") + total("PURGE")
    print("\nOverall atom balances (in / out):")
    for el in ("N", "H", "Ar"):
        i, o = atoms(el)
        print(f"  {el:>2}: {i:8.3f} / {o:8.3f}")


if __name__ == "__main__":
    main()
