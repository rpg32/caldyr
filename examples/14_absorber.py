"""Gas absorption columns: the sum-rates Absorber and the ReboiledAbsorber.

Two worked examples from Hameed, *Chemical Process Simulations using Aspen
Hysys* (Wiley, 2025), ch. 9:

1. **The SO2 absorber of sec. 9.1.1** — 206 kmol/h of 3 mol% SO2 in air
   washed counter-currently with 1.3e5 kg/h of pure water, both at 20 C and
   1 atm, 20 theoretical stages, Peng-Robinson.

   Book (HYSYS-PR) answer: Gas_Out = 204.6 kmol/h at 20.36 C carrying
   0.0576 kmol/h SO2 (99.07% removed); packed-column design D = 1.285 m.

   Honest caveat, stated up front: HYSYS's PR ships a built-in SO2/H2O
   binary interaction parameter giving K_SO2 ~ 30 at 20 C; stock `thermo`
   PR has no kij for the pair (SO2 dissolving in water is chemical —
   hydrolysis to H2SO3 — which no plain cubic EOS represents) and predicts
   K_SO2 ~ 500, so only ~7% of the SO2 absorbs here. The *machinery* is
   validated quantitatively against the Kremser equation — the book's own
   eq. (9.1) — printed below, and against the book's reboiled absorber
   (case 2) where PR thermodynamics is excellent.

2. **The reboiled absorber (stripping tower) of sec. 9.3.5** — 400 kmol/h
   of saturated liquid 55/45 n-pentane/n-heptane, top at 101.3 kPa,
   reboiler at 110 kPa, 8 stages + reboiler, overhead rate 220 kmol/h.

   Book answer: x(nC5) in the overhead vapor = 0.899, x(nC7) in the
   bottoms = 0.876, boilup ratio ~ 1.093 — reproduced here to ~0.1%.

    python examples/14_absorber.py
"""
import math
import sys
from pathlib import Path

_ENGINE = Path(__file__).resolve().parent.parent / "engine"
if _ENGINE.is_dir() and str(_ENGINE) not in sys.path:
    sys.path.insert(0, str(_ENGINE))

from caldyr.core import Component, Flowsheet  # noqa: E402
from caldyr.economics.sizing import SizingOptions, size_flowsheet  # noqa: E402
from caldyr.thermo import make_package  # noqa: E402
from caldyr.unitops import Absorber, ReboiledAbsorber  # noqa: E402

P_ATM = 101325.0
KMOLH = 1000.0 / 3600.0          # kmol/h -> mol/s


def profile_table(design, comps) -> str:
    rows = [f"{'stage':>5} {'T [C]':>8} {'L':>9} {'V':>9} "
            + " ".join(f"x({c[:6]})" for c in comps)]
    for j in range(design["n_stages"]):
        x = design["x_profile"][j]
        rows.append(
            f"{j + 1:>5} {design['T_profile'][j] - 273.15:>8.2f} "
            f"{design['L_profile'][j]:>9.1f} {design['V_profile'][j]:>9.1f} "
            + " ".join(f"{x.get(c, 0.0):>8.5f}" for c in comps))
    return "\n".join(rows)


def so2_absorber() -> None:
    print("=" * 76)
    print("1. SO2 absorber -- Hameed 2025 sec. 9.1.1 (book: 99.07% removal,")
    print("   Gas_Out 204.6 kmol/h at 20.36 C, 0.0576 kmol/h SO2; D = 1.285 m)")
    print("=" * 76)
    fs = Flowsheet(components=[Component("water"), Component("nitrogen"),
                               Component("oxygen"),
                               Component("sulfur dioxide")],
                   property_package="thermo:PR")
    fs.add(Absorber("ABS", {"n_stages": 20, "P": P_ATM}))
    # 'air' is carried as N2/O2 79/21 (not a resolvable pure species).
    fs.feed("GAS", "ABS:gas_in", T=293.15, P=P_ATM, molar_flow=206.0 * KMOLH,
            z={"sulfur dioxide": 0.03, "nitrogen": 0.97 * 0.79,
               "oxygen": 0.97 * 0.21})
    fs.feed("WATER", "ABS:liquid_in", T=293.15, P=P_ATM,
            molar_flow=1.3e5 / 18.01528 * KMOLH, z={"water": 1.0})
    fs.connect("GASOUT", "ABS:vapor_out", None)
    fs.connect("LIQOUT", "ABS:liquid_out", None)
    rep = fs.solve()
    assert rep.converged

    d = fs.units["ABS"].design
    gas_out = fs.streams["GASOUT"]
    so2_out = gas_out.molar_flow * gas_out.z["sulfur dioxide"] / KMOLH
    print(f"\nGas out : {gas_out.molar_flow / KMOLH:7.1f} kmol/h at "
          f"{gas_out.T - 273.15:5.2f} C   (book: 204.6 kmol/h at 20.36 C)")
    print(f"SO2 slip: {so2_out:7.4f} kmol/h "
          f"({d['absorbed']['sulfur dioxide']:.1%} absorbed)   "
          f"(book/HYSYS-PR: 0.0576 kmol/h, 99.07%)")
    print("  -> stock PR has no SO2/H2O binary (K_SO2 ~ 500 vs HYSYS ~ 30):")
    print("     the slip is K-limited, not a column-model difference -- see")
    print("     the Kremser check below for the quantitative validation.")

    print("\nConverged stage profiles (sum-rates / Burningham-Otto, "
          f"{d['iterations']} iterations):")
    print(profile_table(d, ["water", "sulfur dioxide"]))

    # Tray hydraulic design (Fair flooding at 80%) vs the book's packed run.
    pp = make_package("thermo:PR", fs.component_ids)
    tower = next(s for s in size_flowsheet(fs, rep, pp, SizingOptions())
                 if s.unit_id == "ABS")
    print(f"\nTray design: D = {tower.diameter_m:.2f} m "
          f"(book packed design at 80% capacity: 1.285 m)")
    for note in tower.notes:
        print(f"  {note}")


def kremser_check() -> None:
    print()
    print("=" * 76)
    print("Kremser cross-check -- the book's eq. (9.1): fraction absorbed")
    print("  = (A^(N+1) - A)/(A^(N+1) - 1),  A = L/(K V)")
    print("=" * 76)
    n = 6
    fs = Flowsheet(components=[Component("nitrogen"), Component("n-pentane"),
                               Component("n-decane")],
                   property_package="thermo:PR")
    fs.add(Absorber("ABS", {"n_stages": n, "P": P_ATM}))
    fs.feed("GAS", "ABS:gas_in", T=298.15, P=P_ATM, molar_flow=100.0,
            z={"nitrogen": 0.99, "n-pentane": 0.01})
    fs.feed("OIL", "ABS:liquid_in", T=298.15, P=P_ATM, molar_flow=40.0,
            z={"n-decane": 1.0})
    fs.connect("GASOUT", "ABS:vapor_out", None)
    fs.connect("LIQOUT", "ABS:liquid_out", None)
    assert fs.solve().converged
    d = fs.units["ABS"].design
    pp = make_package("thermo:PR", fs.component_ids)
    a = []
    for j in (0, n - 1):
        K = pp.k_values(d["T_profile"][j], d["P_profile"][j],
                        d["x_profile"][j], d["y_profile"][j])["n-pentane"]
        a.append(d["L_profile"][j] / (K * d["V_profile"][j]))
    ae = math.sqrt(a[0] * a[1])
    kremser = (ae ** (n + 1) - ae) / (ae ** (n + 1) - 1.0)
    mesh = d["absorbed"]["n-pentane"]
    print(f"n-pentane from N2 into n-decane, N={n}, A_eff={ae:.3f}:")
    print(f"  MESH    : {mesh:.4f} absorbed")
    print(f"  Kremser : {kremser:.4f} absorbed   "
          f"(delta {abs(mesh - kremser) / kremser:.1%})")


def reboiled_absorber() -> None:
    print()
    print("=" * 76)
    print("2. Reboiled absorber -- Hameed 2025 sec. 9.3.5 (book: overhead")
    print("   220 kmol/h with x(nC5)=0.899; bottoms x(nC7)=0.876; boilup 1.093)")
    print("=" * 76)
    pp = make_package("thermo:PR", ["n-pentane", "n-heptane"])
    z = {"n-pentane": 0.55, "n-heptane": 0.45}
    bub, _ = pp.bubble_dew(110000.0, z)          # saturated-liquid feed

    fs = Flowsheet(components=[Component("n-pentane"),
                               Component("n-heptane")],
                   property_package="thermo:PR")
    fs.add(ReboiledAbsorber("RA", {"n_stages": 9, "P": 101300.0,
                                   "P_bottom": 110000.0,
                                   "vapor_rate": 220.0 * KMOLH}))
    fs.feed("FEED", "RA:feed", T=bub, P=110000.0, molar_flow=400.0 * KMOLH,
            z=z)
    fs.connect("VAP", "RA:vapor_out", None)
    fs.connect("BOT", "RA:bottoms", None)
    fs.connect("QR", "RA:reboiler_duty", None)
    rep = fs.solve()
    assert rep.converged

    vap, bot = fs.streams["VAP"], fs.streams["BOT"]
    d = fs.units["RA"].design
    print(f"\nOverhead: {vap.molar_flow / KMOLH:6.1f} kmol/h, "
          f"x(nC5) = {vap.z['n-pentane']:.4f} at {vap.T - 273.15:.2f} C "
          f"(book: 0.899)")
    print(f"Bottoms : {bot.molar_flow / KMOLH:6.1f} kmol/h, "
          f"x(nC7) = {bot.z['n-heptane']:.4f} at {bot.T - 273.15:.2f} C "
          f"(book: 0.876)")
    print(f"Boilup ratio = {d['boilup_ratio']:.3f} (book ~1.093); "
          f"Q_reb = {d['Q_reboiler'] / 1e6:.2f} MW")
    print("\nStage profiles (bubble-point inner loop, stage 9 = reboiler):")
    print(profile_table(d, ["n-pentane", "n-heptane"]))


if __name__ == "__main__":
    so2_absorber()
    kremser_check()
    reboiled_absorber()
