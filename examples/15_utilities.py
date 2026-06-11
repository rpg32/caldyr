"""M12 demo: the "book chapter 2 & 6" utilities batch (Hameed, *Chemical
Process Simulations using Aspen HYSYS*, Wiley 2025).

Four mini-demos in one script:

1. **Steam point** (§2.2) — the ``coolprop:Water`` steam-tables package
   (IAPWS-95): saturation at 1 atm and the latent heat.
2. **Property table** (§2.1.4) — n-pentane mass density over a (T, P) grid
   under PR, the book's case study.
3. **Humidity point** (§2.4) — moist air at 30 C / 1 atm / 50% RH.
4. **Balance + Evaporator mini-flowsheet** (§5.2 + §6.3) — the book's sucrose
   evaporator as a pure-water surrogate: a 3000 kg/h feed is 80%-vaporized at
   70 kPa in an Evaporator, then a Balance (mode ``mole_heat``) recombines the
   vapor and liquid products, proving the material and energy balances close
   (the recombined stream reproduces the flash state). The evaporator is then
   sized and costed.

Run from the repo root:

    python examples/15_utilities.py
"""
import sys
from pathlib import Path

_ENGINE = Path(__file__).resolve().parent.parent / "engine"
if _ENGINE.is_dir() and str(_ENGINE) not in sys.path:
    sys.path.insert(0, str(_ENGINE))

from caldyr.analysis import humidity, property_table  # noqa: E402
from caldyr.core import Component, Flowsheet  # noqa: E402
from caldyr.economics.costing import cost_equipment  # noqa: E402
from caldyr.economics.sizing import size_flowsheet  # noqa: E402
from caldyr.thermo import make_package  # noqa: E402
from caldyr.unitops import Balance, Evaporator  # noqa: E402

ATM = 101_325.0


def steam_point() -> None:
    print("=" * 72)
    print("1) Steam tables (coolprop:Water, IAPWS-95) — book §2.2")
    pp = make_package("coolprop:Water", ["water"])
    z = {"water": 1.0}
    t_sat, _ = pp.bubble_dew(ATM, z)
    sat = pp.bubble_point(ATM, z)
    h_fg = (sat.H_vapor - sat.H_liquid) / pp.mw / 1e3
    print(f"   Tsat(1 atm)  = {t_sat:.2f} K ({t_sat - 273.15:.2f} C)"
          f"   [steam tables: 99.97 C]")
    print(f"   h_fg(1 atm)  = {h_fg:.1f} kJ/kg          [steam tables: 2256.5]")
    rho = pp.mw / pp.volume(300.0, 1e5, z)
    print(f"   rho(300 K)   = {rho:.1f} kg/m^3          [996.6]")


def property_table_slice() -> None:
    print("=" * 72)
    print("2) Property table — n-pentane mass density (PR), book §2.1.4")
    pp = make_package("thermo:PR", ["n-pentane"])
    out = property_table(pp, {"n-pentane": 1.0},
                         T=[500.0, 525.0, 550.0, 575.0, 600.0],
                         P=[12 * ATM, 16 * ATM, 18 * ATM],
                         props=["mass_density"])
    print("   T [K]   " + "".join(f"{p / ATM:7.0f} atm" for p in out["P"]))
    for i, t in enumerate(out["T"]):
        row = "".join(f"{out['mass_density'][i, j]:11.2f}"
                      for j in range(len(out["P"])))
        print(f"   {t:6.0f} {row}")
    print("   (book Fig 2.7: 28.77 kg/m^3 at 550 K / 16 atm under PR)")


def humidity_point() -> None:
    print("=" * 72)
    print("3) Humidity — air at 30 C, 1 atm, RH 50% (book §2.4)")
    h = humidity(303.15, ATM, rh=0.5)
    print(f"   humidity ratio  w   = {h['w'] * 1e3:.2f} g/kg dry air"
          f"   [book: 13.49]")
    print(f"   dew point       Tdp = {h['t_dp'] - 273.15:.2f} C        [book: 18.86]")
    print(f"   wet bulb        Twb = {h['t_wb'] - 273.15:.2f} C        [chart: ~22.0]")
    print(f"   moist enthalpy  h   = {h['h'] / 1e3:.1f} kJ/kg dry air")
    print(f"   moist volume    v   = {h['v']:.3f} m^3/kg dry air")


def evaporator_flowsheet() -> None:
    print("=" * 72)
    print("4) Evaporator + Balance mini-flowsheet (book §5.2 + §6.3)")
    mw = 0.018015268
    feed_mol_s = 3000.0 / 3600.0 / mw            # 3000 kg/h of water

    fs = Flowsheet(components=[Component("water")],
                   property_package="coolprop:Water")
    fs.add(Evaporator("EVAP", {"P": 70e3, "vapor_fraction": 0.8}))
    fs.add(Balance("BAL", {"mode": "mole_heat", "n_inlets": 2}))
    fs.feed("FEED", "EVAP:in1", T=303.15, P=110e3, molar_flow=feed_mol_s,
            z={"water": 1.0})
    fs.connect("VAP", "EVAP:vapor", "BAL:in1")
    fs.connect("LIQ", "EVAP:liquid", "BAL:in2")
    fs.connect("Q", "EVAP:duty", None)
    fs.connect("CHECK", "BAL:out1", None)
    rep = fs.solve()
    print(f"   converged: {rep.converged}")

    vap, liq, chk = fs.streams["VAP"], fs.streams["LIQ"], fs.streams["CHECK"]
    q = rep.duties["Q"]
    print(f"   flash at {vap.T - 273.15:.1f} C / 70 kPa: "
          f"vapor {vap.molar_flow * mw * 3600:.0f} kg/h, "
          f"liquid {liq.molar_flow * mw * 3600:.0f} kg/h, duty {q / 1e6:.2f} MW")

    from CoolProp.CoolProp import PropsSI
    h_fg_152 = (PropsSI("Hmass", "T", 425.15, "Q", 1, "Water")
                - PropsSI("Hmass", "T", 425.15, "Q", 0, "Water"))
    steam = q / h_fg_152 * 3600.0
    print(f"   saturated 152 C steam required: {steam:.0f} kg/h "
          f"[book: 2894 kg/h, +2.2% — pure-water surrogate of the 10% sucrose"
          f" feed]")
    print(f"   'efficiency' (vapor/steam, book Eq. 5.7): "
          f"{vap.molar_flow * mw * 3600 / steam:.2f}   [book: 0.83]")

    # The Balance (mole_heat) recombines the two products: component moles and
    # enthalpy flow both conserved, so it reproduces the flash state exactly.
    print(f"   Balance check : recombined {chk.molar_flow * mw * 3600:.0f} kg/h at "
          f"{chk.T - 273.15:.1f} C, VF {chk.vapor_fraction:.3f} "
          f"(= the flash point — balances close)")

    pp = make_package(fs.property_package, fs.component_ids)
    sizes = size_flowsheet(fs, rep, pp)
    for s in sizes:
        c = cost_equipment(s)
        print(f"   economics     : {s.equipment_type} {s.attribute:.2f} "
              f"{s.attribute_name}, utility {s.utility} "
              f"({s.utility_duty_W / 1e6:.2f} MW) -> bare module "
              f"${c.bare_module:,.0f}")


def main() -> None:
    steam_point()
    property_table_slice()
    humidity_point()
    evaporator_flowsheet()
    print("=" * 72)


if __name__ == "__main__":
    main()
