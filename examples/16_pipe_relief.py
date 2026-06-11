"""M13 demo: the "book chapter 7" piping & safety batch (Hameed, *Chemical
Process Simulations using Aspen HYSYS*, Wiley 2025).

Two worked examples from the book, reproduced end-to-end:

1. **Pipe network (sec. 7.1, Figs. 7.1/7.13)** — 700 kg/h of water at 60 F and
   1 atm is throttled (VLV-1, dP = 31.7 kPa from the book's Cv sizing), runs
   down 6 ft and across 15 ft of 1-in schedule-40 cast-iron pipe (PIPE-100),
   is pumped by a 2-hp / 50%-efficient pump, then climbs 300 ft through 450 ft
   of 1-in pipe and 375 ft of 0.5-in pipe to point D (PIPE-101). Book results
   (Fig. 7.13, HYSYS + ASME Steam):
       P at point C (pump discharge)  = 568.284 psia = 3918.4 kPa
       P at point D (delivery)        = 410.837 psia = 2832.6 kPa
   This script reproduces both to ~0.1% with PipeSegment (Churchill friction,
   IAPWS-95 water) and answers the book's part (b): the pump hp needed to hold
   point D at 3.2 atm (the book leaves the number to the reader).

2. **PSV sizing (sec. 7.3, Figs. 7.36-7.41)** — a blocked-outlet relief of
   3000 kg/h of steam at 224.5 C, relieving pressure 4.950 barG, backpressure
   0.450 barG, Kd = 0.975. Book results (HYSYS):
       required orifice  = 10.17 cm^2
       selected orifice  = API 526 "K" (11.86 cm^2)
       capacity used     = 85.74%
   Reproduced with `caldyr.analysis.relief_vapor` (API 520 critical vapor
   flow) using IAPWS-95 Z and the ideal-gas k of steam.

Run from the repo root:

    python examples/16_pipe_relief.py
"""
import sys
from pathlib import Path

_ENGINE = Path(__file__).resolve().parent.parent / "engine"
if _ENGINE.is_dir() and str(_ENGINE) not in sys.path:
    sys.path.insert(0, str(_ENGINE))

from caldyr.analysis import relief_vapor  # noqa: E402
from caldyr.core import Component, Flowsheet  # noqa: E402
from caldyr.economics.costing import cost_equipment  # noqa: E402
from caldyr.economics.sizing import size_flowsheet  # noqa: E402
from caldyr.thermo import make_package  # noqa: E402
from caldyr.unitops import PipeSegment, Pump, Valve  # noqa: E402

ATM = 101_325.0
FT = 0.3048
IN = 0.0254
PSI = 6894.757
HP = 745.699872            # W per horsepower
M_H2O = 0.01801528         # kg/mol

T60F = (60.0 - 32.0) / 1.8 + 273.15      # 60 F
N_FEED = (700.0 / 3600.0) / M_H2O        # 700 kg/h water, mol/s

# Schedule-40 cast iron (book Fig. 7.7): ID 1.049 in / 0.622 in, roughness
# 8.497e-4 ft; Crane TP-410 K factors for the fittings of Fig. 7.1.
D_1IN, D_HALF = 1.049 * IN, 0.622 * IN
ROUGH = 8.497e-4 * FT
K_ELBOW_1IN = 30 * 0.023                 # 90-deg std elbow, K = 30 fT
K_GATE_HALF = 160 * 0.027                # half-open gate valve, K = 160 fT
K_UNION = 0.08


def build_pipe_network(pump_p_out: float) -> Flowsheet:
    """The A -> D system of book Fig. 7.1 (valve, PIPE-100, pump, PIPE-101)."""
    fs = Flowsheet(components=[Component("water")],
                   property_package="coolprop:Water")
    # VLV-1: the book sizes the valve from Cv = 2.878 USGPM at 50% opening and
    # HYSYS computes dP = 31.7 kPa (book Fig. 7.3) — used here as the spec.
    fs.add(Valve("VLV1", {"dP": 31.7e3}))
    # PIPE-100: elbow, 6 ft down, elbow, 15 ft horizontal (1-in).
    fs.add(PipeSegment("PIPE100", {
        "length": 21 * FT, "diameter": D_1IN, "roughness": ROUGH,
        "elevation_change": -6 * FT, "fittings_K": 2 * K_ELBOW_1IN,
        "segments": 5}))
    # P-100: the book specs duty = 2 hp at 50% efficiency; our Pump is
    # P_out-specified, so the caller passes the discharge pressure that makes
    # the *fluid* power 1 hp (Eq. 7.2) — the shaft duty then reads back 2 hp.
    fs.add(Pump("P100", {"P_out": pump_p_out, "eta": 0.50}))
    # PIPE-101, 1-in part: 150 ft horizontal + 300 ft riser, 2 elbows + union.
    fs.add(PipeSegment("PIPE101A", {
        "length": 450 * FT, "diameter": D_1IN, "roughness": ROUGH,
        "elevation_change": 300 * FT,
        "fittings_K": 2 * K_ELBOW_1IN + K_UNION, "segments": 15}))
    # PIPE-101, 0.5-in last section: 375 ft + VLV-2 (gate valve, half open).
    fs.add(PipeSegment("PIPE101B", {
        "length": 375 * FT, "diameter": D_HALF, "roughness": ROUGH,
        "fittings_K": K_GATE_HALF, "segments": 13}))

    fs.feed("A", "VLV1:in1", T=T60F, P=ATM, molar_flow=N_FEED, z={"water": 1.0})
    fs.connect("A1", "VLV1:out", "PIPE100:in1")
    fs.connect("B", "PIPE100:out", "P100:in1")
    fs.connect("C", "P100:out", "PIPE101A:in1")
    fs.connect("MID", "PIPE101A:out", "PIPE101B:in1")
    fs.connect("D", "PIPE101B:out", None)
    fs.connect("W_PUMP", "P100:work", None)
    return fs


def pipe_network() -> None:
    print("=" * 72)
    print("1) Pipe network — book sec. 7.1 (700 kg/h water, Fig. 7.1)")
    pp = make_package("coolprop:Water", ["water"])

    # Pass 1 with a placeholder discharge pressure to learn the pump suction
    # state, then set P_out so the fluid power is exactly 1 hp (= 50% of the
    # book's 2-hp driver; Eq. 7.2: Pw = rho g h Q = dP * Q).
    fs = build_pipe_network(pump_p_out=40e5)
    fs.solve()
    b = fs.streams["B"]
    q_b = b.molar_flow * pp.volume(b.T, b.P, b.normalized_z())   # m^3/s
    dp_pump = 1.0 * HP / q_b
    fs = build_pipe_network(pump_p_out=b.P + dp_pump)
    rep = fs.solve()
    assert rep.converged

    p_c, p_d = fs.streams["C"].P, fs.streams["D"].P
    shaft_hp = rep.duties["W_PUMP"] / HP
    print(f"   pump suction (B)   = {fs.streams['B'].P / 1e3:7.1f} kPa")
    print(f"   pump shaft duty    = {shaft_hp:7.2f} hp        [book: 2 hp]")
    print(f"   point C pressure   = {p_c / PSI:7.1f} psia     "
          f"[book Fig. 7.13: 568.284]")
    print(f"   point D pressure   = {p_d / PSI:7.1f} psia     "
          f"[book Fig. 7.13: 410.837]  "
          f"({(p_d / (410.837 * PSI) - 1) * 100:+.2f}%)")

    d1 = fs.units["PIPE101A"].design
    d2 = fs.units["PIPE101B"].design
    print(f"   PIPE-101 breakdown: friction {d1['dP_friction'] / 1e3:.1f} "
          f"(1-in) + {d2['dP_friction'] / 1e3:.1f} (0.5-in) kPa, "
          f"elevation {d1['dP_elevation'] / 1e3:.1f} kPa, "
          f"fittings {(d1['dP_fittings'] + d2['dP_fittings']) / 1e3:.2f} kPa")
    print(f"   regimes: 1-in Re={d1['Re']:.0f} f={d1['friction_factor']:.4f}, "
          f"0.5-in Re={d2['Re']:.0f} f={d2['friction_factor']:.4f} "
          f"(both {d1['flow_regime']})")

    # Part (b): the pump hp that holds point D at 3.2 atm. The flow is fixed,
    # so every loss in the line is unchanged and the pump dP simply shrinks by
    # the excess delivery pressure.
    dp_b = dp_pump - (p_d - 3.2 * ATM)
    print(f"   (b) for P(D) = 3.2 atm: fluid power = {dp_b * q_b / HP:.3f} hp "
          f"-> shaft = {dp_b * q_b / 0.5 / HP:.2f} hp at 50% efficiency")

    # Installed pipe cost (Sinnott C&R Vol. 6 correlation via the sizer).
    sizes = [s for s in size_flowsheet(fs, rep, pp)
             if s.equipment_type == "pipe"]
    total = sum(cost_equipment(s).bare_module for s in sizes)
    print(f"   installed cost of the {sum(1 for _ in sizes)} pipe runs "
          f"(~{(21 + 450 + 375) * FT:.0f} m): ${total:,.0f} (2023)")


def psv_sizing() -> None:
    print("=" * 72)
    print("2) PSV sizing — book sec. 7.3 (blocked-outlet steam relief)")
    # Relieving conditions (book Fig. 7.40): 3000 kg/h steam, T = 224.5 C,
    # P1 = 4.950 barG (set 4.5 barG + 10%), backpressure 0.450 barG, Kd 0.975.
    t_rel = 224.5 + 273.15
    p1 = 4.950e5 + ATM
    p2 = 0.450e5 + ATM
    w = 3000.0 / 3600.0

    # Z from the IAPWS-95 steam tables at relieving conditions; k = Cp0/Cv0 of
    # steam at 497.65 K (ideal-gas heat-capacity ratio).
    pp = make_package("coolprop:Water", ["water"])
    z_comp = p1 * pp.volume(t_rel, p1, {"water": 1.0}) / (8.314462618 * t_rel)
    from CoolProp.CoolProp import PropsSI
    cp0 = PropsSI("Cp0molar", "T", t_rel, "P", p1, "Water")
    k = cp0 / (cp0 - 8.314462618)

    res = relief_vapor(w, t_rel, M_H2O, z_comp, k, p1,
                       backpressure=p2, Kd=0.975)
    print(f"   steam at relieving conditions: Z = {z_comp:.3f}, k = {k:.3f}")
    print(f"   critical-flow pressure = {res.details['P_critical_Pa'] / 1e3:.0f} kPa "
          f"(> backpressure {p2 / 1e3:.0f} kPa -> choked)")
    print(f"   required orifice  = {res.area_m2 * 1e4:6.2f} cm^2   "
          f"[book: 10.17]  ({(res.area_m2 * 1e4 / 10.17 - 1) * 100:+.2f}%)")
    print(f"   selected orifice  = {res.orifice} "
          f"({res.orifice_area_m2 * 1e4:.2f} cm^2)  [book: K (11.86)]")
    print(f"   capacity used     = {res.capacity_used * 100:6.2f} %      "
          f"[book: 85.74]")


if __name__ == "__main__":
    pipe_network()
    psv_sizing()
