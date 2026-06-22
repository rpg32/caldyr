"""M17 demo: a Claus sulfur-recovery unit (Hameed §10.2.3).

The acid gas an amine sweetening unit strips off (rich in H2S — see
examples/23-25) is fed to a Claus plant, which burns a third of the H2S with air
and reacts the rest to elemental sulfur:

    acid gas (H2S) + air
            |
            v
     [THERMAL FURNACE]  (adiabatic, ~1300 C; sulfur as S2)
            |
            v
       [CONDENSER 1] --> liquid sulfur
            |
        (reheat) -> [CATALYTIC CONVERTER 1] (~530 K; sulfur as S8)
            |
       [CONDENSER 2] --> liquid sulfur
            |
        (reheat) -> [CATALYTIC CONVERTER 2] (~480 K)
            |
       [CONDENSER 3] --> liquid sulfur
            |
            v
         tail gas

The whole train carries the sulfur allotropes (S2 above ~800 K, S8 below) on the
NASA ideal-gas property package (``nasa:claus``) — the cubic EOS cannot, because
``chemicals`` has no critical constants for the S2 dimer. Each reactor takes its
equilibrium from Cantera; each condenser knocks out liquid sulfur down to the
liquid-sulfur saturation pressure.

Key result: the air rate is set to oxidise exactly one-third of the H2S, which
drives the catalytic gas to the stoichiometric H2S:SO2 = 2:1 and *maximises*
sulfur recovery — the Claus air-demand principle. This equilibrium model recovers
~98 % of the feed sulfur (a slight over-prediction of a real, kinetically- and
sub-dewpoint-limited plant) and closes both the sulfur-atom and energy balances
exactly.

    python examples/26_claus_sulfur_recovery.py
"""
import sys
from pathlib import Path

_ENGINE = Path(__file__).resolve().parent.parent / "engine"
if _ENGINE.is_dir() and str(_ENGINE) not in sys.path:
    sys.path.insert(0, str(_ENGINE))

from caldyr.core import Component, Flowsheet  # noqa: E402
from caldyr.unitops import ClausReactor, Heater, SulfurCondenser  # noqa: E402

COMPS = ["hydrogen sulfide", "sulfur dioxide", "S2", "S8", "water", "nitrogen",
         "oxygen", "carbon dioxide", "carbon monoxide", "hydrogen"]

ACID_FLOW = 100.0       # mol/s acid gas
H2S_FRAC = 0.90         # acid gas is 90 % H2S, 10 % CO2
# stoichiometric air: oxidise 1/3 of the H2S to SO2 (H2S + 3/2 O2 -> SO2 + H2O)
O2 = ACID_FLOW * H2S_FRAC / 3.0 * 1.5
N2 = O2 * 79.0 / 21.0


def build() -> Flowsheet:
    fs = Flowsheet(components=[Component(c) for c in COMPS],
                   property_package="nasa:claus")
    fs.add(ClausReactor("FURN", {}))                       # adiabatic furnace
    fs.add(SulfurCondenser("C1", {"T": 450.0}))
    fs.add(Heater("RH1", {"T_out": 530.0}))
    fs.add(ClausReactor("CV1", {"T": 530.0}))              # catalytic bed 1
    fs.add(SulfurCondenser("C2", {"T": 440.0}))
    fs.add(Heater("RH2", {"T_out": 480.0}))
    fs.add(ClausReactor("CV2", {"T": 480.0}))              # catalytic bed 2
    fs.add(SulfurCondenser("C3", {"T": 405.0}))

    fs.feed("ACID", "FURN:in1", T=320.0, P=1.6e5, molar_flow=ACID_FLOW,
            z={"hydrogen sulfide": H2S_FRAC, "carbon dioxide": 1.0 - H2S_FRAC})
    fs.feed("AIR", "FURN:in2", T=480.0, P=1.6e5, molar_flow=O2 + N2,
            z={"oxygen": O2 / (O2 + N2), "nitrogen": N2 / (O2 + N2)})

    fs.connect("g0", "FURN:out", "C1:in1")
    fs.connect("q0", "FURN:duty", None)
    fs.connect("L1", "C1:liquid", None)
    fs.connect("qc1", "C1:duty", None)
    fs.connect("g1", "C1:gas", "RH1:in1")
    fs.connect("qr1", "RH1:duty", None)
    fs.connect("g1b", "RH1:out", "CV1:in1")
    fs.connect("qcv1", "CV1:duty", None)
    fs.connect("g2", "CV1:out", "C2:in1")
    fs.connect("L2", "C2:liquid", None)
    fs.connect("qc2", "C2:duty", None)
    fs.connect("g2b", "C2:gas", "RH2:in1")
    fs.connect("qr2", "RH2:duty", None)
    fs.connect("g2c", "RH2:out", "CV2:in1")
    fs.connect("qcv2", "CV2:duty", None)
    fs.connect("g3", "CV2:out", "C3:in1")
    fs.connect("L3", "C3:liquid", None)
    fs.connect("qc3", "C3:duty", None)
    fs.connect("tail", "C3:gas", None)
    return fs


def main() -> None:
    fs = build()
    rep = fs.solve()
    print(f"converged: {rep.converged}\n")

    furn = fs.streams["g0"]
    print(f"thermal furnace flame temperature: {furn.T - 273.15:6.1f} C")
    print(f"  furnace sulfur is S2 (dimer): y_S2 = {furn.z['S2']:.3f}, "
          f"y_S8 = {furn.z.get('S8', 0.0):.2e}\n")

    s_feed = ACID_FLOW * H2S_FRAC
    print("liquid sulfur drained per condenser (mol/s of S atoms):")
    total = 0.0
    for L in ("L1", "L2", "L3"):
        s = fs.streams[L].molar_flow * 8
        total += s
        print(f"  {L}: {s:7.2f}")
    print(f"  total recovered: {total:7.2f}  of {s_feed:.2f} feed S "
          f"=> recovery {100 * total / s_feed:5.2f} %\n")

    tail = fs.streams["tail"]
    ratio = tail.z["hydrogen sulfide"] / tail.z["sulfur dioxide"]
    print(f"tail-gas H2S:SO2 = {ratio:.2f}  (air-demand optimum is 2.00)")
    print(f"tail-gas residual sulfur vapour (S8): {tail.z.get('S8', 0.0):.2e}\n")

    cooling = sum(d for d in rep.duties.values() if d < 0) / 1e6
    h_in = sum(fs.streams[s].molar_flow * fs.streams[s].H for s in ("ACID", "AIR"))
    h_out = sum(fs.streams[s].molar_flow * fs.streams[s].H
                for s in ("L1", "L2", "L3", "tail"))
    q = sum(rep.duties.values())
    print(f"total cooling duty: {cooling:7.2f} MW")
    print(f"overall energy balance residual: {abs(h_in + q - h_out):.3e} W")


if __name__ == "__main__":
    main()
