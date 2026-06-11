"""M0 demo: build the simplest flowsheet, solve it, print a stream table.

Two feeds -> Mixer -> Heater -> product. The energy balance closes to machine
precision (see tests/test_m0_mixer_heater.py). Run from the repo root with the
engine installed (``pip install -e engine``) or with ``engine/`` on PYTHONPATH:

    python examples/01_mixer_heater.py
"""
import sys
from pathlib import Path

# Allow running straight from the repo without an editable install.
_ENGINE = Path(__file__).resolve().parent.parent / "engine"
if _ENGINE.is_dir() and str(_ENGINE) not in sys.path:
    sys.path.insert(0, str(_ENGINE))

from caldyr.core import Component, Flowsheet  # noqa: E402
from caldyr.unitops import Heater, Mixer  # noqa: E402


def build() -> Flowsheet:
    fs = Flowsheet(
        components=[Component("water"), Component("ethanol")],
        property_package="thermo:PR",
    )
    fs.add(Mixer("MIX1", {"dP": 0.0}))
    fs.add(Heater("H1", {"T_out": 350.0, "dP": 0.0}))

    # Feeds (boundary streams carry a full spec).
    fs.feed("S1", "MIX1:in1", T=298.15, P=101325.0, molar_flow=10.0,
            z={"water": 0.6, "ethanol": 0.4})
    fs.feed("S2", "MIX1:in2", T=320.0, P=101325.0, molar_flow=5.0,
            z={"water": 1.0})
    # Internal + product streams (computed).
    fs.connect("S3", "MIX1:out", "H1:in1")
    fs.connect("S4", "H1:out", None)
    fs.connect("Q1", "H1:duty", None)
    return fs


def stream_table(fs: Flowsheet) -> str:
    head = f"{'stream':>7} {'T/K':>8} {'P/kPa':>8} {'n/(mol/s)':>10} {'phase':>7} {'VF':>6} {'H/(J/mol)':>12}"
    rows = [head, "-" * len(head)]
    for sid in ("S1", "S2", "S3", "S4"):
        s = fs.streams[sid]
        rows.append(
            f"{sid:>7} {s.T:>8.2f} {s.P/1e3:>8.2f} {s.molar_flow:>10.3f} "
            f"{s.phase:>7} {s.vapor_fraction:>6.3f} {s.H:>12.1f}"
        )
    return "\n".join(rows)


def main() -> None:
    fs = build()
    report = fs.solve()
    assert report.converged, report.messages

    print("Solve order:", " -> ".join(report.order))
    print(stream_table(fs))
    print(f"\nHeater H1 duty: {report.duties['Q1'] / 1e3:.2f} kW")

    # Energy balance: Σ(feed n·H) + Q == product n·H
    h_in = sum(fs.streams[s].molar_flow * fs.streams[s].H for s in ("S1", "S2"))
    h_out = fs.streams["S4"].molar_flow * fs.streams["S4"].H
    closure = h_in + report.duties["Q1"] - h_out
    print(f"Energy balance residual: {closure:.3e} W (closes to ~machine precision)")


if __name__ == "__main__":
    main()
