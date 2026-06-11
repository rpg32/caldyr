"""Book Ch 15.1 (Hameed 2025): the cyclohexane production plant, end to end.

Benzene hydrogenation (Bz + 3H2 -> C6H12, 98% conversion, adiabatic) with a
feed-effluent heat exchanger, water-free workup (cooler), letdown flash, 70%
vapor recycle with recompression, and an 11-stage column recovering cyclohexane
at the bottom (book: x_cyclohexane ~ 0.998 at ~184 C / 1450 kPa).

Two interlocking recycles: the H2 recycle loop and the feed-effluent HX
thermal loop — a good stress of the tear solver on a real plant topology.
"""
from caldyr.core import Component, Flowsheet
from caldyr.solver import balance_report

KMOLH = 1000.0 / 3600.0   # kmol/h -> mol/s


def build_plant() -> Flowsheet:
    fs = Flowsheet(
        components=[Component(c) for c in
                    ("benzene", "cyclohexane", "hydrogen", "nitrogen", "methane")],
        property_package="thermo:PR",
    )
    from caldyr.unitops import (
        Compressor, FlashDrum, HeatExchanger, Heater, Mixer, Splitter, Valve,
    )
    from caldyr.unitops.conversion_reactor import ConversionReactor
    from caldyr.unitops.rigorous_column import RigorousColumn

    fs.add(Mixer("MIXF", {"dP": 0.0}))      # fresh benzene + hydrogen
    fs.add(Mixer("MIX", {"dP": 0.0}))       # fresh mix + H2 recycle
    fs.add(HeatExchanger("FEHE", {"T_cold_out": 423.15,
                                  "dP_hot": 50e3, "dP_cold": 30e3}))
    fs.add(ConversionReactor("RXN", {
        "reaction": {"stoich": {"benzene": -1, "hydrogen": -3, "cyclohexane": 1},
                     "key": "benzene"},
        "conversion": 0.98, "dP": 100e3,
    }))
    fs.add(Heater("COOL", {"T_out": 323.15, "dP": 20e3}))
    fs.add(Valve("VLV", {"dP": 100e3}))
    fs.add(FlashDrum("SEP", {}))
    fs.add(Splitter("TEE", {"split": 0.7}))
    fs.add(Compressor("K1", {"P_out": 2310e3, "eta": 0.75}))
    # Degassing train: the high-pressure flash liquid still holds ~1% H2 +
    # ~1.8% CH4 — enough that it has NO strict bubble point at column pressure
    # (PR: VF > 0 at every T), which a bubble-point MESH column cannot accept.
    # Letting it down to 300 kPa flashes the lights off; the degassed liquid is
    # pumped back up to the column. (Standard practice; the book's HYSYS column
    # tolerates the lights because its algorithm does not need stage bubble
    # points — flagged as a future column upgrade.)
    fs.add(Valve("VLV2", {"P_out": 300e3}))
    fs.add(FlashDrum("DEGAS", {"T": 333.15, "P": 300e3}))  # heated: drive lights off
    from caldyr.unitops import Pump
    fs.add(Pump("P1", {"P_out": 1450e3, "eta": 0.7}))
    fs.add(RigorousColumn("T100", {
        "n_stages": 11, "feed_stage": 6, "reflux_ratio": 0.9,
        "distillate_rate": 3.4 * KMOLH, "P": 1380e3, "partial_condenser": True,
    }))

    fs.feed("BZ_FEED", "MIXF:in1", T=311.15, P=2310e3, molar_flow=100 * KMOLH,
            z={"benzene": 1.0})
    fs.feed("H2_FEED", "MIXF:in2", T=311.15, P=2310e3, molar_flow=310 * KMOLH,
            z={"hydrogen": 0.98, "nitrogen": 0.005, "methane": 0.015})

    fs.connect("F0", "MIXF:out", "MIX:in1")
    fs.connect("F1", "MIX:out", "FEHE:cold_in")
    fs.connect("TO_REACT", "FEHE:cold_out", "RXN:in1")
    fs.connect("VAP1", "RXN:out", "FEHE:hot_in")
    fs.connect("VAP", "FEHE:hot_out", "COOL:in1")
    fs.connect("TO_VLV", "COOL:out", "VLV:in1")
    fs.connect("TO_SEP", "VLV:out", "SEP:in1")
    fs.connect("REC", "SEP:vapor", "TEE:in1")
    fs.connect("S1", "TEE:out1", "K1:in1")
    fs.connect("H2_REC", "K1:out", "MIX:in2")
    fs.connect("PURGE", "TEE:out2", None)
    fs.connect("LIQ_HP", "SEP:liquid", "VLV2:in1")
    fs.connect("LIQ_LP", "VLV2:out", "DEGAS:in1")
    fs.connect("OFFGAS", "DEGAS:vapor", None)
    fs.connect("LIQ_DEG", "DEGAS:liquid", "P1:in1")
    fs.connect("TO_DIST", "P1:out", "T100:in1")
    fs.connect("OVHD", "T100:distillate", None)
    fs.connect("CYCLOHEXANE", "T100:bottoms", None)

    # The FEHE couples the reactor effluent back onto the feed train, so one of
    # the thermal-loop streams gets torn — and a HeatExchanger (unlike a Mixer)
    # cannot start from an empty tear. Seed whichever stream the solver tears
    # with a rough engineering estimate (the HYSYS-recycle-block equivalent).
    fs.solver_hints = {"tear_guesses": {
        "VAP1": {"T": 520.0, "P": 2.21e6, "molar_flow": 330 * KMOLH,
                 "z": {"cyclohexane": 0.21, "hydrogen": 0.55, "benzene": 0.005,
                       "nitrogen": 0.035, "methane": 0.20}},
        "TO_REACT": {"T": 423.15, "P": 2.28e6, "molar_flow": 530 * KMOLH,
                     "z": {"benzene": 0.18, "hydrogen": 0.70, "nitrogen": 0.02,
                           "methane": 0.10}},
        "H2_REC": {"T": 350.0, "P": 2.31e6, "molar_flow": 120 * KMOLH,
                   "z": {"hydrogen": 0.82, "methane": 0.14, "nitrogen": 0.04}},
    }}
    return fs


if __name__ == "__main__":
    fs = build_plant()
    report = fs.solve(tol=1e-7)
    print(f"converged: {report.converged}  ({report.method}, "
          f"{report.iterations} iters; tears: {report.tear_streams})")
    prod = fs.streams["CYCLOHEXANE"]
    z = prod.normalized_z()
    print(f"cyclohexane product: {prod.molar_flow * 3.6:.1f} kmol/h at "
          f"{prod.T - 273.15:.1f} C, x_C6H12 = {z['cyclohexane']:.4f} "
          f"(book: ~0.998 at ~184 C)")
    ovhd = fs.streams["OVHD"]
    print(f"column overhead: {ovhd.molar_flow * 3.6:.2f} kmol/h (spec 3.4)")
    bal = balance_report(fs)
    print(f"plant-wide mass closure: {bal['overall']['mass_rel']:.2e}")
