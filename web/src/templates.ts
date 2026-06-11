// Template gallery: ready-to-solve flowsheets for the Projects dialog.
import { ammoniaLoop, mixerHeater } from "./examples";
import type { FlowDoc } from "./types";

export const benzeneColumn: FlowDoc = {
  schema: "caldyr.flow/1",
  meta: { ui: { product: "benzene" } },
  components: [{ id: "benzene" }, { id: "toluene" }],
  property_package: "thermo:PR",
  units: [
    { id: "COL", type: "ShortcutColumn",
      params: { light_key: "benzene", heavy_key: "toluene", recovery_light: 0.95,
                recovery_heavy: 0.95, rr_factor: 1.3, P: 101325 },
      xy: [420, 200] },
  ],
  streams: [
    { id: "FEED", from: null, to: "COL:in1",
      spec: { T: 365, P: 101325, molar_flow: 100, z: { benzene: 0.5, toluene: 0.5 } } },
    { id: "DISTILLATE", from: "COL:distillate", to: null },
    { id: "BOTTOMS", from: "COL:bottoms", to: null },
  ],
};

export const smrHydrogen: FlowDoc = {
  schema: "caldyr.flow/1",
  meta: { ui: { product: "hydrogen" } },
  components: [
    { id: "methane" }, { id: "water" }, { id: "hydrogen" },
    { id: "carbon monoxide" }, { id: "carbon dioxide" },
  ],
  property_package: "thermo:PR",
  units: [
    { id: "SMR", type: "GibbsReactor", params: { T: 1100.0, P: 101325.0 }, xy: [320, 180] },
    { id: "COOL", type: "Heater", params: { T_out: 473.15 }, xy: [560, 180] },
    { id: "SEP", type: "Flash", params: { T: 313.15, P: 101325.0 }, xy: [800, 180] },
  ],
  streams: [
    { id: "FEED", from: null, to: "SMR:in1",
      spec: { T: 800.0, P: 101325.0, molar_flow: 40.0,
              z: { methane: 0.25, water: 0.75 } } },
    { id: "S1", from: "SMR:out", to: "COOL:in1" },
    { id: "S2", from: "COOL:out", to: "SEP:in1" },
    { id: "SYNGAS", from: "SEP:vapor", to: null },
    { id: "CONDENSATE", from: "SEP:liquid", to: null },
  ],
};

// Book ch. 15.1 (Hameed 2025): benzene hydrogenation with a feed-effluent HX,
// H2 recycle with recompression, degassing train, and an 11-stage column.
// The solver_hints seed the FEHE/recycle tears (HXs can't start from empty).
export const cyclohexanePlant: FlowDoc = {
  schema: "caldyr.flow/1",
  meta: { ui: { product: "cyclohexane" } },
  components: [
    { id: "benzene" }, { id: "cyclohexane" }, { id: "hydrogen" },
    { id: "nitrogen" }, { id: "methane" },
  ],
  property_package: "thermo:PR",
  units: [
    { id: "MIXF", type: "Mixer", params: { dP: 0 }, xy: [60, 200] },
    { id: "MIX", type: "Mixer", params: { dP: 0 }, xy: [240, 200] },
    { id: "FEHE", type: "HeatExchanger",
      params: { T_cold_out: 423.15, dP_hot: 50e3, dP_cold: 30e3 }, xy: [430, 200] },
    { id: "RXN", type: "ConversionReactor",
      params: { reaction: { stoich: { benzene: -1, hydrogen: -3, cyclohexane: 1 },
                            key: "benzene" },
                conversion: 0.98, dP: 100e3 }, xy: [640, 120] },
    { id: "COOL", type: "Heater", params: { T_out: 323.15, dP: 20e3 }, xy: [640, 300] },
    { id: "VLV", type: "Valve", params: { dP: 100e3 }, xy: [830, 300] },
    { id: "SEP", type: "Flash", params: {}, xy: [1010, 300] },
    { id: "TEE", type: "Splitter", params: { split: 0.7 }, xy: [1010, 110] },
    { id: "K1", type: "Compressor", params: { P_out: 2310e3, eta: 0.75 }, xy: [620, 20] },
    { id: "VLV2", type: "Valve", params: { P_out: 300e3 }, xy: [1010, 460] },
    { id: "DEGAS", type: "Flash", params: { T: 333.15, P: 300e3 }, xy: [1200, 460] },
    { id: "P1", type: "Pump", params: { P_out: 1450e3, eta: 0.7 }, xy: [1390, 460] },
    { id: "T100", type: "RigorousColumn",
      params: { n_stages: 11, feed_stage: 6, reflux_ratio: 0.9,
                distillate_rate: 0.9444, P: 1380e3, partial_condenser: true },
      xy: [1580, 460] },
  ],
  streams: [
    { id: "BZ_FEED", from: null, to: "MIXF:in1",
      spec: { T: 311.15, P: 2310e3, molar_flow: 27.778, z: { benzene: 1.0 } } },
    { id: "H2_FEED", from: null, to: "MIXF:in2",
      spec: { T: 311.15, P: 2310e3, molar_flow: 86.111,
              z: { hydrogen: 0.98, nitrogen: 0.005, methane: 0.015 } } },
    { id: "F0", from: "MIXF:out", to: "MIX:in1" },
    { id: "F1", from: "MIX:out", to: "FEHE:cold_in" },
    { id: "TO_REACT", from: "FEHE:cold_out", to: "RXN:in1" },
    { id: "VAP1", from: "RXN:out", to: "FEHE:hot_in" },
    { id: "VAP", from: "FEHE:hot_out", to: "COOL:in1" },
    { id: "TO_VLV", from: "COOL:out", to: "VLV:in1" },
    { id: "TO_SEP", from: "VLV:out", to: "SEP:in1" },
    { id: "REC", from: "SEP:vapor", to: "TEE:in1" },
    { id: "S1", from: "TEE:out1", to: "K1:in1" },
    { id: "H2_REC", from: "K1:out", to: "MIX:in2" },
    { id: "PURGE", from: "TEE:out2", to: null },
    { id: "LIQ_HP", from: "SEP:liquid", to: "VLV2:in1" },
    { id: "LIQ_LP", from: "VLV2:out", to: "DEGAS:in1" },
    { id: "OFFGAS", from: "DEGAS:vapor", to: null },
    { id: "LIQ_DEG", from: "DEGAS:liquid", to: "P1:in1" },
    { id: "TO_DIST", from: "P1:out", to: "T100:in1" },
    { id: "OVHD", from: "T100:distillate", to: null },
    { id: "CYCLOHEXANE", from: "T100:bottoms", to: null },
  ],
  solver_hints: {
    tear_guesses: {
      VAP1: { T: 520, P: 2.21e6, molar_flow: 91.7,
              z: { cyclohexane: 0.21, hydrogen: 0.55, benzene: 0.005,
                   nitrogen: 0.035, methane: 0.2 } },
      TO_REACT: { T: 423.15, P: 2.28e6, molar_flow: 147.2,
                  z: { benzene: 0.18, hydrogen: 0.7, nitrogen: 0.02, methane: 0.1 } },
      H2_REC: { T: 350, P: 2.31e6, molar_flow: 33.3,
                z: { hydrogen: 0.82, methane: 0.14, nitrogen: 0.04 } },
    },
  },
};

// Book ch. 15.6 (Hameed 2025): vinyl chloride monomer by EDC cracking,
// with the EDC recycle closed. 500 kmol/h EDC + recycle -> furnace-style
// heater -> isothermal cracker (55% per pass) -> air cooler -> HCl tower
// (25 bar) -> VCM tower (15 bar) -> EDC bottoms pumped back to the mixer.
// Deviations from the book, both required by the closed recycle: (1) the
// column draws are ratio specs (distillate_to_feed) instead of fixed rates,
// so the recycle self-stabilizes (the book's fixed 275 / 272 kmol/h draws
// are per-pass numbers — closing the loop roughly doubles production to
// ~500 kmol/h of each product); (2) the HCl tower runs RR 3 (book: 1.5) —
// the bubble-point MESH needs the extra reflux to keep near-critical HCl
// (Tc = 324.7 K) out of the hot stage liquids.
export const vcmPlant: FlowDoc = {
  schema: "caldyr.flow/1",
  meta: { ui: { product: "vinyl chloride" } },
  components: [
    { id: "1,2-dichloroethane" }, { id: "hydrogen chloride" },
    { id: "vinyl chloride" },
  ],
  property_package: "thermo:SRK",
  units: [
    { id: "MIX", type: "Mixer", params: { dP: 0 }, xy: [60, 200] },
    { id: "HEAT", type: "Heater",
      params: { T_out: 443.15, dP: 50e3 }, xy: [250, 200] },
    { id: "RXN", type: "ConversionReactor",
      params: { reaction: { stoich: { "1,2-dichloroethane": -1,
                                      "hydrogen chloride": 1,
                                      "vinyl chloride": 1 },
                            key: "1,2-dichloroethane" },
                conversion: 0.55, T_out: 443.15, dP: 50e3 },
      xy: [440, 200] },
    { id: "ACOOL", type: "AirCooler",
      params: { T_out: 343.15, dP: 50e3, t_air_in: 293.15 }, xy: [640, 200] },
    { id: "T100", type: "RigorousColumn",
      params: { n_stages: 8, feed_stage: 4, reflux_ratio: 3.0,
                distillate_to_feed: 0.3475, P: 25e5 },
      xy: [840, 200] },
    { id: "T101", type: "RigorousColumn",
      params: { n_stages: 8, feed_stage: 4, reflux_ratio: 1.5,
                distillate_to_feed: 0.545, P: 15e5 },
      xy: [1040, 320] },
    { id: "P1", type: "Pump", params: { P_out: 30e5, eta: 0.75 }, xy: [540, 470] },
  ],
  streams: [
    { id: "EDC_FEED", from: null, to: "MIX:in1",
      spec: { T: 293.15, P: 30e5, molar_flow: 138.89,
              z: { "1,2-dichloroethane": 1.0 } } },
    { id: "S1", from: "MIX:out", to: "HEAT:in1" },
    { id: "S2", from: "HEAT:out", to: "RXN:in1" },
    { id: "S3", from: "RXN:out", to: "ACOOL:in1" },
    { id: "S4", from: "ACOOL:out", to: "T100:in1" },
    { id: "HCL", from: "T100:distillate", to: null },
    { id: "S5", from: "T100:bottoms", to: "T101:in1" },
    { id: "VCM", from: "T101:distillate", to: null },
    { id: "S6", from: "T101:bottoms", to: "P1:in1" },
    { id: "EDC_REC", from: "P1:out", to: "MIX:in2" },
  ],
  solver_hints: {
    tear_guesses: {
      EDC_REC: { T: 460, P: 30e5, molar_flow: 113.6,
                 z: { "1,2-dichloroethane": 1.0 } },
    },
  },
};

// Book ch. 15.2 (Hameed 2025), simplified: DME by methanol dehydration.
// Fresh methanol (260 kmol/h) is pumped to 13 atm, mixed with the methanol
// recycle, vaporized/heated to 250 C, and converted adiabatically
// (2 MeOH -> DME + H2O, 30% per pass, outlet ~295 C). The effluent is cooled
// to 100 C, let down to 10 atm for the DME tower (DME overhead, ~47 C — a
// cooling-water condenser), and the bottoms let down to 7 atm for the
// methanol tower (MeOH overhead, recycled; water bottoms >99.9%). Property
// package is thermo:PR — ChemSep carries no NRTL parameters for the
// DME/methanol and DME/water pairs (NRTL would silently treat those pairs as
// ideal), and the MeOH/water/DME system has no azeotrope for PR to miss.
// Column draws are ratio specs so the recycle self-stabilizes.
export const dmePlant: FlowDoc = {
  schema: "caldyr.flow/1",
  meta: { ui: { product: "dimethyl ether" } },
  components: [{ id: "methanol" }, { id: "dimethyl ether" }, { id: "water" }],
  property_package: "thermo:PR",
  units: [
    { id: "P1", type: "Pump", params: { P_out: 1317225, eta: 0.75 }, xy: [60, 200] },
    { id: "MIX", type: "Mixer", params: { dP: 0 }, xy: [240, 200] },
    { id: "HEAT", type: "Heater",
      params: { T_out: 523.15, dP: 30e3 }, xy: [430, 200] },
    { id: "RXN", type: "ConversionReactor",
      params: { reaction: { stoich: { methanol: -2, "dimethyl ether": 1, water: 1 },
                            key: "methanol" },
                conversion: 0.30, dP: 50e3 },
      xy: [620, 200] },
    { id: "COOL", type: "Heater",
      params: { T_out: 373.15, dP: 20e3 }, xy: [810, 200] },
    { id: "VLV", type: "Valve", params: { P_out: 1013250 }, xy: [990, 200] },
    { id: "T100", type: "RigorousColumn",
      params: { n_stages: 14, feed_stage: 7, reflux_ratio: 2.5,
                distillate_to_feed: 0.145, P: 1013250, max_iter: 600 },
      xy: [1180, 200] },
    { id: "VLV2", type: "Valve", params: { P_out: 709275 }, xy: [1180, 420] },
    { id: "T101", type: "RigorousColumn",
      params: { n_stages: 12, feed_stage: 6, reflux_ratio: 2.0,
                distillate_to_feed: 0.838, P: 709275, max_iter: 600 },
      xy: [1370, 420] },
    { id: "P2", type: "Pump", params: { P_out: 1317225, eta: 0.75 }, xy: [620, 420] },
  ],
  streams: [
    { id: "MEOH_FEED", from: null, to: "P1:in1",
      spec: { T: 298.15, P: 101325, molar_flow: 72.22, z: { methanol: 1.0 } } },
    { id: "S0", from: "P1:out", to: "MIX:in1" },
    { id: "S1", from: "MIX:out", to: "HEAT:in1" },
    { id: "S2", from: "HEAT:out", to: "RXN:in1" },
    { id: "S3", from: "RXN:out", to: "COOL:in1" },
    { id: "S4", from: "COOL:out", to: "VLV:in1" },
    { id: "S5", from: "VLV:out", to: "T100:in1" },
    { id: "DME", from: "T100:distillate", to: null },
    { id: "S6", from: "T100:bottoms", to: "VLV2:in1" },
    { id: "S7", from: "VLV2:out", to: "T101:in1" },
    { id: "MEOH_REC", from: "T101:distillate", to: "P2:in1" },
    { id: "WATER", from: "T101:bottoms", to: null },
    { id: "REC", from: "P2:out", to: "MIX:in2" },
  ],
  solver_hints: {
    tear_guesses: {
      REC: { T: 410, P: 1317225, molar_flow: 168.5,
             z: { methanol: 0.98, water: 0.02 } },
    },
  },
};

export interface Template {
  name: string;
  blurb: string;
  product: string;
  flow: FlowDoc;
}

export const TEMPLATES: Template[] = [
  { name: "Mixer + Heater", blurb: "The simplest flowsheet — two feeds, one heater.",
    product: "ethanol", flow: mixerHeater },
  { name: "Ammonia loop", blurb: "Haber-Bosch synthesis loop with recycle and purge.",
    product: "ammonia", flow: ammoniaLoop },
  { name: "Benzene/toluene column", blurb: "Shortcut (FUG) distillation of an equimolar feed.",
    product: "benzene", flow: benzeneColumn },
  { name: "SMR hydrogen", blurb: "Steam-methane reforming via Gibbs equilibrium, cooled and flashed.",
    product: "hydrogen", flow: smrHydrogen },
  { name: "Cyclohexane plant", blurb: "Book ch. 15.1: benzene hydrogenation with feed-effluent HX, H2 recycle, degasser, and column (two interlocked recycles — takes a minute).",
    product: "cyclohexane", flow: cyclohexanePlant },
  { name: "VCM plant", blurb: "Book ch. 15.6: EDC cracked to vinyl chloride + HCl, separated in two pressure columns with the unconverted-EDC recycle closed (ratio draws and RR 3 on the HCl tower keep the closed loop solvable — see the doc comment).",
    product: "vinyl chloride", flow: vcmPlant },
  { name: "DME plant", blurb: "Book ch. 15.2 (simplified): methanol dehydrated to dimethyl ether over an adiabatic reactor, with a two-column separation and the methanol recycle closed (PR — ChemSep has no NRTL parameters for the DME pairs).",
    product: "dimethyl ether", flow: dmePlant },
];
