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
];
