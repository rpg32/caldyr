import type { FlowDoc } from "./types";

// Prebuilt .flow documents so build -> solve -> cost is one click away.

export const mixerHeater: FlowDoc = {
  schema: "caldyr.flow/1",
  components: [{ id: "water" }, { id: "ethanol" }],
  property_package: "thermo:PR",
  units: [
    { id: "MIX", type: "Mixer", params: { dP: 0.0 }, xy: [200, 160] },
    { id: "H", type: "Heater", params: { T_out: 350.0, dP: 0.0 }, xy: [440, 160] },
  ],
  streams: [
    { id: "S1", from: null, to: "MIX:in1",
      spec: { T: 298.15, P: 101325, molar_flow: 10, z: { water: 0.6, ethanol: 0.4 } } },
    { id: "S2", from: null, to: "MIX:in2",
      spec: { T: 320.0, P: 101325, molar_flow: 5, z: { water: 1.0 } } },
    { id: "S3", from: "MIX:out", to: "H:in1" },
    { id: "S4", from: "H:out", to: null },
    { id: "Q", from: "H:duty", to: null },
  ],
};

export const ammoniaLoop: FlowDoc = {
  schema: "caldyr.flow/1",
  components: [{ id: "nitrogen" }, { id: "hydrogen" }, { id: "ammonia" }, { id: "argon" }],
  property_package: "thermo:PR",
  units: [
    { id: "MIX", type: "Mixer", params: { dP: 0.0 }, xy: [160, 200] },
    { id: "PREHEAT", type: "Heater", params: { T_out: 673.15 }, xy: [380, 200] },
    { id: "RXN", type: "EquilibriumReactor",
      params: { reaction: { stoich: { nitrogen: -1, hydrogen: -3, ammonia: 2 }, key: "nitrogen" }, T: 673.15 },
      xy: [600, 200] },
    { id: "COOL", type: "Heater", params: { T_out: 250.0 }, xy: [820, 200] },
    { id: "SEP", type: "Flash", params: { T: 250.0, P: 2e7 }, xy: [1040, 200] },
    { id: "SPLIT", type: "Splitter", params: { split: 0.9 }, xy: [1040, 380] },
  ],
  streams: [
    { id: "MAKEUP", from: null, to: "MIX:in1",
      spec: { T: 300.0, P: 2e7, molar_flow: 100.0,
              z: { nitrogen: 0.2475, hydrogen: 0.7425, ammonia: 0.0, argon: 0.01 } } },
    { id: "S1", from: "MIX:out", to: "PREHEAT:in1" },
    { id: "S2", from: "PREHEAT:out", to: "RXN:in1" },
    { id: "S3", from: "RXN:out", to: "COOL:in1" },
    { id: "S4", from: "COOL:out", to: "SEP:in1" },
    { id: "PRODUCT", from: "SEP:liquid", to: null },
    { id: "VAP", from: "SEP:vapor", to: "SPLIT:in1" },
    { id: "RECYCLE", from: "SPLIT:out1", to: "MIX:in2" },
    { id: "PURGE", from: "SPLIT:out2", to: null },
  ],
};

export const examples: { name: string; product: string; flow: FlowDoc }[] = [
  { name: "Mixer + Heater", product: "ethanol", flow: mixerHeater },
  { name: "Ammonia loop", product: "ammonia", flow: ammoniaLoop },
];
