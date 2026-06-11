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
];
