// Display metadata for engine unit-op parameters: label, SI unit, bounds.
// Keys match the names the engine reads from `unit.params` (see engine/caldyr/unitops).
// Unknown keys still render as generic numeric fields.

import type { Dim } from "./units";

// Widget kind for a param. Defaults to "number" when omitted.
export type ParamType = "number" | "boolean" | "select";

export interface ParamMeta {
  label: string;
  unit: string;       // SI unit label; shown as a suffix for dimensionless params
  dim?: Dim;          // physical dimension → unit-system conversion applies
  min?: number;
  max?: number;
  hint?: string;
  type?: ParamType;
  options?: string[];                 // allowed values for type "select"
  // Show this param only when other params match (e.g. condenser_T applies only
  // when decant_condenser is true); the GUI hides it otherwise.
  requires?: Record<string, unknown>;
}

export const PARAM_META: Record<string, ParamMeta> = {
  T: { label: "Temperature", unit: "K", dim: "temperature", min: 1 },
  T_out: { label: "Outlet temperature", unit: "K", dim: "temperature", min: 1 },
  T_hot_out: { label: "Hot outlet T", unit: "K", dim: "temperature", min: 1 },
  T_cold_out: { label: "Cold outlet T", unit: "K", dim: "temperature", min: 1 },
  P: { label: "Pressure", unit: "Pa", dim: "pressure", min: 0 },
  P_out: { label: "Outlet pressure", unit: "Pa", dim: "pressure", min: 0 },
  dP: { label: "Pressure drop", unit: "Pa", dim: "pressure", min: 0 },
  dP_hot: { label: "Hot-side ΔP", unit: "Pa", dim: "pressure", min: 0 },
  dP_cold: { label: "Cold-side ΔP", unit: "Pa", dim: "pressure", min: 0 },
  Q: { label: "Duty", unit: "W", dim: "power" },
  duty: { label: "Duty", unit: "W", dim: "power" },
  UA: { label: "UA", unit: "W/K", dim: "UA", min: 0 },
  eta: { label: "Isentropic efficiency", unit: "–", min: 0, max: 1 },
  split: { label: "Split fraction → out1", unit: "–", min: 0, max: 1 },
  conversion: { label: "Conversion", unit: "–", min: 0, max: 1 },
  molar_flow: { label: "Molar flow", unit: "mol/s", dim: "molar_flow", min: 0 },
  // -- distillation column --------------------------------------------------
  n_stages: { label: "Stages", unit: "–", min: 3 },
  reflux_ratio: { label: "Reflux ratio", unit: "–", min: 0 },
  distillate_rate: { label: "Distillate rate", unit: "mol/s", dim: "molar_flow", min: 0 },
  method: {
    label: "Solver method", unit: "", type: "select",
    options: ["bubble_point", "sum_rates", "inside_out", "naphtali_sandholm"],
  },
  reboiled: { label: "Reboiled", unit: "", type: "boolean" },
  partial_condenser: { label: "Partial condenser", unit: "", type: "boolean" },
  // -- integrated decanting condenser (heteroazeotropic entrainer column) ----
  decant_condenser: {
    label: "Decant condenser", unit: "", type: "boolean",
    hint: "Integrated decanting condenser: the overhead settles into an organic + aqueous layer; the organic layer is refluxed in full and the aqueous layer is the distillate (anhydrous-ethanol entrainer columns).",
  },
  condenser_T: {
    label: "Condenser T", unit: "K", dim: "temperature", min: 150, max: 1500,
    requires: { decant_condenser: true },
    hint: "Temperature the overhead is condensed + decanted at. Required when the decant condenser is on.",
  },
  reflux_layer: {
    label: "Reflux layer", unit: "", type: "select",
    options: ["organic", "aqueous"], requires: { decant_condenser: true },
    hint: "Which settled layer is refluxed in full (organic = the entrainer-rich layer).",
  },
};

export const metaFor = (key: string): ParamMeta =>
  PARAM_META[key] ?? { label: key, unit: "" };

/** Physical dimension of a unit-op param (for unit conversion), or undefined
 *  for dimensionless/unknown params. */
export const dimFor = (key: string): Dim | undefined => PARAM_META[key]?.dim;

/** Dimension of a metric/spec target value (Optimize objective + constraints,
 *  Adjust spec). Covers both the MetricSpec types and the logical-op spec types. */
export function dimForMetric(type: string | undefined): Dim | undefined {
  switch (type) {
    case "T": return "temperature";
    case "P": return "pressure";
    case "flow": case "molar_flow": case "component_rate": return "molar_flow";
    case "duty": return "power";
    default: return undefined;
  }
}

/** Whether a param applies given the current params (its `requires` predicate). */
export function paramApplies(
  key: string,
  params: Record<string, unknown>,
): boolean {
  const req = PARAM_META[key]?.requires;
  if (!req) return true;
  return Object.entries(req).every(([k, v]) => params[k] === v);
}

/** Out-of-bounds check; returns a human message or null when valid. */
export function validateParam(key: string, v: number): string | null {
  if (Number.isNaN(v)) return "not a number";
  const m = PARAM_META[key];
  if (!m) return null;
  if (m.min !== undefined && v < m.min) return `must be ≥ ${m.min}`;
  if (m.max !== undefined && v > m.max) return `must be ≤ ${m.max}`;
  return null;
}

/** Sum of feed mole fractions; valid when within tol of 1. */
export function compositionSum(z: Record<string, number> | undefined): number {
  return Object.values(z ?? {}).reduce((a, b) => a + (Number(b) || 0), 0);
}

export const ID_RE = /^[A-Za-z0-9_-]+$/;
