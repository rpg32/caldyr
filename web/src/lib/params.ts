// Display metadata for engine unit-op parameters: label, SI unit, bounds.
// Keys match the names the engine reads from `unit.params` (see engine/caldyr/unitops).
// Unknown keys still render as generic numeric fields.

// Widget kind for a param. Defaults to "number" when omitted.
export type ParamType = "number" | "boolean" | "select";

export interface ParamMeta {
  label: string;
  unit: string;       // shown as a suffix; engine is SI throughout
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
  T: { label: "Temperature", unit: "K", min: 1 },
  T_out: { label: "Outlet temperature", unit: "K", min: 1 },
  T_hot_out: { label: "Hot outlet T", unit: "K", min: 1 },
  T_cold_out: { label: "Cold outlet T", unit: "K", min: 1 },
  P: { label: "Pressure", unit: "Pa", min: 0 },
  P_out: { label: "Outlet pressure", unit: "Pa", min: 0 },
  dP: { label: "Pressure drop", unit: "Pa", min: 0 },
  dP_hot: { label: "Hot-side ΔP", unit: "Pa", min: 0 },
  dP_cold: { label: "Cold-side ΔP", unit: "Pa", min: 0 },
  Q: { label: "Duty", unit: "W" },
  duty: { label: "Duty", unit: "W" },
  UA: { label: "UA", unit: "W/K", min: 0 },
  eta: { label: "Isentropic efficiency", unit: "–", min: 0, max: 1 },
  split: { label: "Split fraction → out1", unit: "–", min: 0, max: 1 },
  conversion: { label: "Conversion", unit: "–", min: 0, max: 1 },
  molar_flow: { label: "Molar flow", unit: "mol/s", min: 0 },
  // -- distillation column --------------------------------------------------
  n_stages: { label: "Stages", unit: "–", min: 3 },
  reflux_ratio: { label: "Reflux ratio", unit: "–", min: 0 },
  distillate_rate: { label: "Distillate rate", unit: "mol/s", min: 0 },
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
    label: "Condenser T", unit: "K", min: 150, max: 1500,
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
