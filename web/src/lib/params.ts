// Display metadata for engine unit-op parameters: label, SI unit, bounds.
// Keys match the names the engine reads from `unit.params` (see engine/caldyr/unitops).
// Unknown keys still render as generic numeric fields.

export interface ParamMeta {
  label: string;
  unit: string;       // shown as a suffix; engine is SI throughout
  min?: number;
  max?: number;
  hint?: string;
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
};

export const metaFor = (key: string): ParamMeta =>
  PARAM_META[key] ?? { label: key, unit: "" };

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
