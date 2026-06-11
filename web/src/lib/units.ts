// Display unit sets. The engine (and .flow files) are always SI (K, Pa, mol/s);
// these convert for DISPLAY ONLY at the last moment before rendering.

export type UnitSet = "SI" | "Metric" | "Field";
export const UNIT_SETS: UnitSet[] = ["SI", "Metric", "Field"];

interface Quantity {
  unit: Record<UnitSet, string>;
  to: Record<UnitSet, (si: number) => number>;
}

const QUANTITIES: Record<"T" | "P" | "flow" | "power", Quantity> = {
  T: {
    unit: { SI: "K", Metric: "°C", Field: "°F" },
    to: {
      SI: (k) => k,
      Metric: (k) => k - 273.15,
      Field: (k) => (k - 273.15) * 9 / 5 + 32,
    },
  },
  P: {
    unit: { SI: "kPa", Metric: "bar", Field: "psia" },
    to: {
      SI: (pa) => pa / 1e3,
      Metric: (pa) => pa / 1e5,
      Field: (pa) => pa / 6894.757,
    },
  },
  flow: {
    unit: { SI: "mol/s", Metric: "kmol/h", Field: "lbmol/h" },
    to: {
      SI: (mols) => mols,
      Metric: (mols) => mols * 3.6,
      Field: (mols) => mols * 3600 / 453.59237,
    },
  },
  power: {
    unit: { SI: "kW", Metric: "kW", Field: "MMBtu/h" },
    to: {
      SI: (w) => w / 1e3,
      Metric: (w) => w / 1e3,
      Field: (w) => w * 3.412142e-6, // W -> MMBtu/h
    },
  },
};

export const convert = (q: keyof typeof QUANTITIES, si: number, set: UnitSet): number =>
  QUANTITIES[q].to[set](si);

export const unitOf = (q: keyof typeof QUANTITIES, set: UnitSet): string =>
  QUANTITIES[q].unit[set];

export const fmtQty = (
  q: keyof typeof QUANTITIES, si: number | null | undefined, set: UnitSet, digits = 2,
): string =>
  si == null ? "—"
    : convert(q, si, set).toLocaleString(undefined, { maximumFractionDigits: digits });
