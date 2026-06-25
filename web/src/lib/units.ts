// Unit-system registry. The engine (and .flow files) are ALWAYS SI internally
// (K, Pa, mol/s, kg/s, W, W/K); these convert for the I/O edge only — display
// of outputs and entry of inputs. One app-wide UnitSet (SI | Metric | Field)
// picks a default display unit per physical dimension; per-field overrides
// (a future step) can pick any unit listed for the dimension.
//
// Unit choices and conversion factors mirror Aspen HYSYS's unit-set editor so
// flowsheets read the same as in HYSYS. SI/Metric/Field roughly map to HYSYS's
// SI/EuroSI/Field default sets.

export type UnitSet = "SI" | "Metric" | "Field";
export const UNIT_SETS: UnitSet[] = ["SI", "Metric", "Field"];

// Physical dimensions we convert. (Fractions, counts and dimensionless params
// are not converted — they carry no unit.)
export type Dim =
  | "temperature" | "pressure" | "molar_flow" | "mass_flow" | "power" | "UA";

interface UnitDef {
  label: string;
  toSI: (x: number) => number;   // display value -> SI base
  fromSI: (si: number) => number; // SI base -> display value
}

interface DimDef {
  base: string;                       // SI base unit label
  units: UnitDef[];                   // first entry is the SI base
  system: Record<UnitSet, string>;    // default unit label per unit set
}

/** A linear unit: SI = x * factor. */
const lin = (label: string, factor: number): UnitDef => ({
  label,
  toSI: (x) => x * factor,
  fromSI: (si) => si / factor,
});

const LB = 0.45359237;        // kg per lb
const BTU = 1055.056;         // J per IT Btu
const KCAL = 4186.8;          // J per IT kcal

export const DIMENSIONS: Record<Dim, DimDef> = {
  temperature: {
    base: "K",
    units: [
      { label: "K", toSI: (x) => x, fromSI: (si) => si },
      { label: "°C", toSI: (x) => x + 273.15, fromSI: (si) => si - 273.15 },
      { label: "°F", toSI: (x) => (x + 459.67) * 5 / 9, fromSI: (si) => si * 9 / 5 - 459.67 },
      { label: "R", toSI: (x) => x * 5 / 9, fromSI: (si) => si * 9 / 5 },
    ],
    system: { SI: "K", Metric: "°C", Field: "°F" },
  },
  pressure: {
    base: "Pa",
    units: [
      lin("Pa", 1), lin("kPa", 1e3), lin("bar", 1e5), lin("MPa", 1e6),
      lin("psia", 6894.757), lin("atm", 101325), lin("mmHg", 133.322),
      lin("kg/cm²", 98066.5),
    ],
    system: { SI: "kPa", Metric: "bar", Field: "psia" },
  },
  molar_flow: {
    base: "mol/s",
    units: [
      lin("mol/s", 1), lin("mol/h", 1 / 3600), lin("kmol/h", 1000 / 3600),
      lin("kmol/s", 1000), lin("lbmol/h", LB * 1000 / 3600), lin("lbmol/s", LB * 1000),
    ],
    system: { SI: "mol/s", Metric: "kmol/h", Field: "lbmol/h" },
  },
  mass_flow: {
    base: "kg/s",
    units: [
      lin("kg/s", 1), lin("kg/h", 1 / 3600), lin("kg/min", 1 / 60),
      lin("g/s", 1e-3), lin("tonne/h", 1000 / 3600),
      lin("lb/h", LB / 3600), lin("lb/s", LB), lin("lb/min", LB / 60),
    ],
    system: { SI: "kg/s", Metric: "kg/h", Field: "lb/h" },
  },
  power: {
    base: "W",
    units: [
      lin("W", 1), lin("kW", 1e3), lin("MW", 1e6), lin("kJ/h", 1000 / 3600),
      lin("Btu/h", BTU / 3600), lin("MMBtu/h", 1e6 * BTU / 3600),
      lin("kcal/h", KCAL / 3600), lin("hp", 745.6999),
    ],
    system: { SI: "kW", Metric: "kW", Field: "MMBtu/h" },
  },
  UA: {
    base: "W/K",
    units: [
      lin("W/K", 1), lin("kW/K", 1e3), lin("kJ/C-h", 1000 / 3600),
      lin("Btu/F-h", BTU / (5 / 9) / 3600), lin("kcal/C-h", KCAL / 3600),
    ],
    system: { SI: "W/K", Metric: "W/K", Field: "Btu/F-h" },
  },
};

/** All selectable unit labels for a dimension (for a per-field picker). */
export const unitsForDim = (dim: Dim): string[] =>
  DIMENSIONS[dim].units.map((u) => u.label);

/** The default display unit label for a dimension in a unit set. */
export const defaultUnit = (dim: Dim, set: UnitSet): string =>
  DIMENSIONS[dim].system[set];

function unitDef(dim: Dim, label: string): UnitDef {
  const d = DIMENSIONS[dim];
  return d.units.find((u) => u.label === label) ?? d.units[0];
}

/** SI base value -> display value (in `unit`, or the set default). */
export const toDisplay = (dim: Dim, si: number, set: UnitSet, unit?: string): number =>
  unitDef(dim, unit ?? defaultUnit(dim, set)).fromSI(si);

/** Display value -> SI base value (from `unit`, or the set default). */
export const toSI = (dim: Dim, value: number, set: UnitSet, unit?: string): number =>
  unitDef(dim, unit ?? defaultUnit(dim, set)).toSI(value);

/** Format an SI value in the set's display unit (no unit suffix). "—" if null. */
export function fmtDim(
  dim: Dim, si: number | null | undefined, set: UnitSet, digits = 3, unit?: string,
): string {
  if (si == null || !Number.isFinite(si)) return "—";
  return toDisplay(dim, si, set, unit)
    .toLocaleString(undefined, { maximumFractionDigits: digits });
}
