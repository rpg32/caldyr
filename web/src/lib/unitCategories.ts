// Display taxonomy for the unit-ops palette. The engine registry stays the
// single source of truth for WHICH units exist; this only decides how they are
// grouped for the palette. Any registered unit not listed here falls into
// "Other", so newly-added units still appear without touching this file.

export interface UnitCategory {
  label: string;
  types: string[];
}

export const UNIT_CATEGORIES: UnitCategory[] = [
  { label: "Feeds & mixing", types: ["Source", "Makeup", "Mixer", "Splitter"] },
  { label: "Reactors", types: ["CSTR", "PFR", "ConversionReactor", "EquilibriumReactor", "GibbsReactor", "ClausReactor"] },
  { label: "Columns", types: ["ShortcutColumn", "RigorousColumn", "Absorber", "ReboiledAbsorber", "ExtractionColumn"] },
  { label: "Flash & phase separation", types: ["Flash", "ThreePhaseSeparator", "Decanter", "ComponentSplitter", "Evaporator", "Saturator"] },
  { label: "Heat transfer", types: ["Heater", "HeatExchanger", "AirCooler", "FiredHeater", "MultiStreamExchanger", "SulfurCondenser"] },
  { label: "Pressure change", types: ["Pump", "Compressor", "Expander", "Valve", "PipeSegment"] },
  { label: "Solids handling", types: ["Cyclone", "BaghouseFilter", "RotaryVacuumFilter"] },
  { label: "Logical & utility", types: ["Balance"] },
];

const OTHER = "Other";

/** Group the registry's unit types into ordered display categories. Unknown
 * types land in a trailing "Other" section; within each section, types keep the
 * taxonomy's order (known) or the registry's order (Other). */
export function groupUnitTypes<T extends { type: string }>(
  units: T[],
): { label: string; units: T[] }[] {
  const byType = new Map(units.map((u) => [u.type, u]));
  const claimed = new Set<string>();
  const groups: { label: string; units: T[] }[] = [];

  for (const cat of UNIT_CATEGORIES) {
    const found = cat.types
      .map((t) => byType.get(t))
      .filter((u): u is T => u !== undefined);
    found.forEach((u) => claimed.add(u.type));
    if (found.length) groups.push({ label: cat.label, units: found });
  }

  const leftover = units.filter((u) => !claimed.has(u.type));
  if (leftover.length) groups.push({ label: OTHER, units: leftover });

  return groups;
}
