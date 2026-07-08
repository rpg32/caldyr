import { describe, expect, it } from "vitest";
import { groupUnitTypes, UNIT_CATEGORIES } from "./unitCategories";

const u = (type: string) => ({ type });

describe("groupUnitTypes", () => {
  it("places every unit exactly once and preserves taxonomy order", () => {
    const units = [u("HeatExchanger"), u("Mixer"), u("CSTR"), u("Flash")];
    const groups = groupUnitTypes(units);
    const flat = groups.flatMap((g) => g.units.map((x) => x.type));
    expect(flat.sort()).toEqual(["CSTR", "Flash", "HeatExchanger", "Mixer"].sort());
    // labels come out in taxonomy order
    const labels = groups.map((g) => g.label);
    expect(labels).toEqual(["Feeds & mixing", "Reactors", "Flash & phase separation", "Heat transfer"]);
  });

  it("routes unknown/newly-registered types to a trailing Other section", () => {
    const groups = groupUnitTypes([u("Mixer"), u("BrandNewUnit")]);
    const last = groups[groups.length - 1];
    expect(last.label).toBe("Other");
    expect(last.units.map((x) => x.type)).toEqual(["BrandNewUnit"]);
  });

  it("omits empty categories", () => {
    const groups = groupUnitTypes([u("Mixer")]);
    expect(groups).toHaveLength(1);
    expect(groups[0].label).toBe("Feeds & mixing");
  });

  it("has no duplicate type across taxonomy categories", () => {
    const seen = new Set<string>();
    for (const c of UNIT_CATEGORIES)
      for (const t of c.types) {
        expect(seen.has(t)).toBe(false);
        seen.add(t);
      }
  });
});
