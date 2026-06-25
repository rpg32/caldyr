import { describe, expect, it } from "vitest";
import { componentOrder, compositionRows, fmtFrac, streamMassFlow } from "./composition";

describe("compositionRows", () => {
  it("normalizes, sorts descending, and derives per-component flow", () => {
    const rows = compositionRows({ propane: 0.618, "n-butane": 0.382 }, 55.06);
    expect(rows.map((r) => r.comp)).toEqual(["propane", "n-butane"]);
    expect(rows[0].frac).toBeCloseTo(0.618, 6);
    expect(rows[0].flow).toBeCloseTo(55.06 * 0.618, 4);
  });

  it("renormalizes a non-unit sum", () => {
    const rows = compositionRows({ a: 2, b: 2 }, null);
    expect(rows[0].frac).toBeCloseTo(0.5, 9);
    expect(rows[0].flow).toBeNull();
  });

  it("returns empty for missing or zero composition", () => {
    expect(compositionRows(undefined, 1)).toEqual([]);
    expect(compositionRows({}, 1)).toEqual([]);
    expect(compositionRows({ a: 0 }, 1)).toEqual([]);
  });
});

describe("mass basis", () => {
  const MW = { propane: 0.044096, "n-butane": 0.058122 }; // kg/mol

  it("computes mass fraction + mass flow when molar masses are present", () => {
    const rows = compositionRows({ propane: 0.6177, "n-butane": 0.3823 }, 55.06, MW);
    const mixM = 0.6177 * MW.propane + 0.3823 * MW["n-butane"]; // kg/mol
    const c3 = rows.find((r) => r.comp === "propane")!;
    expect(c3.massFrac).toBeCloseTo((0.6177 * MW.propane) / mixM, 6);
    expect(c3.massFlow).toBeCloseTo(55.06 * 0.6177 * MW.propane, 5);
    // mass fractions sum to 1
    expect(rows.reduce((a, r) => a + (r.massFrac ?? 0), 0)).toBeCloseTo(1, 9);
  });

  it("leaves mass fields null when any molar mass is missing", () => {
    const rows = compositionRows({ propane: 0.5, mystery: 0.5 }, 10, MW);
    expect(rows.every((r) => r.massFrac === null && r.massFlow === null)).toBe(true);
  });

  it("streamMassFlow = molar_flow · Σ x·M", () => {
    const m = streamMassFlow({ propane: 0.5, "n-butane": 0.5 }, 100, MW);
    expect(m).toBeCloseTo(100 * (0.5 * MW.propane + 0.5 * MW["n-butane"]), 6);
    expect(streamMassFlow({ propane: 0.5, mystery: 0.5 }, 100, MW)).toBeNull();
    expect(streamMassFlow({ propane: 1 }, null, MW)).toBeNull();
  });
});

describe("componentOrder", () => {
  it("unions component keys across streams, first-seen order", () => {
    const order = componentOrder([
      { z: { propane: 0.6, "n-butane": 0.4 } },
      { z: { "n-butane": 0.7, ethane: 0.3 } },
    ]);
    expect(order).toEqual(["propane", "n-butane", "ethane"]);
  });
});

describe("fmtFrac", () => {
  it("uses fixed notation above 1e-3 and exponential below", () => {
    expect(fmtFrac(0.618)).toBe("0.6180");
    expect(fmtFrac(0)).toBe("0");
    expect(fmtFrac(1.2e-4)).toBe("1.2e-4");
  });
});
