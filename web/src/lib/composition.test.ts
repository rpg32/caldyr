import { describe, expect, it } from "vitest";
import { componentOrder, compositionRows, fmtFrac } from "./composition";

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
