import { describe, expect, it } from "vitest";
import type { SolveResponse } from "../types";
import { evaluateCalcs } from "./calc";

const SOLVE: SolveResponse = {
  report: {
    converged: true, iterations: 1, residual: 0, method: "direct",
    tear_streams: [], duties: { H_duty: 5000 }, messages: [], history: [],
  },
  streams: {
    S1: { id: "S1", T: 350, P: 2e5, molar_flow: 10, z: { a: 0.5, b: 0.5 },
          H: -1000, phase: "liquid", vapor_fraction: 0 },
  },
};

describe("evaluateCalcs", () => {
  it("evaluates accessors, Math, and chained rows", () => {
    const out = evaluateCalcs([
      { name: "tC", expr: 'T("S1") - 273.15' },
      { name: "rate_a", expr: 'n("S1") * z("S1", "a")' },
      { name: "specific", expr: 'duty("H_duty") / rate_a' },
      { name: "rooted", expr: "Math.sqrt(rate_a)" },
    ], SOLVE);
    expect(out[0].value).toBeCloseTo(76.85, 2);
    expect(out[1].value).toBeCloseTo(5, 9);
    expect(out[2].value).toBeCloseTo(1000, 9);
    expect(out[3].value).toBeCloseTo(Math.sqrt(5), 9);
    expect(out.every((r) => r.error === null)).toBe(true);
  });

  it("reports clear errors for unknown streams, duties, and names", () => {
    const out = evaluateCalcs([
      { name: "x", expr: 'T("NOPE")' },
      { name: "y", expr: 'duty("missing")' },
      { name: "zz", expr: "undefined_name + 1" },
    ], SOLVE);
    expect(out[0].error).toMatch(/no solved stream/);
    expect(out[1].error).toMatch(/no duty/);
    expect(out[2].error).toMatch(/unknown name/);
  });

  it("rejects dangerous syntax outright", () => {
    for (const expr of [
      "globalThis.x", "this.constructor", "(()=>1)()", "a; b",
      "Function('return 1')()", "T('S1') = 5", "`tpl`",
    ]) {
      const [r] = evaluateCalcs([{ name: "bad", expr }], SOLVE);
      expect(r.error).not.toBeNull();
    }
  });

  it("supports the full operator grammar with lazy branches", () => {
    const out = evaluateCalcs([
      { name: "pow", expr: "2 ** 3 + 4 % 3" },
      { name: "cmp", expr: '(T("S1") > 300 ? 1 : 0) + (T("S1") < 300 ? 10 : 20)' },
      // the untaken branch must stay unevaluated (duty("missing") would throw)
      { name: "lazy", expr: '1 > 0 ? 2 : duty("missing")' },
      { name: "andor", expr: '(1 > 2 && duty("missing")) || 7' },
    ], SOLVE);
    expect(out[0].value).toBeCloseTo(9, 9);
    expect(out[1].value).toBeCloseTo(21, 9);
    expect(out[2].value).toBeCloseTo(2, 9);
    expect(out[3].value).toBeCloseTo(7, 9);
    expect(out.every((r) => r.error === null)).toBe(true);
  });

  it("restricts member access to Math's own properties", () => {
    for (const expr of ["T.toString(2)", "Math.toString(2)", "Math.hasOwnProperty(1)"]) {
      const [r] = evaluateCalcs([{ name: "bad", expr }], SOLVE);
      expect(r.error).not.toBeNull();
    }
  });

  it("handles unsolved state gracefully", () => {
    const [r] = evaluateCalcs([{ name: "a", expr: 'T("S1")' }], null);
    expect(r.error).toMatch(/no solved stream/);
  });

  it("rejects non-identifier row names and non-finite results", () => {
    const out = evaluateCalcs([
      { name: "bad name", expr: "1" },
      { name: "inf", expr: "1/0" },
    ], SOLVE);
    expect(out[0].error).toMatch(/identifier/);
    expect(out[1].error).toMatch(/finite/);
  });
});
