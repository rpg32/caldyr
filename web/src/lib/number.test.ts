import { describe, expect, it } from "vitest";
import { commitNumericDraft, isNumericDraft } from "./number";

describe("isNumericDraft", () => {
  it("accepts partial decimals typed on the way to a value", () => {
    for (const t of ["", "-", "+", ".", ".5", "0.5", "5.", "12", "1e", "1e-", "1.5e3", "05"]) {
      expect(isNumericDraft(t), t).toBe(true);
    }
  });

  it("rejects non-numeric junk and multiple dots", () => {
    for (const t of ["abc", "0.5.5", "1,5", "$3", "5px"]) {
      expect(isNumericDraft(t), t).toBe(false);
    }
  });
});

describe("commitNumericDraft", () => {
  it("commits .5 as 0.5 (the bug: it used to mangle to 05)", () => {
    expect(commitNumericDraft(".5")).toBe(0.5);
  });

  it("commits whole and signed values", () => {
    expect(commitNumericDraft("12")).toBe(12);
    expect(commitNumericDraft("-3.5")).toBe(-3.5);
    expect(commitNumericDraft("1.5e3")).toBe(1500);
  });

  it("returns null for incomplete drafts so the prior value is kept", () => {
    for (const t of ["", "-", ".", "+"]) {
      expect(commitNumericDraft(t), t).toBeNull();
    }
  });
});
