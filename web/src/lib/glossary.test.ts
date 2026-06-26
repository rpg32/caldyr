import { describe, expect, it } from "vitest";
import { GLOSSARY, glossaryDef } from "./glossary";

describe("glossaryDef", () => {
  it("looks up terms case-insensitively", () => {
    expect(glossaryDef("LCOP")).toMatch(/Levelized/i);
    expect(glossaryDef("lcop")).toMatch(/Levelized/i);
  });

  it("resolves alternate spellings (aka)", () => {
    expect(glossaryDef("VF")).toBe(glossaryDef("vapor fraction"));
    expect(glossaryDef("dt_min")).toBe(glossaryDef("approach"));
    expect(glossaryDef("CBM")).toBe(glossaryDef("Cbm"));
  });

  it("returns undefined for unknown terms", () => {
    expect(glossaryDef("not-a-term")).toBeUndefined();
  });

  it("every entry has a non-trivial definition", () => {
    for (const e of GLOSSARY) {
      expect(e.def.length, e.term).toBeGreaterThan(20);
    }
  });
});
