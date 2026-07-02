import { describe, expect, it } from "vitest";
import { humanizeType } from "./format";

describe("humanizeType", () => {
  it("splits CamelCase into spaced words", () => {
    expect(humanizeType("EquilibriumReactor")).toBe("Equilibrium Reactor");
    expect(humanizeType("HeatExchanger")).toBe("Heat Exchanger");
    expect(humanizeType("ThreePhaseSeparator")).toBe("Three Phase Separator");
    expect(humanizeType("RotaryVacuumFilter")).toBe("Rotary Vacuum Filter");
  });

  it("keeps all-caps acronyms intact", () => {
    expect(humanizeType("CSTR")).toBe("CSTR");
    expect(humanizeType("PFR")).toBe("PFR");
  });

  it("splits an acronym run from a trailing word", () => {
    expect(humanizeType("CSTRReactor")).toBe("CSTR Reactor");
  });

  it("leaves plain lowercase words untouched", () => {
    expect(humanizeType("feed")).toBe("feed");
    expect(humanizeType("product")).toBe("product");
  });
});
