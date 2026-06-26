import { describe, expect, it } from "vitest";
import type { ReactionEditorOpts } from "../types";
import {
  emptyDraft,
  normalizeReactions,
  previewReaction,
  serializeReactions,
  summarizeReactions,
  validateReactions,
} from "./reactions";

const EQUIL: ReactionEditorOpts = { kind: "stoichiometric", multiple: false, conversion: false, key_required: false };
const CONV: ReactionEditorOpts = { kind: "stoichiometric", multiple: true, conversion: true, key_required: true };
const KIN: ReactionEditorOpts = { kind: "kinetic", multiple: true, conversion: false, key_required: true, reversible: true };

const COMPS = ["nitrogen", "hydrogen", "ammonia", "argon"];

describe("normalizeReactions", () => {
  it("splits a stoich map into positive reactant/product rows", () => {
    const drafts = normalizeReactions({
      reaction: { stoich: { nitrogen: -1, hydrogen: -3, ammonia: 2 }, key: "nitrogen" },
    });
    expect(drafts).toHaveLength(1);
    const d = drafts[0];
    expect(d.reactants).toEqual([
      { comp: "nitrogen", coeff: 1 },
      { comp: "hydrogen", coeff: 3 },
    ]);
    expect(d.products).toEqual([{ comp: "ammonia", coeff: 2 }]);
    expect(d.key).toBe("nitrogen");
  });

  it("pulls unit-level conversion onto a singular reaction", () => {
    const drafts = normalizeReactions({
      reaction: { stoich: { benzene: -1, hydrogen: -3, cyclohexane: 1 }, key: "benzene" },
      conversion: 0.98,
    });
    expect(drafts[0].conversion).toBe(0.98);
  });

  it("reads the list form and per-reaction conversion", () => {
    const drafts = normalizeReactions({
      reactions: [
        { stoich: { a: -1, b: 1 }, key: "a", conversion: 0.5 },
        { stoich: { b: -1, c: 1 }, key: "b", conversion: 0.3 },
      ],
    });
    expect(drafts).toHaveLength(2);
    expect(drafts[1].conversion).toBe(0.3);
  });

  it("reads kinetic fields including reversible", () => {
    const drafts = normalizeReactions({
      reactions: [
        { stoich: { a: -1, b: 1 }, key: "a", k0: 5, Ea: 1000, orders: { a: 2 }, k0_rev: 0.1, Ea_rev: 500 },
      ],
    });
    const d = drafts[0];
    expect(d.k0).toBe(5);
    expect(d.Ea).toBe(1000);
    expect(d.orders).toEqual([{ comp: "a", coeff: 2 }]);
    expect(d.reversible).toBe(true);
    expect(d.k0_rev).toBe(0.1);
  });

  it("returns one empty draft when nothing is set", () => {
    const drafts = normalizeReactions({});
    expect(drafts).toHaveLength(1);
    expect(drafts[0].reactants[0].comp).toBe("");
  });
});

describe("serializeReactions", () => {
  it("round-trips the ammonia equilibrium reaction to signed stoich", () => {
    const drafts = normalizeReactions({
      reaction: { stoich: { nitrogen: -1, hydrogen: -3, ammonia: 2 }, key: "nitrogen" },
    });
    const out = serializeReactions(drafts, EQUIL);
    expect(out).toEqual([{ stoich: { nitrogen: -1, hydrogen: -3, ammonia: 2 }, key: "nitrogen" }]);
  });

  it("writes per-reaction conversion for ConversionReactor", () => {
    const drafts = normalizeReactions({ reaction: { stoich: { a: -1, b: 1 }, key: "a" }, conversion: 0.8 });
    const out = serializeReactions(drafts, CONV);
    expect(out[0]).toEqual({ stoich: { a: -1, b: 1 }, key: "a", conversion: 0.8 });
  });

  it("writes kinetic + reversible fields", () => {
    const d = emptyDraft();
    d.reactants = [{ comp: "a", coeff: 1 }];
    d.products = [{ comp: "b", coeff: 1 }];
    d.key = "a";
    d.k0 = 5;
    d.Ea = 1000;
    d.orders = [{ comp: "a", coeff: 2 }];
    d.reversible = true;
    d.k0_rev = 0.1;
    d.Ea_rev = 500;
    const out = serializeReactions([d], KIN);
    expect(out[0]).toEqual({
      stoich: { a: -1, b: 1 },
      key: "a",
      k0: 5,
      Ea: 1000,
      orders: { a: 2 },
      k0_rev: 0.1,
      Ea_rev: 500,
    });
  });

  it("is a fixed point: normalize ∘ serialize preserves a kinetic reaction", () => {
    const original = { reactions: [{ stoich: { a: -2, b: 1, c: 1 }, key: "a", k0: 3.2, Ea: 5e4, orders: { a: 2 } }] };
    const drafts = normalizeReactions(original);
    const out = serializeReactions(drafts, KIN);
    expect(out).toEqual(original.reactions);
  });

  it("omits reverse fields when not reversible", () => {
    const d = emptyDraft();
    d.reactants = [{ comp: "a", coeff: 1 }];
    d.products = [{ comp: "b", coeff: 1 }];
    d.key = "a";
    const out = serializeReactions([d], KIN) as unknown as Record<string, unknown>[];
    expect(out[0]).not.toHaveProperty("k0_rev");
    expect(out[0]).not.toHaveProperty("orders");
  });
});

describe("previewReaction / summarizeReactions", () => {
  it("omits unit coefficients and shows the arrow", () => {
    const d = normalizeReactions({
      reaction: { stoich: { nitrogen: -1, hydrogen: -3, ammonia: 2 } },
    })[0];
    expect(previewReaction(d)).toBe("nitrogen + 3 hydrogen → 2 ammonia");
  });

  it("summarizes none vs one vs many", () => {
    expect(summarizeReactions({})).toBe("— none —");
    expect(summarizeReactions({ reaction: { stoich: { a: -1, b: 1 } } })).toBe("1 reaction: a → b");
    expect(
      summarizeReactions({ reactions: [{ stoich: { a: -1, b: 1 } }, { stoich: { b: -1, c: 1 } }] }),
    ).toMatch(/^2 reactions: a → b/);
  });
});

describe("validateReactions", () => {
  const ammonia = () =>
    normalizeReactions({ reaction: { stoich: { nitrogen: -1, hydrogen: -3, ammonia: 2 }, key: "nitrogen" } });

  it("passes a well-formed equilibrium reaction", () => {
    expect(validateReactions(ammonia(), EQUIL, COMPS)).toEqual([]);
  });

  it("flags a missing product", () => {
    const d = emptyDraft();
    d.reactants = [{ comp: "nitrogen", coeff: 1 }];
    d.products = [];
    expect(validateReactions([d], EQUIL, COMPS).join(" ")).toMatch(/at least one product/);
  });

  it("flags an unknown component", () => {
    const drafts = normalizeReactions({ reaction: { stoich: { unobtainium: -1, ammonia: 1 } } });
    expect(validateReactions(drafts, EQUIL, COMPS).join(" ")).toMatch(/not a flowsheet component/);
  });

  it("requires a key reactant for ConversionReactor", () => {
    const drafts = normalizeReactions({ reaction: { stoich: { nitrogen: -1, ammonia: 1 } }, conversion: 0.8 });
    drafts[0].key = "";
    expect(validateReactions(drafts, CONV, COMPS).join(" ")).toMatch(/key reactant/);
  });

  it("rejects a key that is not a reactant", () => {
    const drafts = ammonia();
    drafts[0].key = "ammonia";
    expect(validateReactions(drafts, CONV, COMPS).join(" ")).toMatch(/must be a reactant/);
  });

  it("bounds conversion to (0,1]", () => {
    const drafts = ammonia();
    drafts[0].conversion = 1.5;
    expect(validateReactions(drafts, CONV, COMPS).join(" ")).toMatch(/conversion must be/);
  });

  it("requires positive k0 for kinetic", () => {
    const d = emptyDraft();
    d.reactants = [{ comp: "nitrogen", coeff: 1 }];
    d.products = [{ comp: "ammonia", coeff: 1 }];
    d.key = "nitrogen";
    d.k0 = 0;
    expect(validateReactions([d], KIN, COMPS).join(" ")).toMatch(/k0 must be/);
  });

  it("flags a component that cancels across both sides", () => {
    const d = emptyDraft();
    d.reactants = [{ comp: "nitrogen", coeff: 1 }];
    d.products = [{ comp: "nitrogen", coeff: 1 }];
    expect(validateReactions([d], EQUIL, COMPS).join(" ")).toMatch(/cancels/);
  });
});
