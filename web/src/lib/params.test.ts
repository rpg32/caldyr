// Param-metadata + conditional-applicability tests — the typed/decant param UI.
import { describe, expect, it } from "vitest";
import { PARAM_META, metaFor, paramApplies } from "./params";

describe("decant-condenser param metadata", () => {
  it("types the integrated-decant params (boolean / number / select)", () => {
    expect(PARAM_META.decant_condenser.type).toBe("boolean");
    expect(metaFor("condenser_T").unit).toBe("K");
    expect(PARAM_META.condenser_T.type).toBeUndefined(); // number is the default
    expect(PARAM_META.reflux_layer.type).toBe("select");
    expect(PARAM_META.reflux_layer.options).toEqual(["organic", "aqueous"]);
    // the solver method is a select including the one decant mode requires
    expect(PARAM_META.method.options).toContain("naphtali_sandholm");
  });

  it("gates condenser_T / reflux_layer on decant_condenser being on", () => {
    // applies only when the decant condenser is enabled
    expect(paramApplies("condenser_T", { decant_condenser: true })).toBe(true);
    expect(paramApplies("reflux_layer", { decant_condenser: true })).toBe(true);
    // hidden otherwise (absent or explicitly false)
    expect(paramApplies("condenser_T", {})).toBe(false);
    expect(paramApplies("condenser_T", { decant_condenser: false })).toBe(false);
    expect(paramApplies("reflux_layer", {})).toBe(false);
    // params without a `requires` predicate always apply
    expect(paramApplies("reflux_ratio", {})).toBe(true);
    expect(paramApplies("decant_condenser", {})).toBe(true);
  });

  it("falls back to a bare label for unknown keys", () => {
    expect(metaFor("totally_unknown_param")).toEqual({
      label: "totally_unknown_param", unit: "",
    });
  });
});
