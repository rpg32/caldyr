import { describe, expect, it } from "vitest";
import {
  DIMENSIONS, defaultUnit, fmtDim, toDisplay, toSI, unitsForDim, type Dim,
} from "./units";

const DIMS = Object.keys(DIMENSIONS) as Dim[];

describe("unit round-trips", () => {
  it("toSI(toDisplay(si)) === si for every unit of every dimension", () => {
    for (const dim of DIMS) {
      for (const unit of unitsForDim(dim)) {
        for (const si of [1, 42, 1234.5, 0.001]) {
          const disp = toDisplay(dim, si, "SI", unit);
          expect(toSI(dim, disp, "SI", unit), `${dim}/${unit}`).toBeCloseTo(si, 6);
        }
      }
    }
  });
});

describe("known HYSYS conversions", () => {
  it("temperature", () => {
    expect(toSI("temperature", 0, "SI", "°C")).toBeCloseTo(273.15, 6);
    expect(toSI("temperature", 32, "SI", "°F")).toBeCloseTo(273.15, 6);
    expect(toSI("temperature", 212, "SI", "°F")).toBeCloseTo(373.15, 6);
    expect(toSI("temperature", 491.67, "SI", "R")).toBeCloseTo(273.15, 4);
    expect(toDisplay("temperature", 373.15, "SI", "°C")).toBeCloseTo(100, 6);
  });

  it("pressure", () => {
    expect(toSI("pressure", 1, "SI", "atm")).toBeCloseTo(101325, 3);
    expect(toSI("pressure", 1, "SI", "bar")).toBeCloseTo(1e5, 3);
    expect(toSI("pressure", 14.696, "SI", "psia")).toBeCloseTo(101325, 0);
  });

  it("molar flow", () => {
    expect(toSI("molar_flow", 1, "SI", "kmol/h")).toBeCloseTo(0.277778, 5);
    expect(toSI("molar_flow", 1, "SI", "lbmol/h")).toBeCloseTo(0.125998, 5);
  });

  it("mass flow", () => {
    expect(toSI("mass_flow", 1, "SI", "lb/h")).toBeCloseTo(1.25998e-4, 9);
    expect(toSI("mass_flow", 3600, "SI", "kg/h")).toBeCloseTo(1, 9);
  });

  it("power / duty", () => {
    expect(toSI("power", 1, "SI", "kW")).toBeCloseTo(1000, 6);
    expect(toSI("power", 1, "SI", "MMBtu/h")).toBeCloseTo(293071.1, 0);
    expect(toSI("power", 1, "SI", "hp")).toBeCloseTo(745.6999, 3);
  });

  it("UA", () => {
    expect(toSI("UA", 1, "SI", "Btu/F-h")).toBeCloseTo(0.527527, 5);
  });
});

describe("set defaults + formatting", () => {
  it("picks the right default unit per set", () => {
    expect(defaultUnit("temperature", "Field")).toBe("°F");
    expect(defaultUnit("pressure", "Metric")).toBe("bar");
    expect(defaultUnit("molar_flow", "Field")).toBe("lbmol/h");
  });

  it("fmtDim converts and tolerates null", () => {
    expect(fmtDim("temperature", 373.15, "Metric")).toBe("100");
    expect(fmtDim("pressure", null, "SI")).toBe("—");
    expect(fmtDim("pressure", 700000, "SI")).toBe("700"); // kPa
  });
});
