// The HYSYS-Spreadsheet analog: named expressions evaluated over the solved
// flowsheet. Expressions call accessor functions — T("S1"), z("S1","benzene"),
// duty("RXN_duty") — plus Math.* and earlier rows by name, so no identifier
// rewriting is needed and evaluation stays sandboxed to a strict whitelist.
import type { SolveResponse } from "../types";

export interface CalcRow {
  name: string;
  expr: string;
}

export interface CalcValue {
  name: string;
  value: number | null;
  error: string | null;
}

const NAME_RE = /^[A-Za-z_][A-Za-z0-9_]*$/;
const ACCESSORS = new Set(["T", "P", "n", "VF", "H", "z", "duty", "Math"]);
const FORBIDDEN = /(__|=>|[;=`]|\b(constructor|prototype|window|globalThis|Function|eval|import|this|document|fetch)\b)/;

/** Tokens left after removing string literals must all be known. */
function validate(expr: string, known: Set<string>): string | null {
  if (FORBIDDEN.test(expr)) return "expression contains forbidden syntax";
  const noStrings = expr.replace(/"[^"]*"|'[^']*'/g, "");
  if (/["']/.test(noStrings)) return "unterminated string literal";
  const idents = noStrings.match(/[A-Za-z_][A-Za-z0-9_]*/g) ?? [];
  for (const raw of idents) {
    // allow Math.sqrt etc. — the member after Math. is fine
    if (ACCESSORS.has(raw) || known.has(raw)) continue;
    if (/^\d/.test(raw)) continue;
    // Math member names appear as separate idents (e.g. "sqrt" in Math.sqrt)
    if (typeof (Math as unknown as Record<string, unknown>)[raw] !== "undefined") continue;
    return `unknown name "${raw}" — use T/P/n/VF/H("stream"), z("stream","comp"), duty("id"), Math.*, or an earlier row`;
  }
  if (/[^A-Za-z0-9_."'()+\-*/%,.\s<>?:!&|]/.test(expr)) {
    return "expression contains unsupported characters";
  }
  return null;
}

export function evaluateCalcs(rows: CalcRow[], solve: SolveResponse | null): CalcValue[] {
  const out: CalcValue[] = [];
  const prior: Record<string, number> = {};

  const streams = solve?.streams ?? {};
  const duties = solve?.report.duties ?? {};
  const get = (sid: string) => {
    const s = streams[sid];
    if (!s) throw new Error(`no solved stream "${sid}" — solve first?`);
    return s;
  };
  const scope = {
    T: (sid: string) => get(sid).T ?? NaN,
    P: (sid: string) => get(sid).P ?? NaN,
    n: (sid: string) => get(sid).molar_flow ?? NaN,
    VF: (sid: string) => get(sid).vapor_fraction ?? NaN,
    H: (sid: string) => get(sid).H ?? NaN,
    z: (sid: string, comp: string) => {
      const zz = get(sid).z ?? {};
      const total = Object.values(zz).reduce((a, b) => a + b, 0) || 1;
      return (zz[comp] ?? 0) / total;
    },
    duty: (id: string) => {
      if (!(id in duties)) throw new Error(`no duty "${id}" (known: ${Object.keys(duties).join(", ")})`);
      return duties[id];
    },
  };

  for (const row of rows) {
    const name = row.name.trim();
    if (!NAME_RE.test(name)) {
      out.push({ name: row.name, value: null, error: "name must be a simple identifier" });
      continue;
    }
    const known = new Set(Object.keys(prior));
    const invalid = validate(row.expr, known);
    if (invalid) {
      out.push({ name, value: null, error: invalid });
      continue;
    }
    try {
      const argNames = ["T", "P", "n", "VF", "H", "z", "duty", "Math", ...Object.keys(prior)];
      const argVals = [scope.T, scope.P, scope.n, scope.VF, scope.H, scope.z, scope.duty,
        Math, ...Object.values(prior)];
      // eslint-disable-next-line no-new-func
      const fn = new Function(...argNames, `"use strict"; return (${row.expr});`);
      const v = fn(...argVals);
      if (typeof v !== "number" || !Number.isFinite(v)) {
        out.push({ name, value: null, error: "did not evaluate to a finite number" });
      } else {
        prior[name] = v;
        out.push({ name, value: v, error: null });
      }
    } catch (e) {
      out.push({ name, value: null, error: (e as Error).message });
    }
  }
  return out;
}
