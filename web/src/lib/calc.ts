// The HYSYS-Spreadsheet analog: named expressions evaluated over the solved
// flowsheet. Expressions call accessor functions — T("S1"), z("S1","benzene"),
// duty("RXN_duty") — plus Math.* and earlier rows by name.
//
// Expressions arrive in shared .flow files, so they are untrusted input:
// evaluation is a small recursive-descent interpreter over a whitelisted
// grammar — no dynamic code generation (new Function/eval) anywhere.
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

// -- expression interpreter ---------------------------------------------------
// Grammar (JS-compatible subset): ternary ?:, || &&, < >, + -, * / %, ** ,
// unary ! - +, calls f(a, b), member access on Math only, numbers, single- or
// double-quoted strings, and whitelisted identifiers. Branches of ?: and the
// right side of &&/|| stay lazy, matching JS semantics.

type Token =
  | { kind: "num"; value: number }
  | { kind: "str"; value: string }
  | { kind: "ident"; value: string }
  | { kind: "punct"; value: string };

function tokenize(src: string): Token[] {
  const out: Token[] = [];
  let i = 0;
  while (i < src.length) {
    const c = src[i];
    if (/\s/.test(c)) { i++; continue; }
    if (/[0-9]/.test(c) || (c === "." && /[0-9]/.test(src[i + 1] ?? ""))) {
      let j = i;
      while (j < src.length && /[0-9.]/.test(src[j])) j++;
      const value = Number(src.slice(i, j));
      if (!Number.isFinite(value)) throw new Error(`bad number "${src.slice(i, j)}"`);
      out.push({ kind: "num", value });
      i = j;
      continue;
    }
    if (c === '"' || c === "'") {
      const close = src.indexOf(c, i + 1);
      if (close < 0) throw new Error("unterminated string literal");
      out.push({ kind: "str", value: src.slice(i + 1, close) });
      i = close + 1;
      continue;
    }
    if (/[A-Za-z_]/.test(c)) {
      let j = i;
      while (j < src.length && /[A-Za-z0-9_]/.test(src[j])) j++;
      out.push({ kind: "ident", value: src.slice(i, j) });
      i = j;
      continue;
    }
    const two = src.slice(i, i + 2);
    if (two === "**" || two === "&&" || two === "||") {
      out.push({ kind: "punct", value: two });
      i += 2;
      continue;
    }
    if ("(),.?:!<>+-*/%".includes(c)) {
      out.push({ kind: "punct", value: c });
      i += 1;
      continue;
    }
    throw new Error(`unsupported character "${c}"`);
  }
  return out;
}

type Thunk = () => unknown;
const num = (v: unknown): number => v as number;

function parseExpr(src: string, env: Record<string, unknown>): Thunk {
  const tokens = tokenize(src);
  let pos = 0;
  const isPunct = (v: string): boolean => {
    const t = tokens[pos];
    return t?.kind === "punct" && t.value === v;
  };
  const eat = (v: string): void => {
    if (!isPunct(v)) throw new Error(`expected "${v}"`);
    pos++;
  };

  function primary(): Thunk {
    const t = tokens[pos];
    if (!t) throw new Error("unexpected end of expression");
    if (t.kind === "num") { pos++; return () => t.value; }
    if (t.kind === "str") { pos++; return () => t.value; }
    if (t.kind === "ident") {
      pos++;
      return () => {
        if (Object.prototype.hasOwnProperty.call(env, t.value)) return env[t.value];
        throw new Error(`unknown name "${t.value}"`);
      };
    }
    if (t.value === "(") {
      pos++;
      const inner = ternary();
      eat(")");
      return inner;
    }
    throw new Error(`unexpected "${t.value}"`);
  }

  function postfix(): Thunk {
    let node = primary();
    for (;;) {
      if (isPunct(".")) {
        pos++;
        const t = tokens[pos];
        if (t?.kind !== "ident") throw new Error("expected a name after '.'");
        pos++;
        const base = node;
        node = () => {
          if (base() !== Math || !Object.prototype.hasOwnProperty.call(Math, t.value)) {
            throw new Error(`".${t.value}" — member access is only supported on Math`);
          }
          return (Math as unknown as Record<string, unknown>)[t.value];
        };
      } else if (isPunct("(")) {
        pos++;
        const args: Thunk[] = [];
        if (!isPunct(")")) {
          for (;;) {
            args.push(ternary());
            if (!isPunct(",")) break;
            pos++;
          }
        }
        eat(")");
        const callee = node;
        node = () => {
          const fn = callee();
          if (typeof fn !== "function") throw new Error("not a function");
          return (fn as (...a: unknown[]) => unknown)(...args.map((a) => a()));
        };
      } else {
        return node;
      }
    }
  }

  function unary(): Thunk {
    if (isPunct("!")) { pos++; const u = unary(); return () => !u(); }
    if (isPunct("-")) { pos++; const u = unary(); return () => -num(u()); }
    if (isPunct("+")) { pos++; const u = unary(); return () => +num(u()); }
    return power();
  }

  function power(): Thunk {
    const base = postfix();
    if (!isPunct("**")) return base;
    pos++;
    const exp = unary();
    return () => num(base()) ** num(exp());
  }

  function multiplicative(): Thunk {
    let node = unary();
    for (;;) {
      const op = isPunct("*") ? "*" : isPunct("/") ? "/" : isPunct("%") ? "%" : null;
      if (!op) return node;
      pos++;
      const l = node;
      const r = unary();
      node = op === "*" ? () => num(l()) * num(r())
        : op === "/" ? () => num(l()) / num(r())
        : () => num(l()) % num(r());
    }
  }

  function additive(): Thunk {
    let node = multiplicative();
    for (;;) {
      const op = isPunct("+") ? "+" : isPunct("-") ? "-" : null;
      if (!op) return node;
      pos++;
      const l = node;
      const r = multiplicative();
      node = op === "+" ? () => num(l()) + num(r()) : () => num(l()) - num(r());
    }
  }

  function comparison(): Thunk {
    let node = additive();
    for (;;) {
      const op = isPunct("<") ? "<" : isPunct(">") ? ">" : null;
      if (!op) return node;
      pos++;
      const l = node;
      const r = additive();
      node = op === "<" ? () => num(l()) < num(r()) : () => num(l()) > num(r());
    }
  }

  function logicalAnd(): Thunk {
    let node = comparison();
    while (isPunct("&&")) {
      pos++;
      const l = node;
      const r = comparison();
      node = () => { const v = l(); return v ? r() : v; };
    }
    return node;
  }

  function logicalOr(): Thunk {
    let node = logicalAnd();
    while (isPunct("||")) {
      pos++;
      const l = node;
      const r = logicalAnd();
      node = () => { const v = l(); return v ? v : r(); };
    }
    return node;
  }

  function ternary(): Thunk {
    const cond = logicalOr();
    if (!isPunct("?")) return cond;
    pos++;
    const yes = ternary();
    eat(":");
    const no = ternary();
    return () => (cond() ? yes() : no());
  }

  const root = ternary();
  if (pos < tokens.length) {
    const t = tokens[pos];
    throw new Error(`unexpected "${t.kind === "punct" ? t.value : ("value" in t ? t.value : t)}"`);
  }
  return root;
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
      const env: Record<string, unknown> = { ...scope, Math, ...prior };
      const v = parseExpr(row.expr, env)();
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
