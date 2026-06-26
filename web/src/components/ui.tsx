// Small presentational primitives shared across the app.
import { Loader2 } from "lucide-react";
import { useState, type ButtonHTMLAttributes, type InputHTMLAttributes, type ReactNode } from "react";
import { glossaryDef } from "../lib/glossary";
import { commitNumericDraft, isNumericDraft } from "../lib/number";
import { defaultUnit, toDisplay, toSI, type Dim, type UnitSet } from "../lib/units";

interface ButtonProps extends ButtonHTMLAttributes<HTMLButtonElement> {
  variant?: "default" | "primary" | "ghost";
  busy?: boolean;
  icon?: ReactNode;
}

export function Button({ variant = "default", busy, icon, children, className = "", ...rest }: ButtonProps) {
  const base =
    "inline-flex items-center gap-1.5 rounded-md border px-2.5 py-1.5 text-[13px] " +
    "cursor-pointer transition-colors disabled:opacity-50 disabled:cursor-not-allowed " +
    "focus-visible:outline-2 focus-visible:outline-accent focus-visible:outline-offset-1";
  const look =
    variant === "primary"
      ? "border-accent bg-accent/15 text-accent hover:bg-accent/25"
      : variant === "ghost"
        ? "border-transparent bg-transparent text-muted hover:text-text hover:border-line"
        : "border-line bg-panel2 text-text hover:border-accent";
  return (
    <button className={`${base} ${look} ${className}`} disabled={busy || rest.disabled} {...rest}>
      {busy ? <Loader2 size={14} className="animate-spin" aria-hidden /> : icon}
      {children}
    </button>
  );
}

interface NumberInputProps
  extends Omit<InputHTMLAttributes<HTMLInputElement>, "value" | "onChange" | "type"> {
  value: number;
  onChange: (v: number) => void;
}

/** Numeric text field that accepts leading-dot / trailing-dot / exponent input
 *  mid-typing and commits only finite parses. Drop-in for `type="number"`. */
export function NumberInput({ value, onChange, onBlur, ...rest }: NumberInputProps) {
  const [draft, setDraft] = useState<string | null>(null);
  const display = draft ?? (Number.isFinite(value) ? String(value) : "");
  return (
    <input
      {...rest}
      type="text"
      inputMode="decimal"
      value={display}
      onChange={(e) => {
        const t = e.target.value;
        if (!isNumericDraft(t)) return; // reject non-numeric keystrokes
        setDraft(t);
        const n = commitNumericDraft(t);
        if (n !== null) onChange(n);
      }}
      onBlur={(e) => { setDraft(null); onBlur?.(e); }}
    />
  );
}

interface QuantityInputProps
  extends Omit<InputHTMLAttributes<HTMLInputElement>, "value" | "onChange" | "type"> {
  dim: Dim;
  set: UnitSet;
  unit?: string;              // display-unit override; defaults to the set default
  value: number;              // SI base value
  onChange: (si: number) => void; // receives the SI base value
  // When provided, the unit label becomes a picker of these units (per-field
  // override). onUnitChange(null) means "follow the unit set's default".
  units?: string[];
  onUnitChange?: (unit: string | null) => void;
}

/** A numeric field for a physical quantity: edits in the chosen display unit
 *  (set default unless overridden), converts to/from the engine's SI base, and
 *  shows (or lets you pick) the unit. `value`/`onChange` are always SI. */
export function QuantityInput(
  { dim, set, unit, value, onChange, units, onUnitChange, ...rest }: QuantityInputProps,
) {
  const def = defaultUnit(dim, set);
  const u = unit ?? def;
  // toPrecision(12) tames float noise from non-exact factors (e.g. psia) so the
  // field doesn't show 14.696000001 after a round-trip; user edits override it.
  const disp = Number.isFinite(value) ? Number(toDisplay(dim, value, set, u).toPrecision(12)) : NaN;
  return (
    <span className="flex items-center gap-1.5">
      <NumberInput {...rest} value={disp}
        onChange={(d) => onChange(toSI(dim, d, set, u))} />
      {units && onUnitChange ? (
        <select
          className="w-[58px] shrink-0 rounded-md border border-line bg-panel2 px-1 py-0.5 text-[11px] text-muted"
          value={u} aria-label="Unit"
          title="Display unit for this field (the value is stored in SI)"
          // picking the unit-set default clears the override so the field follows the system
          onChange={(e) => onUnitChange(e.target.value === def ? null : e.target.value)}
        >
          {units.map((opt) => <option key={opt} value={opt}>{opt}</option>)}
        </select>
      ) : (
        <span className="w-[58px] shrink-0 whitespace-nowrap text-[11px] text-muted">{u}</span>
      )}
    </span>
  );
}

/** A numeric field that converts to/from the unit system when `dim` is known
 *  (e.g. a swept temperature param) and falls back to a plain SI number field
 *  for dimensionless / unknown quantities. `value`/`onChange` are always SI. */
export function DimField(
  { dim, set, unit, ...rest }: { dim?: Dim; set: UnitSet; unit?: string }
    & Omit<NumberInputProps, never>,
) {
  return dim
    ? <QuantityInput dim={dim} set={set} unit={unit} {...rest} />
    : <NumberInput {...rest} />;
}

/** Inline jargon term: dotted underline + a hover definition from the glossary.
 *  `k` overrides which glossary key to look up (defaults to the text content). */
export function Term({ k, children }: { k?: string; children: ReactNode }) {
  const key = k ?? (typeof children === "string" ? children : "");
  const def = glossaryDef(key);
  if (!def) return <>{children}</>;
  return (
    <abbr title={def}
      className="cursor-help underline decoration-dotted decoration-muted underline-offset-2">
      {children}
    </abbr>
  );
}

export function PanelTitle({ children }: { children: ReactNode }) {
  return (
    <div className="mb-1 mt-2 text-[11px] uppercase tracking-wider text-muted">{children}</div>
  );
}

export function Hint({ children }: { children: ReactNode }) {
  return <div className="p-3 leading-relaxed text-muted">{children}</div>;
}

export function Badge({ ok, children }: { ok: boolean; children: ReactNode }) {
  return (
    <span
      className={`inline-block rounded-full px-2 py-0.5 text-[11px] ${
        ok ? "bg-ok/15 text-ok" : "bg-bad/15 text-bad"
      }`}
    >
      {children}
    </span>
  );
}

export function StaleNotice({ stale }: { stale: boolean }) {
  if (!stale) return null;
  return (
    <div className="mb-2 rounded-md border border-warn/40 bg-warn/10 px-2 py-1 text-[12px] text-warn">
      Flowsheet edited since these results — re-run to refresh.
    </div>
  );
}
