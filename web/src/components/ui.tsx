// Small presentational primitives shared across the app.
import { Loader2 } from "lucide-react";
import { useState, type ButtonHTMLAttributes, type InputHTMLAttributes, type ReactNode } from "react";
import { commitNumericDraft, isNumericDraft } from "../lib/number";

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
