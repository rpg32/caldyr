// Pure helpers for the NumberInput draft logic (kept here so they're testable
// without a DOM). A controlled `type="number"` input cannot hold the partial
// strings a user types on the way to a value (".", ".5", "5.", "1e-3"); it
// mangles ".5" into "05". We instead keep a text draft, validate keystrokes
// with `isNumericDraft`, and commit only when `parseFloat` is finite.

// Allows: "", "-", "+", ".", ".5", "0.5", "5.", "12", "1e", "1e-", "1.5e3".
// Rejects: letters, multiple dots, stray symbols.
const NUMERIC_DRAFT = /^[-+]?(\d*\.?\d*)(e[-+]?\d*)?$/i;

/** Whether `text` is a valid in-progress numeric draft worth keeping. */
export function isNumericDraft(text: string): boolean {
  return text === "" || NUMERIC_DRAFT.test(text);
}

/** Finite numeric value to commit from a draft, or null if not yet a number. */
export function commitNumericDraft(text: string): number | null {
  const n = parseFloat(text);
  return Number.isFinite(n) ? n : null;
}
