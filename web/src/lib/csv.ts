// Tiny CSV builder + download (no dependency).

const cell = (v: unknown): string => {
  let s = String(v ?? "");
  // Ids/names can come from shared .flow files; a leading =+-@ or tab/CR would
  // execute as a formula when the CSV is opened in Excel/Sheets. Numbers are
  // exempt so negative values stay numeric.
  if (typeof v === "string" && /^[=+\-@\t\r]/.test(s)) s = `'${s}`;
  return /[",\n]/.test(s) ? `"${s.replace(/"/g, '""')}"` : s;
};

export function toCsv(rows: (string | number | null | undefined)[][]): string {
  return rows.map((r) => r.map(cell).join(",")).join("\n");
}

export function downloadCsv(rows: (string | number | null | undefined)[][], filename: string): void {
  const blob = new Blob([toCsv(rows)], { type: "text/csv" });
  const url = URL.createObjectURL(blob);
  const a = document.createElement("a");
  a.href = url;
  a.download = filename;
  a.click();
  URL.revokeObjectURL(url);
}
