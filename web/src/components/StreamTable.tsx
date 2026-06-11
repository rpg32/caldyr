import { Download, Scale } from "lucide-react";
import {
  CartesianGrid, Line, LineChart, ResponsiveContainer, Tooltip, XAxis, YAxis,
} from "recharts";
import { downloadCsv } from "../lib/csv";
import { convert, fmtQty, unitOf, UNIT_SETS, type UnitSet } from "../lib/units";
import { useStore } from "../store";
import { Badge, Button, Hint, PanelTitle, StaleNotice } from "./ui";

function UnitSetPicker() {
  const unitSet = useStore((s) => s.unitSet);
  const setUnitSet = useStore((s) => s.setUnitSet);
  return (
    <label className="ml-auto flex items-center gap-1 text-[11px] text-muted">
      units
      <select
        className="rounded-md border border-line bg-panel2 px-1.5 py-0.5 text-text"
        value={unitSet}
        onChange={(e) => setUnitSet(e.target.value as UnitSet)}
      >
        {UNIT_SETS.map((u) => <option key={u} value={u}>{u}</option>)}
      </select>
    </label>
  );
}

function ConvergencePlot({ history }: { history: number[] }) {
  if (history.length < 2) return null;
  const data = history.map((r, i) => ({ iter: i + 1, residual: r }));
  return (
    <>
      <PanelTitle>Convergence (residual per iteration)</PanelTitle>
      <div style={{ height: 130 }}>
        <ResponsiveContainer>
          <LineChart data={data} margin={{ top: 4, right: 8, bottom: 0, left: 0 }}>
            <CartesianGrid stroke="var(--line)" strokeDasharray="3 3" />
            <XAxis dataKey="iter" tick={{ fill: "var(--muted)", fontSize: 10 }}
              stroke="var(--line)" />
            <YAxis scale="log" domain={["auto", "auto"]}
              tick={{ fill: "var(--muted)", fontSize: 10 }} stroke="var(--line)"
              tickFormatter={(v: number) => v.toExponential(0)} width={44} />
            <Tooltip
              contentStyle={{
                background: "var(--panel)", border: "1px solid var(--line)",
                borderRadius: 6, fontSize: 11,
              }}
              labelStyle={{ color: "var(--muted)" }}
              formatter={(v) => [Number(v).toExponential(2), "residual"]}
            />
            <Line type="monotone" dataKey="residual" stroke="var(--accent)"
              strokeWidth={1.5} dot={false} isAnimationActive={false} />
          </LineChart>
        </ResponsiveContainer>
      </div>
    </>
  );
}

export function StreamTable() {
  const res = useStore((s) => s.solveRes);
  const stale = useStore((s) => s.resultsStale);
  const unitSet = useStore((s) => s.unitSet);
  if (!res) return <Hint>Press Solve to compute the stream table.</Hint>;
  const streams = Object.values(res.streams).filter((s) => s.molar_flow != null);

  const exportCsv = () => downloadCsv(
    [
      ["stream", `T (${unitOf("T", unitSet)})`, `P (${unitOf("P", unitSet)})`,
       `n (${unitOf("flow", unitSet)})`, "phase", "vapor_fraction",
       ...Object.keys(streams[0]?.z ?? {}).map((c) => `z_${c}`)],
      ...streams.map((s) => [
        s.id,
        s.T != null ? convert("T", s.T, unitSet) : null,
        s.P != null ? convert("P", s.P, unitSet) : null,
        s.molar_flow != null ? convert("flow", s.molar_flow, unitSet) : null,
        s.phase, s.vapor_fraction,
        ...Object.values(s.z ?? {}),
      ]),
    ],
    "streams.csv",
  );

  return (
    <div>
      <StaleNotice stale={stale} />
      <div className="flex items-center gap-2">
        <Badge ok={res.report.converged}>
          {res.report.converged ? "converged" : "not converged"} · {res.report.method}
          {res.report.iterations ? ` · ${res.report.iterations} iters` : ""}
        </Badge>
        <UnitSetPicker />
        <button
          className="cursor-pointer rounded-md border border-line bg-panel2 p-1 text-muted hover:border-accent"
          onClick={exportCsv} title="Export stream table as CSV" aria-label="Export CSV"
        >
          <Download size={12} />
        </button>
      </div>
      <table className="data-table">
        <thead>
          <tr>
            <th>stream</th>
            <th>T / {unitOf("T", unitSet)}</th>
            <th>P / {unitOf("P", unitSet)}</th>
            <th>n / {unitOf("flow", unitSet)}</th>
            <th>phase</th>
          </tr>
        </thead>
        <tbody>
          {streams.map((s) => (
            <tr key={s.id}>
              <td>{s.id}</td>
              <td>{fmtQty("T", s.T, unitSet)}</td>
              <td>{fmtQty("P", s.P, unitSet)}</td>
              <td>{fmtQty("flow", s.molar_flow, unitSet, 3)}</td>
              <td>{s.phase ?? "—"}</td>
            </tr>
          ))}
        </tbody>
      </table>
      {Object.keys(res.report.duties).length > 0 && (
        <table className="data-table mt-3">
          <thead><tr><th>duty</th><th>{unitOf("power", unitSet)}</th></tr></thead>
          <tbody>
            {Object.entries(res.report.duties).map(([k, v]) => (
              <tr key={k}><td>{k}</td><td>{fmtQty("power", v, unitSet)}</td></tr>
            ))}
          </tbody>
        </table>
      )}
      {res.report.tear_streams.length > 0 && (
        <div className="mt-2 text-[11px] text-muted">
          tear streams: {res.report.tear_streams.join(", ")}
        </div>
      )}
      <ConvergencePlot history={res.report.history ?? []} />
      <BalanceSection />
    </div>
  );
}

function BalanceSection() {
  const balance = useStore((s) => s.balance);
  const busy = useStore((s) => s.balanceBusy);
  const runBalance = useStore((s) => s.runBalance);
  const fmtRel = (x: number) => (x < 1e-12 ? "0" : x.toExponential(1));
  return (
    <div className="mt-3">
      <PanelTitle>Mass & energy closure</PanelTitle>
      {!balance ? (
        <Button icon={<Scale size={13} />} busy={busy} onClick={() => void runBalance()}>
          Balance check
        </Button>
      ) : (
        <>
          <table className="data-table">
            <thead>
              <tr><th>unit</th><th>mass rel</th><th>energy rel</th><th>duty kW</th></tr>
            </thead>
            <tbody>
              {balance.units.map((u) => (
                <tr key={u.unit_id}>
                  <td>{u.unit_id}</td>
                  <td>{fmtRel(u.mass_rel)}</td>
                  <td>{fmtRel(u.energy_rel)}</td>
                  <td>{(u.duty_W / 1000).toFixed(1)}</td>
                </tr>
              ))}
              <tr className="font-semibold">
                <td>overall</td>
                <td>{fmtRel(balance.overall.mass_rel)}</td>
                <td>{fmtRel(balance.overall.energy_rel)}</td>
                <td>{(balance.overall.duty_W / 1000).toFixed(1)}</td>
              </tr>
            </tbody>
          </table>
          {balance.warnings.map((w, i) => (
            <div key={i} className="mt-1 text-[11px] text-warn">⚠ {w}</div>
          ))}
        </>
      )}
    </div>
  );
}
