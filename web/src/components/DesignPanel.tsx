// Per-unit design results from the engine (unit.design): scalar table for
// things like FUG numbers / fuel duty, and stage-profile charts for columns.
import { useMemo } from "react";
import {
  CartesianGrid, Legend, Line, LineChart, ResponsiveContainer, Tooltip,
  XAxis, YAxis,
} from "recharts";
import { defaultUnit, toDisplay } from "../lib/units";
import { useStore } from "../store";
import { PanelTitle } from "./ui";

const TOOLTIP_STYLE = {
  background: "var(--panel)", border: "1px solid var(--line)",
  borderRadius: 6, fontSize: 11,
} as const;

const LINE_COLORS = ["#38bdf8", "#f97316", "#22c55e", "#a855f7", "#ef4444",
  "#eab308", "#14b8a6", "#f472b6"];

const isNumArray = (v: unknown): v is number[] =>
  Array.isArray(v) && v.length > 1 && v.every((x) => typeof x === "number");
const isMatrix = (v: unknown): v is number[][] =>
  Array.isArray(v) && v.length > 1 && v.every(isNumArray);
const isDictRows = (v: unknown): v is Record<string, number>[] =>
  Array.isArray(v) && v.length > 1 && v.every(
    (row) => row !== null && typeof row === "object" && !Array.isArray(row)
      && Object.values(row as object).every((x) => typeof x === "number"));

function StageChart({ data, lines, yLabel }: {
  data: Record<string, number>[];
  lines: { key: string; color: string }[];
  yLabel: string;
}) {
  return (
    <div style={{ height: 160 }}>
      <ResponsiveContainer>
        <LineChart data={data} margin={{ top: 4, right: 10, left: 0, bottom: 0 }}>
          <CartesianGrid stroke="var(--line)" strokeDasharray="3 3" />
          <XAxis dataKey="stage" tick={{ fill: "var(--muted)", fontSize: 10 }}
            stroke="var(--line)" />
          <YAxis tick={{ fill: "var(--muted)", fontSize: 10 }} stroke="var(--line)"
            width={48} domain={["auto", "auto"]}
            tickFormatter={(v: number) => Number(v).toPrecision(3)}
            label={{ value: yLabel, fill: "var(--muted)", fontSize: 10,
                     angle: -90, position: "insideLeft" }} />
          <Tooltip contentStyle={TOOLTIP_STYLE}
            labelFormatter={(v) => `stage ${v}`} />
          {lines.length > 1 && <Legend wrapperStyle={{ fontSize: 10 }} />}
          {lines.map((l) => (
            <Line key={l.key} dataKey={l.key} stroke={l.color} strokeWidth={1.5}
              dot={false} isAnimationActive={false} />
          ))}
        </LineChart>
      </ResponsiveContainer>
    </div>
  );
}

export function DesignPanel({ unitId }: { unitId: string }) {
  const design = useStore((s) => s.solveRes?.designs?.[unitId]);
  const components = useStore((s) => s.components);
  const unitSet = useStore((s) => s.unitSet);

  const { scalars, tProfile, xProfile } = useMemo(() => {
    const scalars: [string, number][] = [];
    let tProfile: Record<string, number>[] | null = null;
    let xProfile: { data: Record<string, number>[]; keys: string[] } | null = null;
    if (design) {
      for (const [k, v] of Object.entries(design)) {
        if (typeof v === "number" && Number.isFinite(v)) scalars.push([k, v]);
      }
      const T = design.T_profile ?? design.T;
      if (isNumArray(T)) {
        tProfile = T.map((t, i) => ({ stage: i + 1, T: toDisplay("temperature", t, unitSet) }));
      }
      const x = design.x_profile ?? design.x;
      if (isDictRows(x)) {
        // list of {component: fraction} per stage — keys carry the names
        const keys = Object.keys(x[0]);
        xProfile = {
          keys,
          data: x.map((row, i) => ({ stage: i + 1, ...row })),
        };
      } else if (isMatrix(x)) {
        const nComp = x[0].length;
        const keys = Array.from({ length: nComp },
          (_, c) => components[c] ?? `comp ${c + 1}`);
        xProfile = {
          keys,
          data: x.map((row, i) => ({
            stage: i + 1,
            ...Object.fromEntries(keys.map((key, c) => [key, row[c]])),
          })),
        };
      }
    }
    return { scalars, tProfile, xProfile };
  }, [design, components, unitSet]);

  if (!design) return null;

  return (
    <div className="mt-2">
      <PanelTitle>Design results</PanelTitle>
      {scalars.length > 0 && (
        <table className="data-table">
          <tbody>
            {scalars.map(([k, v]) => (
              <tr key={k}>
                <td>{k}</td>
                <td>{Math.abs(v) >= 1e5 || (Math.abs(v) < 1e-3 && v !== 0)
                  ? v.toExponential(3) : v.toLocaleString(undefined, { maximumFractionDigits: 4 })}</td>
              </tr>
            ))}
          </tbody>
        </table>
      )}
      {tProfile && (
        <>
          <PanelTitle>Temperature profile</PanelTitle>
          <StageChart data={tProfile} lines={[{ key: "T", color: "var(--accent)" }]}
            yLabel={`T / ${defaultUnit("temperature", unitSet)}`} />
        </>
      )}
      {xProfile && (
        <>
          <PanelTitle>Liquid composition profile</PanelTitle>
          <StageChart data={xProfile.data}
            lines={xProfile.keys.map((key, i) => ({
              key, color: LINE_COLORS[i % LINE_COLORS.length],
            }))}
            yLabel="x" />
        </>
      )}
    </div>
  );
}
