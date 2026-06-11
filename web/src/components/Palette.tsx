import { ArrowRightFromLine, ArrowRightToLine } from "lucide-react";
import { useStore } from "../store";
import { PanelTitle } from "./ui";

export function Palette() {
  const unitTypes = useStore((s) => s.unitTypes);
  const addUnit = useStore((s) => s.addUnit);
  const item =
    "w-full rounded-md border border-line bg-panel2 px-2 py-1.5 text-left text-[13px] " +
    "cursor-pointer hover:border-accent transition-colors";
  return (
    <aside className="flex flex-col gap-1 overflow-auto border-r border-line bg-panel p-2">
      <PanelTitle>Boundaries</PanelTitle>
      <button className={item} title="Material feed with a T/P/flow/composition spec"
        onClick={() => addUnit("feed")}>
        <ArrowRightFromLine size={12} className="mr-1 inline text-ok" aria-hidden />
        Feed
      </button>
      <button className={item} title="Product sink (boundary outlet)"
        onClick={() => addUnit("product")}>
        <ArrowRightToLine size={12} className="mr-1 inline text-muted" aria-hidden />
        Product
      </button>

      <PanelTitle>Unit ops</PanelTitle>
      {unitTypes.map((t) => (
        <button key={t.type} className={item} title={t.doc} onClick={() => addUnit("unit", t)}>
          {t.type}
        </button>
      ))}
    </aside>
  );
}
