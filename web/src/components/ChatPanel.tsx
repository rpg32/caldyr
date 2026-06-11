// AI chat sidebar: converse with the flowsheet copilot (local LLM via the
// engine's caldyr.ai layer), watch its tool calls stream by, and review every
// proposed flowsheet edit as a diff before it touches the canvas.
import { Bot, Check, Loader2, Send, Sparkles, Wrench, X } from "lucide-react";
import { useEffect, useRef, useState } from "react";
import type { FlowDiff } from "../lib/diff";
import { useStore } from "../store";
import { Button } from "./ui";

function DiffCard({ diff }: { diff: FlowDiff }) {
  const accept = useStore((s) => s.acceptPending);
  const reject = useStore((s) => s.rejectPending);
  const Row = ({ children }: { children: React.ReactNode }) => (
    <li className="my-0.5">{children}</li>
  );
  return (
    <div className="my-2 rounded-lg border border-accent/50 bg-accent/8 p-2.5 text-[12px]">
      <div className="mb-1 font-semibold text-accent">
        Proposed flowsheet edit — {diff.count} change{diff.count > 1 ? "s" : ""}
      </div>
      <ul className="mb-2 list-none pl-0">
        {diff.addedUnits.map((u) => (
          <Row key={`au${u.id}`}><span className="text-ok">+ unit</span> {u.id} ({u.type})</Row>
        ))}
        {diff.removedUnits.map((id) => (
          <Row key={`ru${id}`}><span className="text-bad">− unit</span> {id}</Row>
        ))}
        {diff.changedParams.map((c, i) => (
          <Row key={`cp${i}`}>
            <span className="text-warn">~</span> {c.unit}.{c.key}:{" "}
            {JSON.stringify(c.from)} → {JSON.stringify(c.to)}
          </Row>
        ))}
        {diff.addedStreams.map((id) => (
          <Row key={`as${id}`}><span className="text-ok">+ stream</span> {id}</Row>
        ))}
        {diff.removedStreams.map((id) => (
          <Row key={`rs${id}`}><span className="text-bad">− stream</span> {id}</Row>
        ))}
        {diff.rewiredStreams.map((id) => (
          <Row key={`ws${id}`}><span className="text-warn">~ rewired</span> {id}</Row>
        ))}
        {diff.addedComponents.map((c) => (
          <Row key={`ac${c}`}><span className="text-ok">+ component</span> {c}</Row>
        ))}
        {diff.removedComponents.map((c) => (
          <Row key={`rc${c}`}><span className="text-bad">− component</span> {c}</Row>
        ))}
        {diff.packageChanged && (
          <Row><span className="text-warn">~ package</span>{" "}
            {diff.packageChanged.from} → {diff.packageChanged.to}</Row>
        )}
        {diff.logicalChanged && <Row><span className="text-warn">~ logical ops changed</span></Row>}
      </ul>
      <div className="flex gap-2">
        <Button variant="primary" icon={<Check size={13} />} onClick={accept}>Accept</Button>
        <Button variant="ghost" icon={<X size={13} />} onClick={reject}>Reject</Button>
      </div>
    </div>
  );
}

export function ChatPanel() {
  const open = useStore((s) => s.chatOpen);
  const messages = useStore((s) => s.chatMessages);
  const busy = useStore((s) => s.chatBusy);
  const pendingDiff = useStore((s) => s.pendingDiff);
  const sendChat = useStore((s) => s.sendChat);
  const quickAction = useStore((s) => s.quickAction);
  const toggleChat = useStore((s) => s.toggleChat);
  const [draft, setDraft] = useState("");
  const scrollRef = useRef<HTMLDivElement>(null);

  useEffect(() => {
    scrollRef.current?.scrollTo({ top: scrollRef.current.scrollHeight });
  }, [messages, pendingDiff, busy]);

  if (!open) return null;

  const send = () => {
    if (draft.trim()) {
      void sendChat(draft);
      setDraft("");
    }
  };

  return (
    <div className="absolute bottom-2 left-2 top-2 z-20 flex w-[330px] flex-col rounded-lg border border-line bg-panel/95 shadow-xl backdrop-blur">
      <div className="flex items-center gap-2 border-b border-line px-3 py-2">
        <Bot size={15} className="text-accent" aria-hidden />
        <b>Copilot</b>
        <span className="text-[11px] text-muted">local LLM</span>
        <button className="ml-auto cursor-pointer text-muted hover:text-text"
          onClick={toggleChat} aria-label="Close chat"><X size={14} /></button>
      </div>

      <div ref={scrollRef} className="min-h-0 flex-1 overflow-auto p-2.5">
        {messages.length === 0 && (
          <div className="p-2 text-[12px] leading-relaxed text-muted">
            Ask the copilot to build, edit, solve, cost, or explain the
            flowsheet — e.g. <em>"add a compressor before the reactor and
            re-solve"</em>. Edits come back as a diff you accept or reject.
          </div>
        )}
        {messages.map((m, i) =>
          m.role === "event" ? (
            <div key={i} className="my-1 flex items-center gap-1.5 text-[11px] text-muted">
              <Wrench size={10} aria-hidden /> {m.text}
            </div>
          ) : (
            <div key={i}
              className={`my-1.5 max-w-[92%] whitespace-pre-wrap rounded-lg px-2.5 py-1.5 text-[12.5px] leading-relaxed ${
                m.role === "user"
                  ? "ml-auto bg-accent/15 text-text"
                  : m.role === "error"
                    ? "border border-bad/40 bg-bad/10 text-bad"
                    : "bg-panel2 text-text"
              }`}>
              {m.text}
            </div>
          ))}
        {pendingDiff && <DiffCard diff={pendingDiff} />}
        {busy && (
          <div className="my-1 flex items-center gap-1.5 text-[11px] text-muted">
            <Loader2 size={11} className="animate-spin" aria-hidden /> thinking…
          </div>
        )}
      </div>

      <div className="border-t border-line p-2">
        <div className="mb-1.5 flex gap-1.5">
          <Button variant="ghost" icon={<Sparkles size={11} />}
            onClick={() => void quickAction("describe")}>Explain flowsheet</Button>
          <Button variant="ghost" icon={<Sparkles size={11} />}
            onClick={() => void quickAction("diagnose")}>Diagnose solve</Button>
        </div>
        <div className="flex items-center gap-1.5">
          <input
            className="min-w-0 flex-1 rounded-md border border-line bg-panel2 px-2 py-1.5 text-text"
            placeholder="Ask the copilot…"
            value={draft}
            aria-label="Chat message"
            onChange={(e) => setDraft(e.target.value)}
            onKeyDown={(e) => e.key === "Enter" && !e.shiftKey && send()}
            disabled={busy}
          />
          <Button variant="primary" icon={<Send size={13} />} disabled={busy || !draft.trim()}
            onClick={send} aria-label="Send" />
        </div>
      </div>
    </div>
  );
}
