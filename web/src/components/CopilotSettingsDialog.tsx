// Bring-your-own-LLM configuration for the copilot. The provider/model/key are
// saved SERVER-side (localhost API) — the key is sent once here and never kept
// in the browser. Supports local Ollama, any OpenAI-compatible server, and
// Anthropic (Claude).
import { Bot, CheckCircle2, X, XCircle } from "lucide-react";
import { useEffect, useState } from "react";
import { api } from "../api";
import { useStore } from "../store";
import { Button, PanelTitle } from "./ui";

type Provider = "ollama" | "openai" | "anthropic";

const PROVIDERS: { id: Provider; label: string; blurb: string }[] = [
  { id: "ollama", label: "Ollama (local)", blurb: "Runs on your machine — no API key, nothing leaves your computer." },
  { id: "openai", label: "OpenAI / compatible", blurb: "OpenAI, or a local server (LM Studio, vLLM, llama.cpp) via its base URL." },
  { id: "anthropic", label: "Anthropic (Claude)", blurb: "Claude models via your Anthropic API key." },
];

export function CopilotSettingsDialog() {
  const open = useStore((s) => s.copilotSettingsOpen);
  const toggle = useStore((s) => s.toggleCopilotSettings);
  const health = useStore((s) => s.copilotHealth);
  const saveConfig = useStore((s) => s.saveCopilotConfig);

  const [provider, setProvider] = useState<Provider>("ollama");
  const [model, setModel] = useState("");
  const [baseUrl, setBaseUrl] = useState("");
  const [apiKey, setApiKey] = useState("");
  const [test, setTest] = useState<{ ok: boolean; detail: string } | null>(null);
  const [busy, setBusy] = useState<"test" | "save" | null>(null);
  const [saved, setSaved] = useState(false);

  // Seed the form from the saved config whenever the dialog opens.
  useEffect(() => {
    if (!open) return;
    const c = health?.config;
    setProvider((c?.provider as Provider) || "ollama");
    setModel(c?.model ?? "");
    setBaseUrl(c?.base_url ?? "");
    setApiKey("");
    setTest(null);
    setSaved(false);
    const onKey = (e: KeyboardEvent) => e.key === "Escape" && toggle();
    window.addEventListener("keydown", onKey);
    return () => window.removeEventListener("keydown", onKey);
  }, [open, health, toggle]);

  if (!open) return null;

  const hasKey = health?.config.has_key ?? false;
  const ollamaModels = health?.providers.ollama.models ?? [];
  const anthropicReady = health?.providers.anthropic.available ?? false;
  const needsKey = provider === "openai" || provider === "anthropic";

  const update = () => ({
    provider,
    model: model.trim(),
    base_url: baseUrl.trim(),
    ...(apiKey.trim() ? { api_key: apiKey.trim() } : {}),
  });

  const runTest = async () => {
    setBusy("test");
    setTest(null);
    try {
      setTest(await api.aiTestConfig(update()));
    } catch (e) {
      setTest({ ok: false, detail: (e as Error).message });
    } finally {
      setBusy(null);
    }
  };

  const save = async () => {
    setBusy("save");
    try {
      await saveConfig(update());
      setSaved(true);
      setApiKey("");
    } finally {
      setBusy(null);
    }
  };

  const clearKey = async () => {
    setBusy("save");
    try {
      await saveConfig({ provider, clear_key: true });
      setApiKey("");
    } finally {
      setBusy(null);
    }
  };

  const field = "w-full rounded-md border border-line bg-panel2 px-2 py-1.5 text-text";

  return (
    <div className="fixed inset-0 z-40 flex items-center justify-center bg-black/40"
      onClick={toggle}>
      <div role="dialog" aria-modal="true" aria-label="Copilot LLM settings"
        className="max-h-[85vh] w-[520px] overflow-auto rounded-xl border border-line bg-panel p-4 shadow-2xl"
        onClick={(e) => e.stopPropagation()}>
        <div className="mb-1 flex items-center gap-2">
          <Bot size={16} className="text-accent" aria-hidden />
          <b className="text-[15px]">Copilot &mdash; AI provider</b>
          <button className="ml-auto cursor-pointer text-muted hover:text-text"
            onClick={toggle} aria-label="Close dialog"><X size={16} /></button>
        </div>
        <p className="mb-1 text-[12px] leading-snug text-muted">
          Choose which LLM powers the copilot. Your key is stored on this machine
          by the local Caldyr server &mdash; it is never saved in the browser.
        </p>

        <PanelTitle>Provider</PanelTitle>
        <div className="flex flex-col gap-1.5">
          {PROVIDERS.map((p) => (
            <label key={p.id}
              className={`flex cursor-pointer items-start gap-2 rounded-lg border p-2.5 ${
                provider === p.id ? "border-accent bg-accent/8" : "border-line bg-panel2"
              }`}>
              <input type="radio" name="provider" className="mt-1" checked={provider === p.id}
                onChange={() => { setProvider(p.id); setTest(null); setSaved(false); }} />
              <span className="min-w-0">
                <span className="flex items-center gap-2 font-semibold">
                  {p.label}
                  {p.id === "ollama" && (
                    <StatusChip ok={health?.providers.ollama.reachable}
                      okText="reachable" badText="not running" />
                  )}
                  {p.id === "anthropic" && (
                    <StatusChip ok={anthropicReady}
                      okText="SDK installed" badText="SDK missing" />
                  )}
                </span>
                <span className="block text-[11.5px] leading-snug text-muted">{p.blurb}</span>
              </span>
            </label>
          ))}
        </div>

        <PanelTitle>Model</PanelTitle>
        {provider === "ollama" && ollamaModels.length > 0 ? (
          <select className={field} value={model} onChange={(e) => setModel(e.target.value)}
            aria-label="Model">
            <option value="">(default)</option>
            {ollamaModels.map((m) => <option key={m} value={m}>{m}</option>)}
          </select>
        ) : (
          <input className={field} value={model} onChange={(e) => setModel(e.target.value)}
            placeholder={provider === "anthropic" ? "claude-sonnet-4-6" : provider === "openai" ? "gpt-4o-mini" : "model name"}
            aria-label="Model" />
        )}

        {(provider === "ollama" || provider === "openai" || provider === "anthropic") && (
          <>
            <PanelTitle>{provider === "ollama" ? "Host" : "Base URL"}{provider === "anthropic" ? " (optional)" : ""}</PanelTitle>
            <input className={field} value={baseUrl} onChange={(e) => setBaseUrl(e.target.value)}
              placeholder={provider === "ollama" ? "http://localhost:11434"
                : provider === "openai" ? "https://api.openai.com/v1 (or your local server)"
                : "(default Anthropic endpoint)"}
              aria-label="Base URL" />
          </>
        )}

        {needsKey && (
          <>
            <PanelTitle>API key</PanelTitle>
            <input className={field} type="password" value={apiKey}
              onChange={(e) => setApiKey(e.target.value)}
              placeholder={hasKey ? "•••••••••• (a key is saved — leave blank to keep it)" : "paste your API key"}
              aria-label="API key" autoComplete="off" />
            {hasKey && (
              <button className="mt-1 cursor-pointer text-[11px] text-muted hover:text-bad"
                onClick={() => void clearKey()}>Remove saved key</button>
            )}
          </>
        )}

        {test && (
          <div className={`mt-3 flex items-start gap-1.5 rounded-md border p-2 text-[12px] ${
            test.ok ? "border-ok/40 bg-ok/10 text-ok" : "border-bad/40 bg-bad/10 text-bad"
          }`}>
            {test.ok ? <CheckCircle2 size={14} className="mt-0.5 shrink-0" />
              : <XCircle size={14} className="mt-0.5 shrink-0" />}
            <span>{test.detail}</span>
          </div>
        )}

        <div className="mt-3 flex items-center gap-2">
          <Button variant="ghost" busy={busy === "test"} onClick={() => void runTest()}>
            Test connection
          </Button>
          <span className="ml-auto" />
          {saved && <span className="text-[12px] text-ok">Saved ✓</span>}
          <Button variant="primary" busy={busy === "save"} onClick={() => void save()}>
            Save
          </Button>
        </div>
      </div>
    </div>
  );
}

function StatusChip({ ok, okText, badText }: { ok?: boolean; okText: string; badText: string }) {
  if (ok === undefined) return null;
  return (
    <span className={`inline-flex items-center gap-1 rounded px-1.5 py-0.5 text-[10px] ${
      ok ? "bg-ok/15 text-ok" : "bg-bad/15 text-bad"
    }`}>
      {ok ? <CheckCircle2 size={9} /> : <XCircle size={9} />}{ok ? okText : badText}
    </span>
  );
}
