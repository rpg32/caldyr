"""Generate docs/img/distillation-reflux-tradeoff.png — the reflux-ratio
capital-vs-utility trade-off chart used in the distillation tutorial.

This one docs image is a *data plot* (matplotlib), not a UI screenshot: it sweeps
the shortcut benzene/toluene column's reflux factor, re-costs each design through
the engine's TEA pipeline, and plots capital (TCI) and annual utilities vs reflux.

    python web/scripts/reflux_chart.py

Requires the caldyr engine importable (matplotlib too). On this dev box the engine
lives in Python 3.11 — run with C:\\Python311\\python.exe. See README.md here.
"""
import sys
from pathlib import Path

# repo root is two levels up from web/scripts/
ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "engine"))
OUT = ROOT / "docs" / "img" / "distillation-reflux-tradeoff.png"

import matplotlib  # noqa: E402
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from caldyr.core import Component, Flowsheet  # noqa: E402
from caldyr.economics import TEAConfig, analyze  # noqa: E402
from caldyr.unitops import ShortcutColumn  # noqa: E402

P_ATM = 101325.0


def build(rr: float) -> Flowsheet:
    fs = Flowsheet(components=[Component("benzene"), Component("toluene")],
                   property_package="thermo:PR")
    fs.add(ShortcutColumn("COL", {
        "light_key": "benzene", "heavy_key": "toluene",
        "recovery_light": 0.99, "recovery_heavy": 0.98,
        "rr_factor": rr, "P": P_ATM,
    }))
    fs.feed("FEED", "COL:in1", T=365.0, P=P_ATM, molar_flow=100.0,
            z={"benzene": 0.5, "toluene": 0.5})
    for sid, port in (("DIST", "distillate"), ("BOT", "bottoms"),
                      ("QC", "condenser_duty"), ("QR", "reboiler_duty")):
        fs.connect(sid, f"COL:{port}", None)
    return fs


def main() -> None:
    rrs = [1.05, 1.10, 1.20, 1.30, 1.50, 1.75, 2.00]
    tci, util = [], []
    for rr in rrs:
        fs = build(rr)
        rep = fs.solve()
        res = analyze(fs, rep, TEAConfig(product_component="benzene",
                                         product_min_fraction=0.9))
        tci.append(res.capital.tci / 1e6)
        util.append(res.opex.utilities / 1e6)
        print(f"rr={rr:>4}  TCI=${tci[-1]:.3f}M  utilities=${util[-1]:.3f}M/yr")

    plt.rcParams.update({"font.size": 11, "figure.dpi": 150})
    fig, ax1 = plt.subplots(figsize=(6.6, 4.0))
    c1, c2 = "#2563eb", "#ea580c"
    ax1.plot(rrs, tci, "o-", color=c1, lw=2, label="Capital (TCI)")
    ax1.set_xlabel("reflux ratio factor  (R / R$_{min}$)")
    ax1.set_ylabel("Total capital investment / \\$M", color=c1)
    ax1.tick_params(axis="y", labelcolor=c1)
    imin = tci.index(min(tci))
    ax1.annotate(f"capital min\n≈{rrs[imin]:.2f}× R$_{{min}}$",
                 xy=(rrs[imin], tci[imin]), xytext=(rrs[imin] + 0.18, tci[imin] + 0.03),
                 color=c1, fontsize=9, arrowprops=dict(arrowstyle="->", color=c1, lw=1.2))
    ax2 = ax1.twinx()
    ax2.plot(rrs, util, "s--", color=c2, lw=2, label="Utilities (opex)")
    ax2.set_ylabel("Annual utilities / \\$M·yr$^{-1}$", color=c2)
    ax2.tick_params(axis="y", labelcolor=c2)
    ax1.set_title("Benzene/toluene column: capital vs utilities vs reflux", fontsize=12)
    l1, lb1 = ax1.get_legend_handles_labels()
    l2, lb2 = ax2.get_legend_handles_labels()
    ax1.legend(l1 + l2, lb1 + lb2, loc="upper center", fontsize=9, frameon=False)
    ax1.grid(True, alpha=0.25)
    fig.tight_layout()
    fig.savefig(OUT, bbox_inches="tight")
    print("saved", OUT)


if __name__ == "__main__":
    main()
