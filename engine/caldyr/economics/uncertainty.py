"""Uncertainty: Monte-Carlo bands (P10/P50/P90) and a one-at-a-time tornado.

Both perturb the cheap financial layer (`evaluate_economics`) under input
uncertainty without re-solving the flowsheet — equipment sizes are fixed; what
moves is correlation error, prices, discount rate, and capacity factor.
"""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from . import data
from .analyze import TEAConfig, evaluate_economics

# Default uncertainty ranges (multiplicative CV for normals; (lo, hi) for uniforms).
DEFAULTS = dict(
    capex_cv=0.30,            # ±30% bare-module correlation error (Turton ~ class 5)
    product_price_cv=0.20,
    feed_price_cv=0.15,
    discount_rate=(0.08, 0.14),
    capacity_factor=(0.85, 0.98),
)


def _feed_components(fs) -> set[str]:
    """Components appearing in the flowsheet's boundary feed streams (no upstream
    unit) — the raw materials whose price the feed-price sensitivity perturbs.
    Derived from the flowsheet (not hard-coded), so the sweep is meaningful for
    any plant; for a separation/dehydration plant the product component is itself
    a feed, so its raw-material cost is (correctly) perturbed for the LCOP."""
    comps: set[str] = set()
    for conn in fs.connections:
        if conn.from_unit is None and conn.to_unit is not None:
            s = fs.streams.get(conn.stream_id)
            if s is not None and s.molar_flow:
                comps.update(s.normalized_z().keys())
    return comps


@dataclass
class MonteCarloResult:
    n: int
    lcop: dict = field(default_factory=dict)     # p10/p50/p90/mean/std
    npv: dict = field(default_factory=dict)
    lcop_samples: np.ndarray | None = field(default=None, repr=False)
    npv_samples: np.ndarray | None = field(default=None, repr=False)


def _stats(samples: np.ndarray) -> dict:
    return {
        "p10": float(np.percentile(samples, 10)),
        "p50": float(np.percentile(samples, 50)),
        "p90": float(np.percentile(samples, 90)),
        "mean": float(np.mean(samples)),
        "std": float(np.std(samples)),
    }


def _eval_lcop_npv(fs, sizes, cfg, *, capex_mult, prices, hours, rate):
    res = evaluate_economics(fs, sizes, cfg, capex_multiplier=capex_mult,
                             prices_per_kg=prices, operating_hours=hours, discount_rate=rate)
    return res.profitability.lcop, res.profitability.npv


def monte_carlo(fs, sizes, cfg: TEAConfig, *, n: int = 2000, seed: int = 0,
                **overrides) -> MonteCarloResult:
    p = {**DEFAULTS, **overrides}
    rng = np.random.default_rng(seed)
    base_prices = {**data.PRICES_PER_KG, **(cfg.prices_per_kg or {})}
    feeds = _feed_components(fs)

    lcops = np.empty(n)
    npvs = np.empty(n)
    for i in range(n):
        capex_mult = float(np.exp(rng.normal(0.0, p["capex_cv"])))     # lognormal ~ 1
        prod_mult = max(0.05, rng.normal(1.0, p["product_price_cv"]))
        feed_mult = max(0.05, rng.normal(1.0, p["feed_price_cv"]))
        rate = float(rng.uniform(*p["discount_rate"]))
        cap_factor = float(rng.uniform(*p["capacity_factor"]))

        prices = dict(base_prices)
        prices[cfg.product_component] = base_prices[cfg.product_component] * prod_mult
        for feed in feeds:
            if feed in prices:
                prices[feed] = base_prices[feed] * feed_mult

        lcops[i], npvs[i] = _eval_lcop_npv(
            fs, sizes, cfg, capex_mult=capex_mult, prices=prices,
            hours=8760.0 * cap_factor, rate=rate)

    return MonteCarloResult(n=n, lcop=_stats(lcops), npv=_stats(npvs),
                            lcop_samples=lcops, npv_samples=npvs)


@dataclass
class TornadoBar:
    variable: str
    low_value: float
    high_value: float
    low_lcop: float
    high_lcop: float

    @property
    def swing(self) -> float:
        return abs(self.high_lcop - self.low_lcop)


def tornado(fs, sizes, cfg: TEAConfig, **overrides) -> list[TornadoBar]:
    """One-at-a-time sensitivity of LCOP. Each variable is driven to its low and
    high while the others stay at base; bars are sorted by swing (largest first)."""
    p = {**DEFAULTS, **overrides}
    base_prices = {**data.PRICES_PER_KG, **(cfg.prices_per_kg or {})}
    base_hours = cfg.operating_hours
    base_rate = cfg.discount_rate
    prod = cfg.product_component
    feeds = _feed_components(fs)

    def lcop(*, capex_mult=1.0, prod_mult=1.0, feed_mult=1.0, hours=base_hours, rate=base_rate):
        prices = dict(base_prices)
        prices[prod] = base_prices[prod] * prod_mult
        for feed in feeds:
            if feed in prices:
                prices[feed] = base_prices[feed] * feed_mult
        return _eval_lcop_npv(fs, sizes, cfg, capex_mult=capex_mult, prices=prices,
                              hours=hours, rate=rate)[0]

    bars = [
        TornadoBar(f"capex +/-{p['capex_cv']:.0%}", 1 - p["capex_cv"], 1 + p["capex_cv"],
                   lcop(capex_mult=1 - p["capex_cv"]), lcop(capex_mult=1 + p["capex_cv"])),
        TornadoBar(f"feed price +/-{p['feed_price_cv']:.0%}",
                   1 - p["feed_price_cv"], 1 + p["feed_price_cv"],
                   lcop(feed_mult=1 - p["feed_price_cv"]),
                   lcop(feed_mult=1 + p["feed_price_cv"])),
        TornadoBar("discount rate", p["discount_rate"][0], p["discount_rate"][1],
                   lcop(rate=p["discount_rate"][0]), lcop(rate=p["discount_rate"][1])),
        TornadoBar("capacity factor", p["capacity_factor"][0], p["capacity_factor"][1],
                   lcop(hours=8760.0 * p["capacity_factor"][0]),
                   lcop(hours=8760.0 * p["capacity_factor"][1])),
    ]
    return sorted(bars, key=lambda b: b.swing, reverse=True)
