#!/usr/bin/env python3
"""Paper figures from results/results.json. One y-scale per axes, labelled axes throughout."""
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
R = json.loads((ROOT / "results" / "results.json").read_text())
S = json.loads((ROOT / "results" / "real_solves.json").read_text())
FIG = ROOT / "results" / "figures"
BLUE, ORANGE, AQUA, YELLOW, MAGENTA = "#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4"
INK, INK2, GRID, SURFACE = "#0b0b0b", "#52514e", "#e6e5e0", "#ffffff"
plt.rcParams.update({
    "figure.dpi": 100, "savefig.dpi": 220, "font.size": 10, "axes.edgecolor": INK2, "axes.labelcolor": INK, "text.color": INK,
    "xtick.color": INK2, "ytick.color": INK2, "axes.grid": True, "grid.color": GRID, "grid.linewidth": 0.8, "axes.axisbelow": True,
    "axes.spines.top": False, "axes.spines.right": False, "figure.facecolor": SURFACE, "axes.facecolor": SURFACE, "legend.frameon": False,
    "lines.linewidth": 2,
})
T = R["techniques"]


def pad(p, n):
    a = np.zeros(n)
    a[: len(p)] = p
    return a


def save(fig, name):
    fig.tight_layout()
    fig.savefig(FIG / f"{name}.png")
    plt.close(fig)


# 1 ---------- the three independent blocks
fig, ax = plt.subplots(1, 3, figsize=(10, 2.9), sharey=True)
ls, oll, pll = R["block_pmfs"]["ls"], R["block_pmfs"]["oll"], R["block_pmfs"]["pll"]
for a, (title, p) in zip(ax, [("Last pair: AUF + F2L alg", ls), ("OLL: AUF + alg", oll), ("PLL: AUF + alg + AUF", pll)]):
    a.bar(np.arange(len(p)), p, color=BLUE, width=0.8)
    a.set_title(title, fontsize=10, loc="left")
    a.set_xlabel("Moves (STM)")
    m = np.arange(len(p)) @ np.array(p)
    a.axvline(m, color=INK2, lw=1, ls="--")
    a.set_ylim(0, 0.62)
    a.text(m, 0.585, f"mean {m:.2f}", color=INK2, fontsize=9, ha="center", backgroundcolor=SURFACE)
    a.xaxis.set_major_locator(matplotlib.ticker.MaxNLocator(integer=True))
ax[0].set_ylabel("Probability")
save(fig, "fig1_blocks")

# 2 ---------- enumeration vs flat AUF convolution
b = np.array(R["baseline"]["pmf"])
f = np.array(R["baseline"]["convolution_vs_enumeration"]["flat_pmf"])
n = max(len(b), len(f))
x = np.arange(n)
fig, ax = plt.subplots(1, 2, figsize=(10, 3.2))
ax[0].bar(x, pad(b, n), color=BLUE, width=0.8, label="Full enumeration")
ax[0].plot(x, pad(f, n), color=ORANGE, marker="o", ms=3.5, label="Convolution with flat AUFs (0.25, 0.75)")
ax[0].set_xlim(8, 52); ax[0].set_xlabel("Total moves T (STM)"); ax[0].set_ylabel("P(T = m)"); ax[0].set_ylim(0, 0.125); ax[0].legend(loc="upper left", fontsize=8.5)
ax[1].bar(x, (pad(f, n) - pad(b, n)) * 100, color=ORANGE, width=0.8)
ax[1].set_xlim(8, 52); ax[1].set_xlabel("Total moves T (STM)"); ax[1].set_ylabel("Flat-AUF error (percentage points)")
ax[1].axhline(0, color=INK2, lw=1)
save(fig, "fig2_flat_auf_error")

# 3 ---------- wraparound
fig, ax = plt.subplots(figsize=(7, 3.1))
W = R["wraparound"]["tests"]
n_small = "32"
w = np.array(W[n_small]["pmf"])
ax.bar(np.arange(len(b)), b, color=BLUE, width=0.8, label=f"True distribution (N = {R['wraparound']['required_N']})")
ax.bar(np.arange(len(w)), w, color=ORANGE, width=0.45, label="DFT with N = 32: mass above 31 wraps to the start")
ax.set_xlabel("Index read as total moves (STM)"); ax.set_ylabel("Probability"); ax.legend(loc="upper left", fontsize=8.5)
ax.set_ylim(0, max(b.max(), w.max()) * 1.35)
save(fig, "fig3_wraparound")

# 4, 5 ------- overlay of pmf and cdf
series = [("Baseline", BLUE, "-"), ("All CFOP techniques, fewest moves", ORANGE, "--"), ("ZBLS + ZBLL", AQUA, "-."), ("1LLL", YELLOW, ":")]
labels = {"Baseline": "Baseline (119 algs)", "All CFOP techniques, fewest moves": "All CFOP techniques, fewest moves", "ZBLS + ZBLL": "ZBLS + ZBLL", "1LLL": "1LLL (one alg for the last layer)"}
for cdf, name, ylabel in ((False, "fig4_overlay_pmf", "P(T = m)"), (True, "fig5_overlay_cdf", "P(T ≤ m)")):
    fig, ax = plt.subplots(figsize=(7.2, 3.5))
    for k, c, st in series:
        p = np.array(T[k]["pmf"])
        y = np.cumsum(p) if cdf else p
        ax.plot(np.arange(len(p)), y, color=c, ls=st, label=f"{labels[k]}: mean {T[k]['mean']:.2f}")
    ax.set_xlim(8, 50); ax.set_xlabel("Total moves m (STM)"); ax.set_ylabel(ylabel)
    ax.set_ylim(0, 1.02 if cdf else 0.19)
    if cdf:
        ax.legend(fontsize=8.5, loc="upper center", bbox_to_anchor=(0.5, -0.2), ncol=2)
    else:
        ax.legend(fontsize=8.5, loc="upper right")
    save(fig, name)

# 6 ---------- moves saved per technique
names = ["U2 insert, always", "Sledge insert, always", "U2 insert, two-step", "Sledge insert, two-step", "U2 + sledge, two-step", "U2 + sledge, full look-ahead", "Cancellations",
         "Angles OLL, speed-gated", "Angles PLL, speed-gated", "Angles both, speed-gated", "Angles both, fewest moves", "Angles OLL, two-step", "COLL, always when available", "COLL, two-step",
         "All CFOP techniques, speed-gated", "All CFOP techniques, fewest moves"]
vals = [T[k]["saved"] for k in names]
fig, ax = plt.subplots(figsize=(7.2, 4.6))
y = np.arange(len(names))[::-1]
ax.barh(y, vals, color=[BLUE if v >= 0 else ORANGE for v in vals], height=0.7)
ax.set_yticks(y); ax.set_yticklabels(names, fontsize=8.5); ax.axvline(0, color=INK2, lw=1)
for yi, v in zip(y, vals):
    ax.text(v + (0.03 if v >= 0 else -0.03), yi, f"{v:+.2f}", va="center", ha="left" if v >= 0 else "right", fontsize=8, color=INK)
ax.set_xlim(min(vals) - 0.35, max(vals) + 0.35)
ax.set_xlabel("Mean moves saved against the baseline (STM)"); ax.grid(axis="y", visible=False)
save(fig, "fig6_moves_saved")

# 7 ---------- trade-off
pts = [("U2 + sledge, two-step", "", (0, 0)), ("Cancellations", "Cancellations; inserts", (7, 5)), ("COLL, two-step", "COLL", (-30, 2)), ("OLLCP, two-step", "OLLCP", (7, 3)), ("ZBLL only (no edge control)", "ZBLL only", (7, 3)), ("ZBLS + ZBLL", "ZBLS + ZBLL", (-30, 8)), ("OLS + PLL", "OLS + PLL", (-20, 9)), ("1LLL", "1LLL", (7, 3)), ("Angles both, fewest moves", "Angles (fewest moves)", (-4, 8))]
fig, ax = plt.subplots(figsize=(7.2, 3.6))
for k, lab, off in pts:
    added = max(T[k]["algs_added"], 0)
    xv = added if added else 1
    ax.plot(xv, T[k]["saved"], "o", color=BLUE, ms=8, mec=SURFACE, mew=2)
    ax.annotate(lab, (xv, T[k]["saved"]), textcoords="offset points", xytext=off, fontsize=8.5, color=INK)
ax.set_xscale("log"); ax.set_xlim(0.8, 25000); ax.set_ylim(0, 9.8)
ax.set_xticks([1, 10, 100, 1000, 10000]); ax.set_xticklabels(["0", "10", "100", "1,000", "10,000"])
ax.set_xlabel("Algorithms added to the 119 baseline algs (log scale; 0 shown at the left edge)"); ax.set_ylabel("Mean moves saved (STM)")
save(fig, "fig7_tradeoff")

# 8 ---------- paired: moves saved state by state
fig, ax = plt.subplots(figsize=(7.2, 3.2))
for k, c, off, lab in (("U2 + sledge, two-step", BLUE, -0.2, "Two-step rule"), ("U2 + sledge, full look-ahead", ORANGE, 0.2, "Full look-ahead")):
    d = {int(a): b for a, b in T[k]["saved_pmf"].items() if a != "0"}
    xs = np.array(sorted(d))
    ax.bar(xs + off, [d[v] for v in xs], width=0.38, color=c, label=lab)
ax.set_xlabel("Moves saved against the baseline in the same state (0 omitted; negative: the choice cost moves)", fontsize=9); ax.set_ylabel("Probability")
ax.legend(); ax.set_xlim(-14, 18)
save(fig, "fig8_paired_saving")

# 9 ---------- naive convolution when independence fails
kname = R["independence_failure_example"]["technique"]
p = np.array(T[kname]["pmf"])
fig, ax = plt.subplots(figsize=(7.2, 3.2))
q = np.array(R["independence_failure_example"]["conv_pmf"])
n = max(len(p), len(q))
ax.bar(np.arange(n), pad(p, n), color=BLUE, width=0.8, label="Full enumeration")
ax.plot(np.arange(n), pad(q, n), color=ORANGE, marker="o", ms=3.5, label="Convolution of its step distributions")
ax.set_xlim(8, 50); ax.set_xlabel("Total moves T (STM)"); ax.set_ylabel("P(T = m)"); ax.legend(fontsize=8.5)
save(fig, "fig9_independence_failure")

# 10 --------- model vs reconstructions, last pair
m = np.array(S["last_pair_stm"]["model_pmf"])
r_all = np.array(S["last_pair_stm"]["real_cfop_all"]["pmf"])
fig, ax = plt.subplots(figsize=(7.2, 3.1))
x = np.arange(20)
ax.bar(x - 0.2, pad(m, 20), width=0.38, color=BLUE, label=f"Model: mean {S['last_pair_stm']['model_mean_including_auf']:.2f}")
ax.bar(x + 0.2, pad(r_all, 20), width=0.38, color=ORANGE, label=f"Pro solves, n = {S['last_pair_stm']['real_cfop_all']['n']}: mean {S['last_pair_stm']['real_cfop_all']['mean']:.2f}")
ax.set_xlabel("Moves in the last pair, AUF included (STM)"); ax.set_ylabel("Probability"); ax.set_xlim(-0.6, 16); ax.legend()
save(fig, "fig10_last_pair_vs_real")
print("figures written")
