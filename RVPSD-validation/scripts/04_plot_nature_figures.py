#!/usr/bin/env python3
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import pearsonr
from pathlib import Path

plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
    "font.size": 7,
    "axes.labelsize": 8,
    "axes.titlesize": 9,
    "xtick.labelsize": 7,
    "ytick.labelsize": 7,
    "legend.fontsize": 7,
    "figure.dpi": 300,
    "savefig.dpi": 300,
    "savefig.bbox": "tight",
    "savefig.pad_inches": 0.02,
    "axes.linewidth": 0.6,
    "xtick.major.width": 0.6,
    "ytick.major.width": 0.6,
    "xtick.minor.width": 0.4,
    "ytick.minor.width": 0.4,
    "xtick.major.size": 3,
    "ytick.major.size": 3,
    "xtick.minor.size": 2,
    "ytick.minor.size": 2,
    "axes.spines.top": False,
    "axes.spines.right": False,
})

NATURE_RED = "#D55E00"
NATURE_ORANGE = "#E69F00"
NATURE_BLUE = "#56B4E9"
NATURE_GREEN = "#009E73"
NATURE_GRAY = "#999999"
NATURE_DARK = "#333333"

OUTDIR = Path("/lustre/home/acct-bioxwz/bioxwz-yqz/comparison_work/results/openclaw")
CSV_PATH = OUTDIR / "validation_with_plddt.csv"

def load_data():
    df = pd.read_csv(CSV_PATH)
    df = df.dropna(subset=["mean_plddt", "tm_score"])
    df["plddt_bin"] = pd.cut(df["mean_plddt"], bins=[0,50,70,90,100], labels=["<50","50–70","70–90","≥90"])
    return df

def format_pvalue(p):
    if p < 1e-300:
        return "< 1 × 10$^{-300}$"
    elif p < 1e-10:
        exp = int(np.floor(np.log10(p)))
        return f"< 1 × 10$^{{{exp}}}$"
    else:
        return f"= {p:.2e}"

def plot_scatter_nature(df):
    fig, ax = plt.subplots(figsize=(8.5/2.54, 8.5/2.54))
    hb = ax.hexbin(df["mean_plddt"], df["tm_score"], gridsize=45, cmap="Blues", mincnt=1, edgecolors="none")
    z = np.polyfit(df["mean_plddt"], df["tm_score"], 1)
    p_line = np.poly1d(z)
    x_line = np.linspace(20, 100, 200)
    ax.plot(x_line, p_line(x_line), color=NATURE_RED, linewidth=1.5, zorder=5)
    ax.axhline(0.8, color=NATURE_GREEN, linestyle="--", linewidth=0.8, alpha=0.7)
    ax.axhline(0.5, color=NATURE_ORANGE, linestyle="--", linewidth=0.8, alpha=0.7)
    ax.axvline(70, color=NATURE_GRAY, linestyle=":", linewidth=0.6, alpha=0.5)
    ax.axvline(90, color=NATURE_GRAY, linestyle=":", linewidth=0.6, alpha=0.5)
    r, pval = pearsonr(df["mean_plddt"], df["tm_score"])
    stats_text = f"$r$ = {r:.3f}\\n$p$ {format_pvalue(pval)}\\n$n$ = {len(df):,}"
    ax.text(0.97, 0.05, stats_text, transform=ax.transAxes, fontsize=7, va="bottom", ha="right",
            color=NATURE_DARK, bbox=dict(boxstyle="round,pad=0.3", facecolor="white", edgecolor=NATURE_GRAY, alpha=0.95, linewidth=0.5))
    ax.text(35, 0.92, "High\\n(TM ≥ 0.8)", fontsize=6, color=NATURE_GREEN, ha="center", va="center", alpha=0.8)
    ax.text(35, 0.55, "Same fold\\n(TM ≥ 0.5)", fontsize=6, color=NATURE_ORANGE, ha="center", va="center", alpha=0.8)
    ax.set_xlabel("Mean pLDDT", labelpad=3)
    ax.set_ylabel("TM-score", labelpad=3)
    ax.set_xlim(15, 102)
    ax.set_ylim(-0.02, 1.05)
    ax.set_xticks([20, 40, 60, 80, 100])
    ax.set_yticks([0.0, 0.2, 0.4, 0.6, 0.8, 1.0])
    cbar = plt.colorbar(hb, ax=ax, shrink=0.6, pad=0.02, aspect=15)
    cbar.set_label("Count", fontsize=6, labelpad=2)
    cbar.ax.tick_params(labelsize=6, width=0.5)
    plt.savefig(OUTDIR / "Fig1_pLDDT_TMscore_Nature.pdf", format="pdf")
    plt.savefig(OUTDIR / "Fig1_pLDDT_TMscore_Nature.png", format="png")
    plt.close()
    print("[OK] Fig1 saved")

def plot_stratified_nature(df):
    fig, ax = plt.subplots(figsize=(8.5/2.54, 7.5/2.54))
    palette = [NATURE_RED, NATURE_ORANGE, NATURE_BLUE, NATURE_GREEN]
    sns.boxplot(data=df, x="plddt_bin", y="tm_score", palette=palette, width=0.6, linewidth=0.8, fliersize=2, whis=1.5, ax=ax)
    medians = df.groupby("plddt_bin", observed=False)["tm_score"].median()
    for i, med in enumerate(medians):
        ax.scatter(i, med, marker="D", s=25, color="white", edgecolor=NATURE_DARK, linewidth=0.8, zorder=5)
    ax.axhline(0.8, color=NATURE_GREEN, linestyle="--", linewidth=0.8, alpha=0.6)
    ax.axhline(0.5, color=NATURE_ORANGE, linestyle="--", linewidth=0.8, alpha=0.6)
    counts = df["plddt_bin"].value_counts().sort_index()
    for i, count in enumerate(counts):
        y_max = df[df["plddt_bin"] == counts.index[i]]["tm_score"].max()
        ax.text(i, min(y_max + 0.04, 1.02), f"$n$={count:,}", ha="center", va="bottom", fontsize=6, color=NATURE_DARK)
    ax.set_xlabel("pLDDT confidence bin", labelpad=3)
    ax.set_ylabel("TM-score", labelpad=3)
    ax.set_ylim(-0.05, 1.12)
    ax.set_yticks([0.0, 0.2, 0.4, 0.6, 0.8, 1.0])
    ax.yaxis.grid(True, linestyle="-", linewidth=0.3, color="#E0E0E0", zorder=0)
    ax.set_axisbelow(True)
    plt.savefig(OUTDIR / "Fig2_pLDDT_stratified_Nature.pdf", format="pdf")
    plt.savefig(OUTDIR / "Fig2_pLDDT_stratified_Nature.png", format="png")
    plt.close()
    print("[OK] Fig2 saved")

def plot_combined_panel(df):
    fig = plt.figure(figsize=(17.8/2.54, 8.0/2.54))
    gs = fig.add_gridspec(1, 2, wspace=0.35)
    r, pval = pearsonr(df["mean_plddt"], df["tm_score"])
    p_str = format_pvalue(pval)

    ax1 = fig.add_subplot(gs[0, 0])
    hb = ax1.hexbin(df["mean_plddt"], df["tm_score"], gridsize=40, cmap="Blues", mincnt=1, edgecolors="none")
    z = np.polyfit(df["mean_plddt"], df["tm_score"], 1)
    p_line = np.poly1d(z)
    x_line = np.linspace(20, 100, 200)
    ax1.plot(x_line, p_line(x_line), color=NATURE_RED, linewidth=1.2, zorder=5)
    ax1.axhline(0.8, color=NATURE_GREEN, linestyle="--", lw=0.7, alpha=0.6)
    ax1.axhline(0.5, color=NATURE_ORANGE, linestyle="--", lw=0.7, alpha=0.6)
    ax1.axvline(70, color=NATURE_GRAY, linestyle=":", lw=0.5, alpha=0.4)
    ax1.axvline(90, color=NATURE_GRAY, linestyle=":", lw=0.5, alpha=0.4)
    ax1.text(0.97, 0.05, f"$r$ = {r:.3f}\\n$p$ {p_str}\\n$n$ = {len(df):,}",
             transform=ax1.transAxes, fontsize=6.5, va="bottom", ha="right",
             bbox=dict(boxstyle="round,pad=0.25", facecolor="white", edgecolor=NATURE_GRAY, alpha=0.95, lw=0.5))
    ax1.set_xlabel("Mean pLDDT", labelpad=2)
    ax1.set_ylabel("TM-score", labelpad=2)
    ax1.set_xlim(15, 102)
    ax1.set_ylim(-0.02, 1.05)
    ax1.set_xticks([20, 40, 60, 80, 100])
    ax1.set_yticks([0.0, 0.2, 0.4, 0.6, 0.8, 1.0])
    cbar = plt.colorbar(hb, ax=ax1, shrink=0.5, pad=0.02, aspect=12)
    cbar.set_label("Count", fontsize=5.5, labelpad=1)
    cbar.ax.tick_params(labelsize=5.5, width=0.4)

    ax2 = fig.add_subplot(gs[0, 1])
    palette = [NATURE_RED, NATURE_ORANGE, NATURE_BLUE, NATURE_GREEN]
    sns.boxplot(data=df, x="plddt_bin", y="tm_score", palette=palette, width=0.55, linewidth=0.7, fliersize=1.5, whis=1.5, ax=ax2)
    medians = df.groupby("plddt_bin", observed=False)["tm_score"].median()
    for i, med in enumerate(medians):
        ax2.scatter(i, med, marker="D", s=20, color="white", edgecolor=NATURE_DARK, linewidth=0.7, zorder=5)
    ax2.axhline(0.8, color=NATURE_GREEN, linestyle="--", lw=0.7, alpha=0.6)
    ax2.axhline(0.5, color=NATURE_ORANGE, linestyle="--", lw=0.7, alpha=0.6)
    counts = df["plddt_bin"].value_counts().sort_index()
    for i, count in enumerate(counts):
        y_max = df[df["plddt_bin"] == counts.index[i]]["tm_score"].max()
        ax2.text(i, min(y_max + 0.03, 1.02), f"$n$={count:,}", ha="center", va="bottom", fontsize=5.5)
    ax2.set_xlabel("pLDDT confidence bin", labelpad=2)
    ax2.set_ylabel("TM-score", labelpad=2)
    ax2.set_ylim(-0.05, 1.12)
    ax2.set_yticks([0.0, 0.2, 0.4, 0.6, 0.8, 1.0])
    ax2.yaxis.grid(True, linestyle="-", linewidth=0.25, color="#E0E0E0", zorder=0)
    ax2.set_axisbelow(True)

    ax1.text(-0.18, 1.05, "a", transform=ax1.transAxes, fontsize=10, fontweight="bold", va="top")
    ax2.text(-0.18, 1.05, "b", transform=ax2.transAxes, fontsize=10, fontweight="bold", va="top")

    plt.savefig(OUTDIR / "Fig3_combined_Nature.pdf", format="pdf")
    plt.savefig(OUTDIR / "Fig3_combined_Nature.png", format="png")
    plt.close()
    print("[OK] Fig3 saved")

def main():
    if not CSV_PATH.exists():
        print(f"[ERROR] {CSV_PATH} not found.")
        return
    df = load_data()
    print(f"[INFO] Loaded {len(df)} pairs.")
    plot_scatter_nature(df)
    plot_stratified_nature(df)
    plot_combined_panel(df)
    print("[ALL DONE]")

if __name__ == "__main__":
    main()