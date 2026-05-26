#!/usr/bin/env python3
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import pearsonr
from pathlib import Path

QUERY_DIR = Path("/lustre/home/acct-bioxsyy/share/yangqiangzhen/ruijin_RNA_virus_project/.../all_best_pdb_files")
VAL_CSV = Path("/lustre/home/acct-bioxwz/bioxwz-yqz/comparison_work/results/openclaw/high_similarity_validation.csv")
OUTDIR = Path("/lustre/home/acct-bioxwz/bioxwz-yqz/comparison_work/results/openclaw")
OUTDIR.mkdir(parents=True, exist_ok=True)

def extract_mean_plddt(pdb_path):
    plddts = []
    try:
        with open(pdb_path, "r") as f:
            for line in f:
                if line.startswith(("ATOM  ", "HETATM")):
                    bfac_str = line[60:66].strip()
                    if bfac_str:
                        try:
                            plddts.append(float(bfac_str))
                        except ValueError:
                            continue
    except Exception as e:
        print(f"  [WARN] Error reading {pdb_path}: {e}")
        return np.nan
    if not plddts:
        return np.nan
    return float(np.mean(plddts))

def find_query_file(query_id, query_dir):
    for cand in [query_dir / f"{query_id}.pdb", query_dir / query_id, query_dir / f"{query_id}.cif"]:
        if cand.exists() and cand.stat().st_size > 1000:
            return cand
    return None

def main():
    val_df = pd.read_csv(VAL_CSV)
    unique_queries = val_df["query_id"].unique()
    print(f"[INFO] Unique queries: {len(unique_queries)}")

    plddt_records = []
    for i, qid in enumerate(unique_queries, 1):
        qfile = find_query_file(qid, QUERY_DIR)
        if qfile is None:
            plddt_records.append({"query_id": qid, "mean_plddt": np.nan, "query_path": None})
            continue
        plddt = extract_mean_plddt(qfile)
        plddt_records.append({"query_id": qid, "mean_plddt": plddt, "query_path": str(qfile)})
        if i % 500 == 0:
            print(f"  [{i}/{len(unique_queries)}] {qid}: pLDDT={plddt:.2f}")

    plddt_df = pd.DataFrame(plddt_records)
    plddt_df.to_csv(OUTDIR / "plddt_per_structure.csv", index=False)

    merged = val_df.merge(plddt_df[["query_id", "mean_plddt"]], on="query_id", how="left")
    merged = merged.dropna(subset=["mean_plddt", "tm_score"])
    r, p = pearsonr(merged["mean_plddt"], merged["tm_score"])

    merged["plddt_bin"] = pd.cut(merged["mean_plddt"], bins=[0,50,70,90,100], labels=["<50","50-70","70-90","≥90"])
    stratified = merged.groupby("plddt_bin", observed=False)["tm_score"].agg(["count","mean","median","std"])

    print("\n========== CORRELATION RESULTS ==========")
    print(f"N (valid pairs)       : {len(merged)}")
    print(f"Pearson r             : {r:.4f}")
    print(f"P-value               : {p:.2e}")
    print(f"R-squared             : {r**2:.4f}")
    print("\n----- Stratified by pLDDT -----")
    print(stratified.to_string())

    # Plot
    fig, ax = plt.subplots(figsize=(6,5))
    sns.regplot(data=merged, x="mean_plddt", y="tm_score",
                scatter_kws={"alpha":0.4,"s":20,"color":"steelblue"},
                line_kws={"color":"crimson"}, ax=ax)
    ax.axhline(0.8, color="green", ls="--", lw=1, label="High quality (TM=0.8)")
    ax.axhline(0.5, color="orange", ls="--", lw=1, label="Same fold (TM=0.5)")
    ax.set_xlabel("Mean pLDDT", fontsize=12)
    ax.set_ylabel("TM-score (vs. Experimental PDB)", fontsize=12)
    ax.set_title(f"N={len(merged)} | r={r:.3f} | p={p:.2e}", fontsize=12)
    ax.legend(loc="lower right")
    plt.tight_layout()
    plt.savefig(OUTDIR / "Fig_pLDDT_vs_TMscore.pdf", dpi=300)
    plt.savefig(OUTDIR / "Fig_pLDDT_vs_TMscore.png", dpi=300)

    merged.to_csv(OUTDIR / "validation_with_plddt.csv", index=False)
    print(f"\n[INFO] Saved: validation_with_plddt.csv")

if __name__ == "__main__":
    main()