#!/usr/bin/env python3
import os, re, argparse, subprocess, pandas as pd
from pathlib import Path
from urllib.request import urlretrieve
from multiprocessing import Pool

def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--m8", required=True)
    p.add_argument("--query_dir", required=True)
    p.add_argument("--outdir", required=True)
    p.add_argument("--usalign", required=True)
    p.add_argument("--prob_thresh", type=float, default=0.5)
    p.add_argument("--max_hits", type=int, default=1)
    p.add_argument("--threads", type=int, default=1)
    return p.parse_args()

def extract_pdb_id(target_full):
    t = target_full.lower()
    m = re.match(r'([0-9][a-z0-9]{3})', t)
    return m.group(1) if m else t.split('_')[0].split('-')[0].split('.')[0]

def get_local_pdb_path(pdb_id, pdb_dl_dir):
    pdb_id = pdb_id.lower()
    for ext in ['pdb', 'cif']:
        f = pdb_dl_dir / f"{pdb_id}.{ext}"
        if f.exists() and f.stat().st_size > 1000:
            return str(f)
    return None

def download_pdb(pdb_id, outdir):
    pdb_id = pdb_id.lower()
    outdir = Path(outdir)
    existing = get_local_pdb_path(pdb_id, outdir)
    if existing:
        return existing
    for ext in ['pdb', 'cif']:
        f = outdir / f"{pdb_id}.{ext}"
        url = f"https://files.rcsb.org/download/{pdb_id}.{ext}"
        try:
            urlretrieve(url, f)
            if f.exists() and f.stat().st_size > 1000:
                return str(f)
            f.unlink(missing_ok=True)
        except Exception:
            f.unlink(missing_ok=True)
    return None

def run_usalign(qpdb, tpdb, usalign_path):
    infmt1 = "1" if str(qpdb).lower().endswith(".cif") else "0"
    infmt2 = "1" if str(tpdb).lower().endswith(".cif") else "0"
    cmd = [usalign_path, str(qpdb), str(tpdb), "-ter", "0", "-infmt1", infmt1, "-infmt2", infmt2]
    try:
        res = subprocess.run(cmd, capture_output=True, text=True, timeout=180)
        out = res.stdout
        tm = re.search(r"TM-score=\s+([0-9.]+)", out)
        rmsd = re.search(r"RMSD=\s+([0-9.]+)", out)
        if tm and rmsd:
            return float(tm.group(1)), float(rmsd.group(1)), out
        return None, None, out[:800]
    except Exception as e:
        return None, None, str(e)

def process_pair(task):
    idx, row, query_dir, usalign_path, outdir, n_total = task
    query_dir = Path(query_dir)
    outdir = Path(outdir)
    logs = outdir / "usalign_logs"
    logs.mkdir(parents=True, exist_ok=True)

    query_name = str(row["query"])
    target_full = str(row["target"])
    pdb_id = extract_pdb_id(target_full)
    chain = ""
    if "_" in target_full:
        parts = target_full.rsplit("_", 1)
        if len(parts) == 2 and len(parts[1]) == 1 and parts[1].isalpha():
            chain = parts[1]

    query_pdb = None
    for cand in [query_dir / f"{query_name}.pdb", query_dir / query_name, query_dir / f"{query_name}.cif"]:
        if cand.exists():
            query_pdb = cand
            break
    if query_pdb is None:
        return None, f"[{idx+1}/{n_total}] SKIP query not found: {query_name}"

    target_pdb = get_local_pdb_path(pdb_id, outdir / "pdb_downloads")
    if target_pdb is None:
        return None, f"[{idx+1}/{n_total}] SKIP target not in local pool: {pdb_id}"

    safe_target = target_full.replace('/', '_').replace('\\', '_')
    logf = logs / f"{query_name}_vs_{safe_target}.txt"
    if logf.exists() and logf.stat().st_size > 100:
        with open(logf) as f:
            content = f.read()
        tm = re.search(r"TM-score=\s+([0-9.]+)", content)
        rmsd = re.search(r"RMSD=\s+([0-9.]+)", content)
        if tm and rmsd:
            record = {
                "query_id": query_name, "query_file": str(query_pdb),
                "pdb100_target_id": target_full, "pdb_id": pdb_id, "chain": chain,
                "experimental_pdb_file": target_pdb,
                "foldseek_prob": round(float(row["prob"]), 4),
                "foldseek_fident": round(float(row["fident"]), 4),
                "foldseek_evalue": row["evalue"], "foldseek_alnlen": int(row["alnlen"]),
                "query_start": int(row["qstart"]), "query_end": int(row["qend"]),
                "target_start": int(row["tstart"]), "target_end": int(row["tend"]),
                "tm_score": round(float(tm.group(1)), 4),
                "rmsd_A": round(float(rmsd.group(1)), 3),
                "foldseek_query_aln_seq": str(row["qaln"]),
                "foldseek_target_aln_seq": str(row["taln"])
            }
            return record, f"[{idx+1}/{n_total}] CACHED {query_name} vs {target_full}: TM={float(tm.group(1)):.3f}"

    tm_score, rmsd, usalign_out = run_usalign(query_pdb, target_pdb, usalign_path)
    if tm_score is None:
        fail_log = logs / f"FAIL_{query_name}_vs_{safe_target}.txt"
        with open(fail_log, "w") as f:
            f.write(usalign_out)
        return None, f"[{idx+1}/{n_total}] FAIL USalign: {query_name} vs {target_full}"

    with open(logf, "w") as f:
        f.write(usalign_out)

    record = {
        "query_id": query_name, "query_file": str(query_pdb),
        "pdb100_target_id": target_full, "pdb_id": pdb_id, "chain": chain,
        "experimental_pdb_file": target_pdb,
        "foldseek_prob": round(float(row["prob"]), 4),
        "foldseek_fident": round(float(row["fident"]), 4),
        "foldseek_evalue": row["evalue"], "foldseek_alnlen": int(row["alnlen"]),
        "query_start": int(row["qstart"]), "query_end": int(row["qend"]),
        "target_start": int(row["tstart"]), "target_end": int(row["tend"]),
        "tm_score": round(tm_score, 4), "rmsd_A": round(rmsd, 3),
        "foldseek_query_aln_seq": str(row["qaln"]),
        "foldseek_target_aln_seq": str(row["taln"])
    }
    return record, f"[{idx+1}/{n_total}] OK {query_name} vs {target_full}: TM={tm_score:.3f} RMSD={rmsd:.2f}A"

def main():
    args = parse_args()
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    pdb_dl = outdir / "pdb_downloads"
    pdb_dl.mkdir(exist_ok=True)

    cols = ["query","target","fident","alnlen","mismatch","gapopen",
            "qstart","qend","tstart","tend","evalue","bits","prob",
            "qlen","tlen","qaln","taln"]
    df = pd.read_csv(args.m8, sep="\t", names=cols, comment="#", low_memory=False)
    df = df[df["prob"] >= args.prob_thresh].copy()
    df = df.sort_values(["query","prob"], ascending=[True, False])
    selected = df.groupby("query").head(args.max_hits).reset_index(drop=True)
    n_total = len(selected)
    print(f"[INFO] Pairs to validate: {n_total}")

    unique_pdbs = selected["target"].apply(extract_pdb_id).unique()
    already_local = sum(1 for pid in unique_pdbs if get_local_pdb_path(pid, pdb_dl))
    to_download = [pid for pid in unique_pdbs if not get_local_pdb_path(pid, pdb_dl)]
    print(f"[INFO] Already local: {already_local} / {len(unique_pdbs)}")
    print(f"[INFO] Need download: {len(to_download)}")

    if to_download:
        print(f"[INFO] Downloading {len(to_download)} missing PDBs ...")
        for i, pid in enumerate(to_download, 1):
            download_pdb(pid, pdb_dl)
            if i % 50 == 0:
                print(f"  ... {i}/{len(to_download)}")

    tasks = [(idx, row, args.query_dir, args.usalign, args.outdir, n_total)
             for idx, row in selected.iterrows()]
    records = []
    if args.threads > 1:
        with Pool(args.threads) as pool:
            for rec, msg in pool.imap_unordered(process_pair, tasks):
                print(msg, flush=True)
                if rec: records.append(rec)
    else:
        for task in tasks:
            rec, msg = process_pair(task)
            print(msg, flush=True)
            if rec: records.append(rec)

    outcsv = outdir / "high_similarity_validation.csv"
    if records:
        pd.DataFrame(records).to_csv(outcsv, index=False)
        print(f"\n[SUCCESS] {len(records)} alignments saved to {outcsv}")
    else:
        print("[WARN] No successful alignments.")

if __name__ == "__main__":
    main()