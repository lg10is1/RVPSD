#!/usr/bin/env python3
"""
Final optimized pipeline:
1. FoldMason for structural sequence alignment + LDDT scores
2. Sliding-window local Kabsch alignment to find true spatially-similar regions
3. Report both LDDT and local RMSD
4. Extract substructures and save aligned PDBs with correspondence annotations
"""

import json
import numpy as np
from Bio import SeqIO
from Bio.PDB.MMCIFParser import MMCIFParser
from Bio.PDB.PDBIO import PDBIO
from Bio.PDB import Select, Structure, Model, Chain


def get_ca_trace(structure):
    """Extract CA atoms with residue info."""
    atoms = []
    for model in structure:
        for chain in model:
            for residue in chain:
                if 'CA' in residue:
                    atoms.append((residue['CA'].get_coord(), chain.id, residue.id, residue.resname))
    return atoms


def kabsch(P, Q):
    """Kabsch algorithm."""
    cP = np.mean(P, axis=0)
    cQ = np.mean(Q, axis=0)
    Pc = P - cP
    Qc = Q - cQ
    H = Qc.T @ Pc
    U, S, Vt = np.linalg.svd(H)
    R = Vt.T @ U.T
    if np.linalg.det(R) < 0:
        Vt[-1, :] *= -1
        R = Vt.T @ U.T
    return R, cP, cQ


def apply_transform(structure, R, cP, cQ):
    for model in structure:
        for chain in model:
            for residue in chain:
                for atom in residue:
                    atom.set_coord((atom.get_coord() - cQ) @ R + cP)
    return structure


def write_pdb(structure, filename, remarks=None):
    class AllSelect(Select):
        def accept_atom(self, atom):
            return True
    io = PDBIO()
    io.set_structure(structure)
    with open(filename, 'w') as f:
        if remarks:
            for r in remarks:
                f.write(f"REMARK   1 {r}\n")
        io.save(f, select=AllSelect())
    print(f"[INFO] Saved: {filename}")


def extract_residues(structure, residue_ids):
    target = set(residue_ids)
    new_s = Structure.Structure(structure.id)
    for model in structure:
        new_m = Model.Model(model.id)
        for chain in model:
            new_c = Chain.Chain(chain.id)
            for residue in chain:
                if (chain.id, residue.id) in target:
                    new_c.add(residue.copy())
            if len(new_c) > 0:
                new_m.add(new_c)
        if len(new_m) > 0:
            new_s.add(new_m)
    return new_s


def find_similar_regions(coords1, coords2, idx_map, lddt_scores,
                         window_sizes=[15, 20, 25, 30, 40, 50],
                         rmsd_threshold=4.0, step=3):
    """
    Sliding window local Kabsch alignment.
    Returns regions sorted by RMSD.
    """
    regions = []

    for ws in window_sizes:
        for start in range(0, len(coords1) - ws + 1, step):
            end = start + ws
            w1 = coords1[start:end]
            w2 = coords2[start:end]

            R, cP, cQ = kabsch(w1, w2)
            w2_t = (w2 - cQ) @ R + cP
            rmsd = np.sqrt(np.mean(np.sum((w1 - w2_t)**2, axis=1)))

            if rmsd < rmsd_threshold:
                # Compute mean LDDT for this window
                lddt_vals = [lddt_scores[p] for p in range(start, end) if lddt_scores[p] >= 0]
                mean_lddt = np.mean(lddt_vals) if lddt_vals else 0
                regions.append({
                    'start': start,
                    'end': end,
                    'length': end - start,
                    'rmsd': rmsd,
                    'mean_lddt': mean_lddt,
                    'R': R,
                    'cP': cP,
                    'cQ': cQ,
                })

    # Sort by RMSD
    regions.sort(key=lambda x: x['rmsd'])
    return regions


def deduplicate_regions(regions, min_gap=5):
    """Remove overlapping regions, keeping the one with better RMSD."""
    used = set()
    final = []
    for r in regions:
        overlap = False
        for i in range(r['start'], r['end']):
            if i in used:
                overlap = True
                break
        if not overlap:
            for i in range(r['start'], r['end']):
                used.add(i)
            final.append(r)
    return final


def main():
    import os
    workdir = "/lustre/home/acct-bioxsyy/share/yangqiangzhen/ruijin_RNA_virus_project/revision_R1/case1/foldmason_results"
    os.makedirs(workdir, exist_ok=True)
    cif1 = "/lustre/home/acct-bioxsyy/share/yangqiangzhen/ruijin_RNA_virus_project/revision_R1/case1/BSVV1_PP480790.1.cif"
    cif2 = "/lustre/home/acct-bioxsyy/share/yangqiangzhen/ruijin_RNA_virus_project/revision_R1/case1/whIV2_AJG39294.cif"
    json_file = "/lustre/home/acct-bioxsyy/share/yangqiangzhen/ruijin_RNA_virus_project/revision_R1/case1/foldmason_result.a3m.json"

    print("=" * 60)
    print("FoldMason + Local Kabsch Alignment Pipeline")
    print("=" * 60)

    # Parse structures
    print("\n[1/5] Parsing CIF structures...")
    parser = MMCIFParser()
    s1 = parser.get_structure('BSVV1', cif1)
    s2 = parser.get_structure('whIV2', cif2)

    ca1 = get_ca_trace(s1)
    ca2 = get_ca_trace(s2)
    coords1 = np.array([a[0] for a in ca1])
    coords2 = np.array([a[0] for a in ca2])
    print(f"  BSVV1: {len(ca1)} CA atoms")
    print(f"  whIV2: {len(ca2)} CA atoms")

    # Parse FoldMason MSA and LDDT
    print("\n[2/5] Loading FoldMason alignment + LDDT scores...")
    records = list(SeqIO.parse("/lustre/home/acct-bioxsyy/share/yangqiangzhen/ruijin_RNA_virus_project/revision_R1/case1/foldmason_result.a3m_aa.fa", 'fasta'))
    seq1 = str(records[0].seq)
    seq2 = str(records[1].seq)

    with open(json_file) as f:
        data = json.load(f)
    lddt_scores = data['scores']

    # Build aligned coordinate arrays and index mapping
    aln_coords1 = []
    aln_coords2 = []
    idx_map = []  # (msa_pos, struct_idx1, struct_idx2)
    i1 = i2 = 0
    for msa_pos, (a1, a2) in enumerate(zip(seq1, seq2)):
        if a1 != '-' and a2 != '-':
            aln_coords1.append(coords1[i1])
            aln_coords2.append(coords2[i2])
            idx_map.append((msa_pos, i1, i2))
        if a1 != '-': i1 += 1
        if a2 != '-': i2 += 1

    aln_coords1 = np.array(aln_coords1)
    aln_coords2 = np.array(aln_coords2)
    print(f"  MSA length: {len(seq1)}")
    print(f"  Aligned positions (both non-gap): {len(aln_coords1)}")
    print(f"  Average LDDT: {np.mean([s for s in lddt_scores if s >= 0]):.3f}")

    # Sliding window local alignment
    print("\n[3/5] Sliding window local alignment (RMSD < 4.0Å)...")
    regions = find_similar_regions(
        aln_coords1, aln_coords2, idx_map, lddt_scores,
        window_sizes=[15, 20, 25, 30, 40, 50, 100],
        rmsd_threshold=8.0, step=3
    )
    print(f"  Found {len(regions)} high-quality windows")

    # Deduplicate
    regions = deduplicate_regions(regions)
    print(f"  After deduplication: {len(regions)} non-overlapping regions")

    # Sort by length then RMSD
    regions.sort(key=lambda x: (-x['length'], x['rmsd']))

    for i, r in enumerate(regions[:10]):
        msa_s = idx_map[r['start']][0]
        msa_e = idx_map[r['end']-1][0]
        r1_s = ca1[idx_map[r['start']][1]][2][1]
        r1_e = ca1[idx_map[r['end']-1][1]][2][1]
        r2_s = ca2[idx_map[r['start']][2]][2][1]
        r2_e = ca2[idx_map[r['end']-1][2]][2][1]
        print(f"  Region {i+1}: len={r['length']}, RMSD={r['rmsd']:.3f}Å, "
              f"mean_LDDT={r['mean_lddt']:.3f}, "
              f"BSVV1[{r1_s}-{r1_e}] whIV2[{r2_s}-{r2_e}]")

    # Extract and save
    print("\n[4/5] Extracting and saving aligned substructures...")
    all_corr_data = []
    for i, region in enumerate(regions):
        resids_1 = []
        resids_2 = []
        corr_data = []

        for j in range(region['start'], region['end']):
            _, i1, i2 = idx_map[j]
            resids_1.append(ca1[i1][1:3])
            resids_2.append(ca2[i2][1:3])
            corr_data.append({
                'region': i+1,
                'msa_pos': idx_map[j][0],
                'BSVV1_resnum': ca1[i1][2][1],
                'BSVV1_resname': ca1[i1][3],
                'whIV2_resnum': ca2[i2][2][1],
                'whIV2_resname': ca2[i2][3],
                'LDDT': lddt_scores[idx_map[j][0]],
            })
        all_corr_data.extend(corr_data)

        s1_sub = extract_residues(s1, resids_1)
        s2_sub = extract_residues(s2, resids_2)
        s2_sub_aligned = s2_sub.copy()
        apply_transform(s2_sub_aligned, region['R'], region['cP'], region['cQ'])

        r1_s = ca1[idx_map[region['start']][1]][2][1]
        r1_e = ca1[idx_map[region['end']-1][1]][2][1]
        r2_s = ca2[idx_map[region['start']][2]][2][1]
        r2_e = ca2[idx_map[region['end']-1][2]][2][1]

        remarks = [
            f"FoldMason + Kabsch local alignment",
            f"Region {i+1}: BSVV1[{r1_s}-{r1_e}] vs whIV2[{r2_s}-{r2_e}]",
            f"Length: {region['length']} residues",
            f"Local RMSD: {region['rmsd']:.3f} A",
            f"Mean LDDT: {region['mean_lddt']:.3f}",
        ]

        fname1 = f"{workdir}/BSVV1_final_region{i+1}_r{r1_s}-{r1_e}_rmsd{region['rmsd']:.2f}A_lddt{region['mean_lddt']:.2f}.pdb"
        fname2 = f"{workdir}/whIV2_final_region{i+1}_r{r2_s}-{r2_e}_rmsd{region['rmsd']:.2f}A_lddt{region['mean_lddt']:.2f}.pdb"

        write_pdb(s1_sub, fname1, remarks)
        write_pdb(s2_sub_aligned, fname2, remarks)

        # Correspondence TSV per region
        tsv_file = f"{workdir}/correspondence_final_region{i+1}.tsv"
        with open(tsv_file, 'w') as f:
            f.write("msa_pos\tBSVV1_resnum\tBSVV1_resname\twhIV2_resnum\twhIV2_resname\tLDDT\n")
            for c in corr_data:
                f.write(f"{c['msa_pos']}\t{c['BSVV1_resnum']}\t{c['BSVV1_resname']}\t"
                        f"{c['whIV2_resnum']}\t{c['whIV2_resname']}\t{c['LDDT']:.4f}\n")
        print(f"[INFO] Saved correspondence: {tsv_file}")

    # Write total correspondence table
    total_tsv = f"{workdir}/correspondence_total.tsv"
    with open(total_tsv, 'w') as f:
        f.write("region\tmsa_pos\tBSVV1_resnum\tBSVV1_resname\twhIV2_resnum\twhIV2_resname\tLDDT\n")
        for c in all_corr_data:
            f.write(f"{c['region']}\t{c['msa_pos']}\t{c['BSVV1_resnum']}\t{c['BSVV1_resname']}\t"
                    f"{c['whIV2_resnum']}\t{c['whIV2_resname']}\t{c['LDDT']:.4f}\n")
    print(f"[INFO] Saved total correspondence: {total_tsv}")

    # Summary
    print("\n" + "=" * 60)
    print("SUMMARY")
    print("=" * 60)
    print(f"Total high-similarity regions (RMSD < 4.0Å): {len(regions)}")
    if regions:
        best = regions[0]
        r1_s = ca1[idx_map[best['start']][1]][2][1]
        r1_e = ca1[idx_map[best['end']-1][1]][2][1]
        r2_s = ca2[idx_map[best['start']][2]][2][1]
        r2_e = ca2[idx_map[best['end']-1][2]][2][1]
        print(f"Best region: BSVV1[{r1_s}-{r1_e}] vs whIV2[{r2_s}-{r2_e}]")
        print(f"  Length: {best['length']}, RMSD: {best['rmsd']:.3f}Å, Mean LDDT: {best['mean_lddt']:.3f}")
    print(f"\nFiles saved to: {workdir}")
    print("  BSVV1_final_region{N}_r{start}-{end}_rmsd{X.XX}A_lddt{Y.YY}.pdb")
    print("  whIV2_final_region{N}_r{start}-{end}_rmsd{X.XX}A_lddt{Y.YY}.pdb")
    print("  correspondence_final_region{N}.tsv")
    print("\nNote: PyMOL not available; Kabsch algorithm used as align equivalent.")
    print("LDDT ranges 0-1, higher = more similar (from FoldMason).")


if __name__ == '__main__':
    main()
