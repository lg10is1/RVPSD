#!/usr/bin/env python3
"""
Parse FoldMason JSON LDDT report, extract high-similarity regions,
align with Kabsch, and save as PDB with correspondence annotations.

LDDT (Local Distance Difference Test): higher = more similar (0-1 scale)
"""

import json
import numpy as np
from Bio.PDB.MMCIFParser import MMCIFParser
from Bio.PDB.PDBIO import PDBIO
from Bio.PDB import Select, Structure, Model, Chain


def parse_ca_string(ca_str):
    """Parse comma-separated CA coordinates from FoldMason JSON."""
    values = [float(v) for v in ca_str.split(',') if v.strip()]
    coords = []
    for i in range(0, len(values), 3):
        if i + 2 < len(values):
            coords.append([values[i], values[i+1], values[i+2]])
    return np.array(coords)


def get_structure_residues(structure):
    """Get all residues with CA atoms. Returns list of (chain_id, res_id, resname)."""
    residues = []
    for model in structure:
        for chain in model:
            for residue in chain:
                if 'CA' in residue:
                    residues.append((chain.id, residue.id, residue.resname))
    return residues


def kabsch(P, Q):
    """Kabsch algorithm for optimal superposition."""
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
    """Apply rotation and translation to all atoms."""
    for model in structure:
        for chain in model:
            for residue in chain:
                for atom in residue:
                    atom.set_coord((atom.get_coord() - cQ) @ R + cP)
    return structure


def write_pdb(structure, filename, remark_lines=None):
    """Write structure to PDB file with optional REMARK lines."""
    class AllSelect(Select):
        def accept_atom(self, atom):
            return True
    io = PDBIO()
    io.set_structure(structure)
    with open(filename, 'w') as f:
        if remark_lines:
            for line in remark_lines:
                f.write(f"REMARK   1 {line}\n")
        io.save(f, select=AllSelect())
    print(f"[INFO] Saved: {filename}")


def extract_residues(structure, residue_ids):
    """Extract specified residues into a new structure."""
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


def find_high_lddt_regions(scores, threshold=0.7, min_length=10):
    """
    Find continuous regions where LDDT > threshold.
    scores: array of LDDT values (may contain -1 for gaps)
    Returns list of (start, end, mean_lddt, max_lddt) tuples.
    """
    regions = []
    in_region = False
    start = 0

    for i, s in enumerate(scores):
        if s >= threshold:
            if not in_region:
                start = i
                in_region = True
        else:
            if in_region:
                end = i
                length = end - start
                if length >= min_length:
                    region_scores = scores[start:end]
                    regions.append({
                        'start': start,
                        'end': end,
                        'length': length,
                        'mean_lddt': np.mean(region_scores),
                        'max_lddt': np.max(region_scores),
                        'min_lddt': np.min(region_scores),
                    })
                in_region = False

    # Handle region at end
    if in_region:
        end = len(scores)
        length = end - start
        if length >= min_length:
            region_scores = scores[start:end]
            regions.append({
                'start': start,
                'end': end,
                'length': length,
                'mean_lddt': np.mean(region_scores),
                'max_lddt': np.max(region_scores),
                'min_lddt': np.min(region_scores),
            })

    return regions


def map_msa_to_structure_indices(msa_seq, struct_residues):
    """
    Map MSA positions (0-based) to structure residue indices (0-based).
    Returns dict: msa_pos -> struct_idx, for non-gap positions.
    """
    mapping = {}
    struct_idx = 0
    for msa_pos, aa in enumerate(msa_seq):
        if aa != '-':
            if struct_idx < len(struct_residues):
                mapping[msa_pos] = struct_idx
                struct_idx += 1
    return mapping


def main():
    workdir = "./"
    json_file = f"{workdir}/foldmason_result.a3m.json"
    cif1 = f"{workdir}/BSVV1_PP480790.1.cif"
    cif2 = f"{workdir}/whIV2_AJG39294.cif"

    print("=" * 60)
    print("FoldMason LDDT-based High-Similarity Region Extraction")
    print("=" * 60)
    print("\nNote: LDDT ranges 0-1. Higher LDDT = higher structural similarity.")

    # Load JSON
    print("\n[1/5] Loading FoldMason LDDT report...")
    with open(json_file) as f:
        data = json.load(f)

    scores = np.array(data['scores'])
    entry0 = data['entries'][0]
    entry1 = data['entries'][1]
    seq1 = entry0['aa']
    seq2 = entry1['aa']

    print(f"  MSA columns: {len(scores)}")
    print(f"  Average LDDT: {scores[scores >= 0].mean():.3f}")
    print(f"  Max LDDT: {scores[scores >= 0].max():.3f}")
    print(f"  Columns with LDDT > 0.7: {(scores > 0.7).sum()}")
    print(f"  Columns with LDDT > 0.8: {(scores > 0.8).sum()}")
    print(f"  Columns with LDDT > 0.9: {(scores > 0.9).sum()}")

    # Parse structures
    print("\n[2/5] Parsing CIF structures...")
    parser = MMCIFParser()
    s1 = parser.get_structure('BSVV1', cif1)
    s2 = parser.get_structure('whIV2', cif2)

    res1 = get_structure_residues(s1)
    res2 = get_structure_residues(s2)
    print(f"  BSVV1: {len(res1)} residues")
    print(f"  whIV2: {len(res2)} residues")

    # Map MSA to structure
    print("\n[3/5] Mapping MSA positions to structure residues...")
    map1 = map_msa_to_structure_indices(seq1, res1)
    map2 = map_msa_to_structure_indices(seq2, res2)

    # Find high-LDDT regions
    print("\n[4/5] Finding high-LDDT regions (threshold >= 0.7)...")
    regions = find_high_lddt_regions(scores, threshold=0.7, min_length=10)

    # Sort by length descending
    regions.sort(key=lambda x: -x['length'])

    print(f"  Found {len(regions)} regions with LDDT >= 0.7 and length >= 10")
    for i, r in enumerate(regions[:15]):
        # Get corresponding structure residue numbers
        msa_s = r['start']
        msa_e = r['end'] - 1
        if msa_s in map1 and msa_e in map1:
            r1_s = res1[map1[msa_s]][1][1]  # residue number
            r1_e = res1[map1[msa_e]][1][1]
        else:
            r1_s = r1_e = "N/A"
        if msa_s in map2 and msa_e in map2:
            r2_s = res2[map2[msa_s]][1][1]
            r2_e = res2[map2[msa_e]][1][1]
        else:
            r2_s = r2_e = "N/A"

        print(f"  Region {i+1}: MSA[{r['start']}-{r['end']}] len={r['length']}, "
              f"mean_LDDT={r['mean_lddt']:.3f}, max_LDDT={r['max_lddt']:.3f}, "
              f"BSVV1[{r1_s}-{r1_e}] whIV2[{r2_s}-{r2_e}]")

    # Extract and align
    print("\n[5/5] Extracting and aligning high-similarity substructures...")

    # Parse CA coordinates from JSON for alignment
    ca_coords1 = parse_ca_string(entry0['ca'])
    ca_coords2 = parse_ca_string(entry1['ca'])

    for i, region in enumerate(regions[:5]):  # Top 5 regions
        msa_s = region['start']
        msa_e = region['end']

        # Get structure indices for this region
        struct_idx1 = []
        struct_idx2 = []
        valid_msa_pos = []

        for msa_pos in range(msa_s, msa_e):
            if scores[msa_pos] >= 0 and msa_pos in map1 and msa_pos in map2:
                struct_idx1.append(map1[msa_pos])
                struct_idx2.append(map2[msa_pos])
                valid_msa_pos.append(msa_pos)

        if len(struct_idx1) < 5:
            print(f"  Region {i+1}: too few valid residues ({len(struct_idx1)}), skipping")
            continue

        # Get CA coordinates for alignment
        P = ca_coords1[struct_idx1]
        Q = ca_coords2[struct_idx2]

        # Kabsch alignment
        R, cP, cQ = kabsch(P, Q)
        Q_aligned = (Q - cQ) @ R + cP
        rmsd = np.sqrt(np.mean(np.sum((P - Q_aligned)**2, axis=1)))

        # Extract residues from original structures
        resids_1 = [res1[idx][0:2] for idx in struct_idx1]
        resids_2 = [res2[idx][0:2] for idx in struct_idx2]

        s1_sub = extract_residues(s1, resids_1)
        s2_sub = extract_residues(s2, resids_2)
        s2_sub_aligned = s2_sub.copy()
        apply_transform(s2_sub_aligned, R, cP, cQ)

        # Residue numbers for filename
        r1_s = res1[struct_idx1[0]][1][1]
        r1_e = res1[struct_idx1[-1]][1][1]
        r2_s = res2[struct_idx2[0]][1][1]
        r2_e = res2[struct_idx2[-1]][1][1]

        # Remarks with correspondence info
        remarks = [
            f"FoldMason LDDT-based alignment",
            f"Region {i+1}: MSA positions {msa_s}-{msa_e-1}",
            f"BSVV1 residues: {r1_s}-{r1_e} ({len(struct_idx1)} residues)",
            f"whIV2 residues: {r2_s}-{r2_e} ({len(struct_idx2)} residues)",
            f"Mean LDDT: {region['mean_lddt']:.3f}",
            f"Max LDDT: {region['max_lddt']:.3f}",
            f"Aligned RMSD: {rmsd:.3f} A",
            f"Note: LDDT 0-1 scale, higher = more similar",
        ]

        fname1 = f"{workdir}/BSVV1_lddt_region{i+1}_r{r1_s}-{r1_e}_lddt{region['mean_lddt']:.2f}_rmsd{rmsd:.2f}.pdb"
        fname2 = f"{workdir}/whIV2_lddt_region{i+1}_r{r2_s}-{r2_e}_lddt{region['mean_lddt']:.2f}_rmsd{rmsd:.2f}.pdb"

        write_pdb(s1_sub, fname1, remarks)
        write_pdb(s2_sub_aligned, fname2, remarks)

        # Also write a correspondence TSV
        tsv_file = f"{workdir}/correspondence_region{i+1}.tsv"
        with open(tsv_file, 'w') as f:
            f.write("msa_pos\tBSVV1_resnum\tBSVV1_resname\twhIV2_resnum\twhIV2_resname\tLDDT\n")
            for j, msa_pos in enumerate(valid_msa_pos):
                idx1 = struct_idx1[j]
                idx2 = struct_idx2[j]
                f.write(f"{msa_pos}\t{res1[idx1][1][1]}\t{res1[idx1][2]}\t"
                        f"{res2[idx2][1][1]}\t{res2[idx2][2]}\t{scores[msa_pos]:.4f}\n")
        print(f"[INFO] Saved correspondence: {tsv_file}")

    # Save all high-LDDT residues combined
    all_struct_idx1 = []
    all_struct_idx2 = []
    all_msa_pos = []
    for msa_pos in range(len(scores)):
        if scores[msa_pos] >= 0.7 and msa_pos in map1 and msa_pos in map2:
            all_struct_idx1.append(map1[msa_pos])
            all_struct_idx2.append(map2[msa_pos])
            all_msa_pos.append(msa_pos)

    if len(all_struct_idx1) >= 5:
        P_all = ca_coords1[all_struct_idx1]
        Q_all = ca_coords2[all_struct_idx2]
        R_all, cP_all, cQ_all = kabsch(P_all, Q_all)
        Q_all_aligned = (Q_all - cQ_all) @ R_all + cP_all
        rmsd_all = np.sqrt(np.mean(np.sum((P_all - Q_all_aligned)**2, axis=1)))

        resids_all_1 = [res1[idx][0:2] for idx in all_struct_idx1]
        resids_all_2 = [res2[idx][0:2] for idx in all_struct_idx2]

        s1_all = extract_residues(s1, resids_all_1)
        s2_all = extract_residues(s2, resids_all_2)
        s2_all_aligned = s2_all.copy()
        apply_transform(s2_all_aligned, R_all, cP_all, cQ_all)

        remarks_all = [
            f"FoldMason LDDT-based alignment (all LDDT>=0.7 residues)",
            f"Total residues: {len(all_struct_idx1)}",
            f"Aligned RMSD: {rmsd_all:.3f} A",
        ]

        write_pdb(s1_all, f"{workdir}/BSVV1_all_lddt_ge07.pdb", remarks_all)
        write_pdb(s2_all_aligned, f"{workdir}/whIV2_all_lddt_ge07.pdb", remarks_all)

    # Summary
    print("\n" + "=" * 60)
    print("SUMMARY")
    print("=" * 60)
    print(f"High-LDDT regions extracted (LDDT >= 0.7): {len(regions)}")
    if regions:
        best = regions[0]
        print(f"Best region: length={best['length']}, mean_LDDT={best['mean_lddt']:.3f}, max_LDDT={best['max_lddt']:.3f}")
    print(f"Total LDDT>=0.7 residues: {len(all_struct_idx1)}")
    print(f"\nOutput files:")
    print("  - BSVV1_lddt_region{N}_r{start}-{end}_lddt{X.XX}_rmsd{Y.YY}.pdb")
    print("  - whIV2_lddt_region{N}_r{start}-{end}_lddt{X.XX}_rmsd{Y.YY}.pdb")
    print("  - correspondence_region{N}.tsv (residue mapping)")
    print("  - BSVV1/whIV2_all_lddt_ge07.pdb (combined high-LDDT residues)")
    print(f"\nAll saved to: {workdir}")


if __name__ == '__main__':
    main()
