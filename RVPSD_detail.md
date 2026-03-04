# RVPSD Database Resources

This repository contains supporting database resources and search configurations used in the construction and deployment of the RNA Viral Protein Structure Database (RVPSD).

RVPSD is a large-scale structural resource integrating AlphaFold2-predicted RNA viral protein structures with taxonomy, sequence, and functional annotations.

---

## Repository Structure

### 1. blastp_db/

BLASTP search database generated using **BLAST+ v2.16.0**.

This directory contains the formatted BLAST database files used for sequence-based retrieval in RVPSD.

#### Example BLASTP command:

```bash
blastp -query input.fasta        -out output.blast        -db db.fasta        -outfmt 6        -evalue 1e-5        -max_target_seqs 10
```

- `-outfmt 6`: tabular format  
- `-evalue 1e-5`: significance threshold  
- `-max_target_seqs 10`: return top 10 matches  

---

### 2. foldseek_db/

Foldseek structural search database generated using **Foldseek**.

This directory contains the structural database used for structure-based retrieval in RVPSD.

#### Example Foldseek command:

```bash
foldseek easy-search input.pdb                      /database_path/db_name                      output.foldseek                      /tmp_folder_path_name                      --format-output query,target,alntmscore,u,t
```

- `alntmscore`: alignment TM-score  
- Output includes query ID, target ID, and structural similarity metrics  

---

### 3. mysql_db/

MySQL database dump files compatible with **MySQL 5.7**.

This directory contains structured metadata tables used in the RVPSD backend.

#### Example SQL query:

```sql
SELECT * FROM virus WHERE cds_id LIKE '%{xxxx}%';
```

---

## Software Versions

- BLAST+ v2.16.0  
- Foldseek (default stable release)  
- MySQL v5.7  

---

## Data Sources

All viral nucleotide and protein sequences were sourced from public NCBI databases (Riboviria, Taxonomy ID: 2559587).  
Taxonomic classification follows ICTV Virus Metadata Resource (VMR_21-221122_MSL37).

---

## License

Code in this repository is released under the MIT License.  
Associated datasets are released under the CC-BY 4.0 License unless otherwise stated.

---

## RVPSD Web Platform

RVPSD is publicly accessible at:

https://virus.9itsg.net/#/home

---

## Citation

If you use RVPSD resources, please cite:

Qiangzhen Yang, Zhongshuai Tian, Tao Hu, Jiangrong Lou, Hengcong Liu, Edward C. Holmes, Yongyong Shi, Juan Li, Weifeng Shi. RVPSD: A comprehensive and user-friendly web database for RNA viral protein structures. bioRxiv 2026.02.06.704141; doi: https://doi.org/10.64898/2026.02.06.704141

---

## Contact

For questions or technical issues, please contact the corresponding author.
