### `scripts/01_foldseek_search.sh`


#!/bin/bash
#SBATCH --job-name=AF2_validate
#SBATCH --nodes=1
#SBATCH --cpus-per-task=40
#SBATCH --mem=200G
#SBATCH --time=24:00:00

set -e
source ~/miniconda3/etc/profile.d/conda.sh
conda activate foldseek

BASE_DIR="/results"
work_dir="${BASE_DIR}/work_dir"
DB_TAR="${BASE_DIR}/pdb100.tar.gz"
DB_DIR="${BASE_DIR}/pdb100"
QUERY_DIR="/pdb_files"

mkdir -p ${work_dir}/tmp ${work_dir}/pdb_downloads ${work_dir}/usalign_logs

# Extract pdb100
if [ ! -f "${DB_DIR}/pdb100.dbtype" ]; then
    echo "[STEP 1/7] Extracting pdb100.tar.gz ..."
    mkdir -p ${DB_DIR}
    tar -xzf ${DB_TAR} -C ${DB_DIR}
else
    echo "[STEP 1/7] pdb100 exists, skip"
fi

# Create query DB
if [ ! -f "${work_dir}/queryDB.dbtype" ]; then
    echo "[STEP 2/7] Creating queryDB ..."
    foldseek createdb ${QUERY_DIR} ${work_dir}/queryDB --threads 40
else
    echo "[STEP 2/7] queryDB exists, skip"
fi

# Search
echo "[STEP 3/7] foldseek search ..."
foldseek search ${work_dir}/queryDB ${DB_DIR}/pdb100 ${work_dir}/alnDB ${work_dir}/tmp/search_tmp \
    -s 7.5 -e 10 --max-seqs 1 -a 1 --threads 40

# Convert to m8
echo "[STEP 4/7] Converting to m8 ..."
foldseek convertalis ${work_dir}/queryDB ${DB_DIR}/pdb100 ${work_dir}/alnDB ${work_dir}/foldseek_results.m8 \
    --format-output "query,target,fident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits,prob,qlen,tlen,qaln,taln" \
    --threads 40

echo "[OK] Raw results: ${work_dir}/foldseek_results.m8"

# Compile US-align
USALIGN_BIN="${work_dir}/USalign"
if [ ! -f "${USALIGN_BIN}" ]; then
    echo "[STEP 5/7] Compiling US-align ..."
    cd ${work_dir}
    wget -q https://zhanggroup.org/US-align/bin/module/USalign.cpp -O USalign.cpp
    g++ -O3 -ffast-math -o USalign USalign.cpp
    chmod +x USalign
    echo "[OK] USalign ready"
else
    echo "[STEP 5/7] USalign exists, skip"
fi

echo "[STEP 6/7] Running post-processing ..."
python3 ${work_dir}/02_process_results.py \
    --m8 ${work_dir}/foldseek_results.m8 \
    --query_dir ${QUERY_DIR} \
    --outdir ${work_dir} \
    --usalign ${USALIGN_BIN} \
    --prob_thresh 0.5 \
    --threads 40

echo "[STEP 7/7] DONE"