import os
import subprocess
from multiprocessing import Pool
import argparse

parser = argparse.ArgumentParser(description='Run AlphaFold for the given directory.')
parser.add_argument('--dir_name', type=str, required=True, help='Input the directory name to process')

args = parser.parse_args()
dir_name = args.dir_name

fasta_dir = os.path.join('$pata_to', dir_name)
output_dir = os.path.join('$pata_to', dir_name)


data_dir = '$path_to/alphafold_dbs_new '
model_preset = 'monomer'
max_template_date = '2023-01-01'
progress_file = dir_name + '-fasta_progress.txt'


def check_pdb_file(fasta_file):
    pdb_file_path = os.path.join(output_dir, fasta_file, 'ranked_0.pdb')
    return os.path.exists(pdb_file_path)

def run_alphafold(fasta_path, gpu_id):
    command = ['$path_to/run_alphafold.sh',
               '-i' + fasta_path,
               '-t' + max_template_date,
               '-p' + model_preset,
               '-d' + data_dir,
               '-o' + output_dir,
               '-u' + str(gpu_id),
               '-r none']

    subprocess.run(' '.join(command), shell=True)

    with open(progress_file, 'a') as f:
        f.write(fasta_path + '\n')

fasta_files = [f for f in os.listdir(fasta_dir) if f.endswith('.fasta')]
fasta_files = sorted(fasta_files)
fasta_files=fasta_files[0:]
if os.path.exists(progress_file):
    with open(progress_file, 'r') as f:
        completed_tasks = f.read().splitlines()
    fasta_files = [f for f in fasta_files if f not in completed_tasks and not check_pdb_file(f)]

num_GPUs = 8
pool = Pool(processes=16)
for i, fasta_file in enumerate(fasta_files):
    gpu_id = i % num_GPUs
    fasta_path = os.path.join(fasta_dir, fasta_file)
    pool.apply_async(run_alphafold, args=(fasta_path, gpu_id))

pool.close()
pool.join()

print("All jobs done！")