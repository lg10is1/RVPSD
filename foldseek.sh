conda activate foldseek
cd $work_dir
mkdir -p tmp

foldseek easy-search $path_to_pdb_files_dir/ \
$compared_database \
$results_name \
tmp --format-output query,target,alntmscore,u,t  
