#!/bin/bash

#SBATCH --nodes=1
#SBATCH --ntasks=20
#SBATCH --partition=amem
#SBATCH --job-name="lme_model"
#SBATCH --output="/scratch/alpine/rsummers@xsede.org/teddy_dnam_analysis/logs/%x_%j.out"
#SBATCH --error="/scratch/alpine/rsummers@xsede.org/teddy_dnam_analysis/logs/%x_%j.err"
#SBATCH --account=amc-general
#SBATCH --time=10:00:00
#SBATCH --mem=400G
#SBATCH --qos=mem

# --- Load Miniforge / Mamba ---
module load miniforge

# Activate environment 
mamba activate myenv

# ---- Paths ----
data_dir="/scratch/alpine/rsummers@xsede.org/teddy_dnam_analysis"

pheno="${data_dir}/pheno_scaled.csv"
matrix_qs2="${data_dir}/matrix.filt.qs2"
out_csv="${data_dir}/results/lme_results.csv"
r_script="${data_dir}/run_lme_model.R"

# Make sure logs/results dirs exist
mkdir -p "${data_dir}/logs" "${data_dir}/results"


# Prevent oversubscription 
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1

which Rscript
Rscript -e 'print(.libPaths())'

Rscript "$r_script" "$pheno" "$matrix_qs2" "$out_csv" 


