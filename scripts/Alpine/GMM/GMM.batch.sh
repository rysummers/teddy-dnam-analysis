#!/bin/bash

#SBATCH --nodes=1
#SBATCH --ntasks=20
#SBATCH --partition=amem
#SBATCH --job-name="GMM"
#SBATCH --output="/scratch/alpine/rsummers@xsede.org/teddy_dnam_analysis/logs/%x_%j.out"
#SBATCH --error="/scratch/alpine/rsummers@xsede.org/teddy_dnam_analysis/logs/%x_%j.err"
#SBATCH --account=amc-general
#SBATCH --time=1:00:00
#SBATCH --mem=400G
#SBATCH --qos=mem

# --- Load Miniforge / Mamba ---
module load miniforge

# Activate environment 
mamba activate myenv

# ---- Paths ----
data_dir="/scratch/alpine/rsummers@xsede.org/teddy_dnam_analysis"

int_BLUPs="${data_dir}/results/lmeR_spline_RE_blup_int.qs"
lmer_int_res="${data_dir}/results/lmeR_spline_RE_summary_inty.qs"
out_prefix="${data_dir}/results/GMM"
r_script="${data_dir}/run_GMM.R"

# Make sure logs/results dirs exist
mkdir -p "${data_dir}/logs" "${data_dir}/results"


# Prevent oversubscription 
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1

which Rscript
Rscript -e 'print(.libPaths())'

Rscript "$r_script" "$int_BLUPs" "$lmer_int_res" "$out_prefix" 


