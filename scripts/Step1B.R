#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(nlme)
  library(lme4)
  library(doParallel)
  library(broom.mixed)
  library(data.table)
  library(lmerTest)
  library(future.apply)
  library(qs2)
  library(tidyverse)
  library(SummarizedExperiment)
  library(variancePartition)   # dream()
  library(BiocParallel)        # SnowParam()
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
  stop("Usage: Step1B.R <pheno> <matrix_qs2> <out_csv> [workers]")
}

pheno  <- args[[1]]
matrix_qs2 <- args[[2]]
out_csv    <- args[[3]]

# workers: number of cores to run in parallel
# careful: too many cores and be too demanding on memory
workers <- 12

message("pheno:  ", pheno)
message("matrix_qs: ", matrix_qs2)
message("out_csv:    ", out_csv)
message("workers:   ", workers)

# Prevent thread oversubscription (common on HPC)
Sys.setenv(
  OMP_NUM_THREADS = "1",
  OPENBLAS_NUM_THREADS = "1",
  MKL_NUM_THREADS = "1")

# load objects
pheno_filt_ct <- read_csv(pheno)

# convert subject ID, site, and gender to factors
pheno_filt_ct <- pheno_filt_ct %>% 
  mutate(across(c(maskid, cc, gender), as.factor))

matrix   <- qs_read(matrix_qs2, nthreads = min(6, workers))

# load R scripts containing needed functions
source("/scratch/alpine/rsummers@xsede.org/teddy_dnam_analysis/getSubSpecIntSlope.R")

# Get optional 7th argument
cov_list <- c("gender", "new_CD8T", "new_CD4T", "new_NK", "new_Bcell", "new_Mono")


# Run models -------------------------------------------------------------------

# Get CpGs to test
cpgs_to_test <- rownames(matrix)

# Start parallel
num_cores <- 12
cl <- makeCluster(num_cores)
registerDoParallel(cl)

# Run function in parallel over probe groups, want to reduce overhead compared to
# running in parallel over each CpG
cpg_groups <- split(
  cpgs_to_test,
  ceiling(seq_along(cpgs_to_test)/length(cpgs_to_test)*num_cores))

# fit model
t0 <- Sys.time()

model_result <- bind_rows(
  foreach(probe_group = cpg_groups,
          .packages = c("nlme", "dplyr", "tidyr", "broom.mixed")) %dopar% {
            group_results <- lapply(probe_group, function(probe) {
              safe_mod_call_rich(probe, matrix, pheno_filt_ct,
                                 "rgName", "maskid", "age_yrs_c", cov_list)
            })
            do.call(bind_rows, group_results)
          }
)

# Stop parallel
stopCluster(cl)

# Remove row names
model_result <- remove_rownames(model_result)

# Save results
save(model_result, file = "Step1B_Model_scaled.Rdata")
qs_save(model_result, "Step1B_Model_scaled.qs2", nthreads = 6)

t1 <- Sys.time()
message("lme runtime: ", 
        round(as.numeric(difftime(t1, t0, units = "mins")), 2), 
        " minutes")



