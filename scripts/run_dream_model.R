#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(qs2)
  library(tidyverse)
  library(SummarizedExperiment)
  library(variancePartition)   # dream()
  library(BiocParallel)        # SnowParam()
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
  stop("Usage: run_dream_model.R <pheno> <matrix_qs2> <out_qs2> [workers]")
}

pheno  <- args[[1]]
matrix_qs2 <- args[[2]]
out_qs2    <- args[[3]]

# workers: number of cores to run in parallel
# careful: too many cores and be too demanding on memory
workers <- 10

message("pheno_qs:  ", pheno)
message("matrix_qs: ", matrix_qs2)
message("out_qs:    ", out_qs2)
message("workers:   ", workers)

# Prevent thread oversubscription (common on HPC)
Sys.setenv(
  OMP_NUM_THREADS = "1",
  OPENBLAS_NUM_THREADS = "1",
  MKL_NUM_THREADS = "1")

# load objects
pheno_filt_ct <- read_csv(pheno)
matrix.filt   <- qs_read(matrix_qs2, nthreads = min(6, workers))

# match pheno to matrix columns
pheno_ordered <- pheno_filt_ct[match(colnames(matrix.filt), pheno_filt_ct$rgName), ]
stopifnot(all(!is.na(pheno_ordered$rgName)))

# SummarizedExperiment
se.filt <- SummarizedExperiment(
  assays = list(Mvalue = matrix.filt),
  colData = as.data.frame(pheno_ordered))

# model formula
form <- ~ age_yrs_c + gender + cc + new_Bcell + new_CD4T + new_CD8T + new_Mono + new_NK +
  (1 + age_yrs_c | maskid)

# BiocParallel param
# workers should match number of cores requested in SLURM --ntasks 
param <- SnowParam(workers = workers, type = "SOCK", progressbar = TRUE)

# pre-check that save function will work before dream call
qs2::qs_save(list(test="ok"), out_qs2)
file.remove(out_qs2) # deletes temporary test file if it succeeds 


# fit model
t0 <- Sys.time()
dream.model <- suppressWarnings(
  dream(
    exprObj  = assay(se.filt, "Mvalue"),
    formula  = form,
    data     = colData(se.filt),
    ddf      = "Satterthwaite",
    BPPARAM  = param
  )
)
t1 <- Sys.time()
message("dream runtime: ", round(as.numeric(difftime(t1, t0, units = "mins")), 2), " minutes")

# save output
qs_save(dream.model, out_qs2, nthreads = min(6, workers))
message("Saved: ", out_qs2)

# backup if qs2 doesn't work
saveRDS(dream.model, sub("\\.qs2$", ".rds", out_qs2))


