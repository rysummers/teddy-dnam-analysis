#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(nlme)
  library(lme4)
  library(data.table)
  library(lmerTest)
  library(future.apply)
  library(qs2)
  library(qs)
  library(tidyverse)
  library(SummarizedExperiment)
  library(variancePartition)   # dream()
  library(BiocParallel)        # SnowParam()
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
  stop("Usage: run_lmeR_model.R <pheno> <matrix_qs2> <out_csv> [workers]")
}

pheno  <- args[[1]]
matrix_qs2 <- args[[2]]
out_csv    <- args[[3]]

# workers: number of cores to run in parallel
# careful: too many cores and be too demanding on memory
workers <- 18

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
  mutate(across(c(maskid, cc, gender, CaseControl), as.factor))

# convert age to years for scalability in modeul running
pheno_filt_ct$age_yrs <- pheno_filt_ct$sample_agedys / 365.25

matrix.filt   <- qread(matrix_qs2, nthreads = min(6, workers))

# load R scripts containing needed functions
source("/scratch/alpine/rsummers@xsede.org/teddy_dnam_analysis/lmeR_randomIntSlopeModel.R")


# fit model
t0 <- Sys.time()

# testing different optimizer settings
ctrl <- lme4::lmerControl(
  optimizer = "bobyqa",
  optCtrl = list(maxfun = 2e5)) #1e5 is default

fit_cpg <- function(probe) {
  mod_function_lmer(
    probe = probe,
    matrix = matrix.filt,
    pheno = pheno_filt_ct,
    sample_var = "rgName",    
    id_var = "maskid",           
    age_var = "age_yrs",   
    covs = c("gender", "cc", "CaseControl", "new_Bcell", "new_CD4T", 
             "new_CD8T", "new_Mono", "new_NK"),    
    return_blups = F, # keep FALSE for intercept only modelling
    control = ctrl,
    random_slope = "age_yrs" # references the age by days (yrs) variable
  )}

# set number of cores
future::plan(future::multicore, workers = 12)
options(future.globals.maxSize = 20 * 1024^3)  # 20 GiBb

test_probes <- rownames(matrix.filt) # dim = 790944x1407
test_res <- future_lapply(test_probes, fit_cpg)
test_res <- data.table::rbindlist(test_res, fill = TRUE)

t1 <- Sys.time()
message("lme runtime: ", 
        round(as.numeric(difftime(t1, t0, units = "mins")), 2), 
        " minutes")

# save output to csv
write.csv(test_res, file = out_csv)
message("Saved: ", out_csv)



