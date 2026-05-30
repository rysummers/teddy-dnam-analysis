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
  library(variancePartition)
  library(BiocParallel)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
  stop("Usage: run_lmeR_model.R <pheno> <matrix_qs2> <out_csv_prefix> [workers]")
}

pheno      <- args[[1]]
matrix_qs2 <- args[[2]]
out_prefix <- args[[3]]

workers <- 20

message("pheno:  ", pheno)
message("matrix_qs: ", matrix_qs2)
message("out_prefix: ", out_prefix)
message("workers:   ", workers)

Sys.setenv(
  OMP_NUM_THREADS = "1",
  OPENBLAS_NUM_THREADS = "1",
  MKL_NUM_THREADS = "1")

pheno_filt_ct <- read_csv(pheno, show_col_types = FALSE)

pheno_filt_ct <- pheno_filt_ct %>%
  mutate(across(c(maskid, cc, gender, CaseControl), as.factor))

pheno_filt_ct$age_yrs <- pheno_filt_ct$sample_agedys / 365.25

matrix.filt <- qread(matrix_qs2, nthreads = min(12, workers))

source("/scratch/alpine/rsummers@xsede.org/teddy_dnam_analysis/lmeR_SplineModel_RE.R")

t0 <- Sys.time()

ctrl <- lme4::lmerControl(
  optimizer = "bobyqa",
  optCtrl = list(maxfun = 2e5))

fit_cpg <- function(probe) {
  mod_function_lmer_ns3(
    probe = probe,
    matrix = matrix.filt,
    pheno = pheno_filt_ct,
    sample_var = "rgName",
    id_var = "maskid",
    age_var = "age_yrs",
    covs = c("gender", "cc", "CaseControl", "new_Bcell", "new_CD4T",
             "new_CD8T", "new_Mono", "new_NK"),
    age_grid = as.numeric(
      quantile(pheno_filt_ct$age_yrs,
               probs = seq(0, 1, length.out = 12),
               na.rm = TRUE)),
    return_predictions = TRUE,
    return_blups = FALSE,
    return_subject_predictions = TRUE,
    subject_pred_age_grid = c(0.2026010, 0.5010267, 0.7501711,0.9938398),
    REML = TRUE,
    control = ctrl,
    random_slope = FALSE)
}

future::plan(future::multicore, workers = workers)
options(future.globals.maxSize = 20 * 1024^3)

test_probes <- rownames(matrix.filt)
test_res <- future_lapply(test_probes, fit_cpg)

summary_tbl <- data.table::rbindlist(
  lapply(test_res, `[[`, "summary_row"),
  fill = TRUE)

blup_tbl <- data.table::rbindlist(
  lapply(test_res, `[[`, "blup_row"),
  fill = TRUE)

subject_pred_tbl <- data.table::rbindlist(
  lapply(test_res, `[[`, "subject_pred_row"),
  fill = TRUE)

t1 <- Sys.time()
message(
  "lme runtime: ",
  round(as.numeric(difftime(t1, t0, units = "mins")), 2),
  " minutes")

qs::qsave(
  summary_tbl,
  file = paste0(out_prefix, "_summary.qs"),
  preset = "balanced",
  nthreads = min(12, workers))

qs::qsave(
  blup_tbl,
  file = paste0(out_prefix, "_blup.qs"),
  preset = "balanced",
  nthreads = min(12, workers))

qs::qsave(
  subject_pred_tbl,
  file = paste0(out_prefix, "_subject_predictions.qs"),
  preset = "balanced",
  nthreads = min(12, workers))


