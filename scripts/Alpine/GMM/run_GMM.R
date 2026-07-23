#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(nlme)
  library(lme4)
  library(lmerTest)
  library(qs2)
  library(qs)
  library(tidyverse)
  library(BiocParallel)
  library(future.apply)
  library(data.table)
  library(mclust)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
  stop("Usage: run_GMM.R <int_BLUPs> <lmer_int_res> <out_prefix> [workers]")
}

int_BLUPs <- args[[1]]
lmer_int_res <- args[[2]]
out_prefix <- args[[3]]

workers <- 20

message("int_BLUPs:  ", int_BLUPs)
message("lmer_int_res: ", lmer_int_res)
message("out_prefix: ", out_prefix)


Sys.setenv(
  OMP_NUM_THREADS = "1",
  OPENBLAS_NUM_THREADS = "1",
  MKL_NUM_THREADS = "1")

# read files in
lmeR_spline_int_BLUP <- qread(
  "/scratch/alpine/rsummers@xsede.org/teddy_dnam_analysis/results/lmeR_spline_RE_blup_ctrls.qs", 
  nthreads = 12)
lmeR_spline_int_res <- qread(
  "/scratch/alpine/rsummers@xsede.org/teddy_dnam_analysis/results/lmeR_spline_RE_summary_ctrls.qs", 
  nthreads = 12)

m2beta <- function(m) {
  2^m / (1 + 2^m)}

# change between time points
# age grid is uneven so we need to divide by age gap
# X.traj = predicted trajectory matrix
# rows = CpGs, columns = predicted ages

# extract trajectory column names
pred_cols <- grep("^pred_age_", names(lmeR_spline_int_res), value = TRUE)
# extract numeric ages
ages <- as.numeric(sub("pred_age_", "", pred_cols))

# order columns by age
ord <- order(ages)
pred_cols <- pred_cols[ord]
ages <- ages[ord]

# predicted M-value matrix: rows = CpGs, columns = ages
Mhat <- as.matrix(lmeR_spline_int_res[, ..pred_cols])
rownames(Mhat) <- lmeR_spline_int_res$CpG
# identify CpGs with complete predicted trajectories
# CpGs with contrast error identified earlier - e.g. class misbalance
good_cpgs <- rowSums(is.na(Mhat)) == 0
Mhat <- Mhat[good_cpgs, , drop = FALSE]

# convert predicted M-values to Beta-values
Betahat <- m2beta(Mhat)

# numerical derivative dBeta/dage between adjacent ages
dBeta_dage <- t(apply(Betahat, 1, function(x) {
  diff(x) / diff(ages)
}))

# name derivative columns by age interval
colnames(dBeta_dage) <- paste0(
  "dBeta_dage_",
  round(ages[-length(ages)], 3),
  "_to_",
  round(ages[-1], 3))


traj_range_beta <- apply(
  Betahat, 1, # remove bad CpGs
  function(x) max(x, na.rm = T) - min(x, na.rm = T))

mean_abs_slope_beta <- rowMeans(abs(dBeta_dage), na.rm = TRUE)

# extract intercept columns
int_cols <- grep("^intercept_", names(lmeR_spline_int_BLUP), value = TRUE)
spline_intercepts <- lmeR_spline_int_BLUP

t0 <- Sys.time()

options(future.globals.maxSize = 5 * 1024^3)
plan(multicore, workers = 12)

gmm_results <- future_lapply(
  seq_len(nrow(spline_intercepts)),
  function(i) {
    
    x <- as.numeric(spline_intercepts[i, ..int_cols])
    x <- x[is.finite(x)]
    
    # --- insufficient data ---
    if (length(x) < 20) {
      return(data.table(
        CpG = spline_intercepts$CpG[i],
        best_G = NA_integer_,
        bic = NA_real_,
        mean_uncertainty = NA_real_,
        max_uncertainty = NA_real_,
        min_prop = NA_real_,
        mean_sep_min = NA_real_,
        sd_min = NA_real_,
        sd_max = NA_real_,
        call_type = "insufficient_data"))
    }
    
    # --- near-zero variance -> trivial group 1 ---
    if (sd(x) < 1e-8) {
      return(data.table(
        CpG = spline_intercepts$CpG[i],
        best_G = 1L,
        bic = NA_real_,
        mean_uncertainty = 0,
        max_uncertainty = 0,
        min_prop = 1,
        mean_sep_min = NA_real_,
        sd_min = sd(x),
        sd_max = sd(x),
        call_type = "zero_variance"))
    }
    
    mod <- tryCatch(
      Mclust(x, G = 1:3, verbose = FALSE),
      error = function(e) NULL)
    
    # --- Mclust failed ---
    if (is.null(mod)) {
      return(data.table(
        CpG = spline_intercepts$CpG[i],
        best_G = NA_integer_,
        bic = NA_real_,
        mean_uncertainty = NA_real_,
        max_uncertainty = NA_real_,
        min_prop = NA_real_,
        mean_sep_min = NA_real_,
        sd_min = NA_real_,
        sd_max = NA_real_,
        call_type = "mclust_error"))
    }
    
    means <- sort(as.numeric(mod$parameters$mean))
    props <- as.numeric(mod$parameters$pro)
    sds <- tryCatch(
      sqrt(as.numeric(mod$parameters$variance$sigmasq)),
      error = function(e) NA_real_)
    
    # --- normal fitted case ---
    data.table(
      CpG = spline_intercepts$CpG[i],
      best_G = mod$G,
      bic = max(mod$BIC, na.rm = TRUE),
      mean_uncertainty = mean(mod$uncertainty, na.rm = TRUE),
      max_uncertainty = max(mod$uncertainty, na.rm = TRUE),
      min_prop = min(props, na.rm = TRUE),
      mean_sep_min = if (mod$G > 1) min(diff(means)) else NA_real_,
      sd_min = min(sds, na.rm = TRUE),
      sd_max = max(sds, na.rm = TRUE),
      call_type = "fitted")
  },
  future.seed = TRUE)

gmm_results <- data.table::rbindlist(gmm_results, fill = TRUE)

plan(sequential)

t1 <- Sys.time()
message(
  "GMM runtime: ",
  round(as.numeric(difftime(t1, t0, units = "mins")), 2),
  " minutes")

gmm_file <- paste0(out_prefix, "_results.qs")


qs::qsave(
  gmm_results,
  file = gmm_file,
  preset = "balanced",
  nthreads = min(12, workers))

message("Saved GMM results: ", gmm_file)


