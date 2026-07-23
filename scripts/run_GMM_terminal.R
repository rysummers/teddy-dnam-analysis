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
  library(progressr)
})


spline_intercepts <- qread(
  "/Users/ryan_summers/GitHub/teddy-dnam-analysis/results/secondary/lmeR_spline_RE_blup_ctrls.qs",
  nthreads=12)

# extract intercept columns
int_cols <- grep("^intercept_", names(spline_intercepts), value = TRUE)

gmm_results <- qread(
  "/Users/ryan_summers/GitHub/teddy-dnam-analysis/results/GMM_ctrls_results.qs", 
  nthreads = 12)

#remove second bic
gmm_results <- gmm_results[, -4] 

# extract CpGs with NA values from first GMM run
cpg_nas <- gmm_results |> filter(is.na(best_G))

# filter to NA values from GMM
#spline_intercepts <- spline_intercepts[spline_intercepts$CpG %in% cpg_nas$CpG]

# M-values to Beta ftn
m2beta <- function(m) {
  2^m / (1 + 2^m)}

t0 <- Sys.time()
options(future.globals.maxSize = 5 * 1024^3)
plan(multicore, workers = 12)

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


qs::qsave(
  gmm_results,
  file = '/Users/ryan_summers/Desktop/gmm_filtered_results.qs',
  preset = "balanced",
  nthreads = 12)

