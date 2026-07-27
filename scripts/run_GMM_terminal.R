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
  library(matrixStats) 
})


spline_intercepts <- qread(
  "/Users/ryan_summers/GitHub/teddy-dnam-analysis/results/secondary/lmeR_spline_RE_blup_ctrls_Mvalue.qs",
  nthreads=12)

pheno_complete <- suppressMessages(read_csv("/Users/ryan_summers/Library/Mobile\ Documents/com~apple~CloudDocs/GitHub/teddy-dnam-analysis/data/processed/pheno_complete.csv") %>% select(-c(1:2)))


# extract intercept columns
int_cols <- grep("^intercept_", names(spline_intercepts), value = TRUE)


# M-values to Beta ftn
m2beta <- function(m) {
  2^m / (1 + 2^m)}

#### subject-level means for comparison

# m_mat <- qread(
#   "/Users/ryan_summers/GitHub/teddy-dnam-analysis/data/processed/matrix.ctrls.qs",
#   nthreads = 12)
# 
# pheno_complete <- pheno_complete %>%
#   filter(rgName %in% colnames(m_mat))  # match samples w/Matrix****
# 
# 
# ## m_mat: CpGs (rows) x samples (columns), matched to pheno_complete via rgName
# ## pheno_complete: has rgName (sample-level) and maskid (subject-level)
# 
# stopifnot(all(colnames(m_mat) == pheno_complete$rgName))  # confirm alignment before proceeding

## Convert once, up front, to beta -- avoid repeated conversion inside a loop
#beta_mat <- m2beta(m_mat)

# ## Get subject grouping aligned to columns of beta_mat
# maskid_vec <- pheno_complete$maskid[match(colnames(beta_mat), pheno_complete$rgName)]
# 
# unique_subjects <- unique(maskid_vec)

# ## Build the CpG x subject matrix of subject-level means
# subj_mean_beta_mat <- sapply(unique_subjects, function(id) {
#   cols <- which(maskid_vec == id)
#   if (length(cols) == 1) {
#     beta_mat[, cols]
#   } else {
#     rowMeans(beta_mat[, cols, drop = FALSE], na.rm = TRUE)
#   }
# })
# 
# colnames(subj_mean_beta_mat) <- unique_subjects
# rownames(subj_mean_beta_mat) <- rownames(beta_mat)
# 
# dim(subj_mean_beta_mat)  # should be n_CpGs x n_subjects (102, per your earlier numbers)
# 
# t0 <- Sys.time()
# options(future.globals.maxSize = 5 * 1024^3)
# plan(multicore, workers = 12)
# 
# gmm_results <- future_lapply(
#   seq_len(nrow(subj_mean_beta_mat)),
#   function(i) {
#     
#     cpg_name <- rownames(subj_mean_beta_mat)[i]
#     x <- as.numeric(subj_mean_beta_mat[i, ])
#     x <- x[is.finite(x)]
#     
#     # --- insufficient data ---
#     if (length(x) < 20) {
#       return(data.table(
#         CpG = cpg_name,
#         best_G = NA_integer_,
#         bic = NA_real_,
#         mean_uncertainty = NA_real_,
#         max_uncertainty = NA_real_,
#         min_prop = NA_real_,
#         mean_sep_min = NA_real_,
#         sd_min = NA_real_,
#         sd_max = NA_real_,
#         call_type = "insufficient_data"))
#     }
#     
#     # --- near-zero variance -> trivial group 1 ---
#     if (sd(x) < 1e-8) {
#       return(data.table(
#         CpG = cpg_name,
#         best_G = 1L,
#         bic = NA_real_,
#         mean_uncertainty = 0,
#         max_uncertainty = 0,
#         min_prop = 1,
#         mean_sep_min = NA_real_,
#         sd_min = sd(x),
#         sd_max = sd(x),
#         call_type = "zero_variance"))
#     }
#     
#     mod <- tryCatch(
#       Mclust(x, G = 1:3, verbose = FALSE),
#       error = function(e) NULL)
#     
#     # --- Mclust failed ---
#     if (is.null(mod)) {
#       return(data.table(
#         CpG = cpg_name,
#         best_G = NA_integer_,
#         bic = NA_real_,
#         mean_uncertainty = NA_real_,
#         max_uncertainty = NA_real_,
#         min_prop = NA_real_,
#         mean_sep_min = NA_real_,
#         sd_min = NA_real_,
#         sd_max = NA_real_,
#         call_type = "mclust_error"))
#     }
#     
#     means <- sort(as.numeric(mod$parameters$mean))
#     props <- as.numeric(mod$parameters$pro)
#     sds <- tryCatch(
#       sqrt(as.numeric(mod$parameters$variance$sigmasq)),
#       error = function(e) NA_real_)
#     
#     # --- normal fitted case ---
#     data.table(
#       CpG = cpg_name,
#       best_G = mod$G,
#       bic = max(mod$BIC, na.rm = TRUE),
#       mean_uncertainty = mean(mod$uncertainty, na.rm = TRUE),
#       max_uncertainty = max(mod$uncertainty, na.rm = TRUE),
#       min_prop = min(props, na.rm = TRUE),
#       mean_sep_min = if (mod$G > 1) min(diff(means)) else NA_real_,
#       sd_min = min(sds, na.rm = TRUE),
#       sd_max = max(sds, na.rm = TRUE),
#       call_type = "fitted")
#   },
#   future.seed = TRUE)
# 
# gmm_results <- data.table::rbindlist(gmm_results, fill = TRUE)


t0 <- Sys.time()
options(future.globals.maxSize = 5 * 1024^3)
plan(multicore, workers = 12)

#### for intercept features #####
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
  file = '/Users/ryan_summers/Desktop/gmm_results_Mvalue.qs',
  preset = "balanced",
  nthreads = 12)

