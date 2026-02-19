# TMS 20260211

# Modified script created by Lauren & Sarah.
# This version replaces nlme::lme() with lme4::lmer()
# to improve computational speed, numerical stability,
# and scalability for high-throughput EWAS (~700k CpGs).

# Arguments:
# probe = single CpG getting summary measures for
# matrix = batch-corrected DNAm M matrix
# pheno = longitudinal phenotype file 
# sample_var = name of column in phenotype file with sample names
# sample name format needs to match that in M matrix
# id_var = name of column in phenotype file with subject names
# age_var = name of column in phenotype file with age at sample collection
# covs (optional) = list of desired covariates in format c("cov1", "cov2")
# return_blups = logical; if TRUE, return subject-level intercepts and slopes
# use_lmerTest = logical; if TRUE, uses lmerTest for Satterthwaite df + p-values
# REML = logical; TRUE (default) for REML estimation, FALSE for ML comparisons

# Model:
#   CpG ~ age_var + covs + (age_var | id_var)

# Updated 20260211:
# - Replaced nlme::lme() with lme4::lmer():
#     * Faster and more memory-efficient for large-scale EWAS
#     * More robust handling of near-singular random-effects structures
#     * Sparse matrix backend improves scalability
#
# - Optional support for lmerTest:
#     * Provides Satterthwaite degrees of freedom
#     * Returns Wald p-values comparable to nlme
#     * If lmerTest not available, asymptotic normal p-values are computed
#
# - Wrapped lmer() in tryCatch + withCallingHandlers to capture both 
#   warnings and errors:
#     * fit_status = "ok" / "warning" / "error"
#     * warn_message stores the first warning message (additional warnings muffled)
#     * error_message stores the error message when fitting fails
#
# - Returns compact per-CpG summary row (no model object stored).
#
# - Diagnostics included for model evaluation:
#     * model_call
#     * nObs
#     * AIC, BIC
#     * -2*logLik (neg2LL)
#     * is_singular flag (via lme4::isSingular)
#     * variance components:
#           - random intercept variance
#           - random slope variance
#           - residual variance
#     * eigenvalues of random-effects variance–covariance matrix
#           (used to detect near-singularity / essentially zero variance)
#
# - Fixed effects output includes:
#     * Estimate (intercept + age slope)
#     * Standard errors
#     * Wald t-statistics
#     * Degrees of freedom (if lmerTest used)
#     * Wald p-values
#
# - Optional: return subject-level BLUPs
#     * subject-specific intercepts and slopes
#     * computed as fixed effect + random effect
#
# Notes:
# - REML=TRUE is recommended for final model estimation.
# - Use REML=FALSE only when performing likelihood ratio tests
#   comparing different fixed-effects structures.
# - For high-dimensionality inference, Wald tests with BH-FDR adjustment
#   across CpGs are typically recommended **Needs further research**.
# - Singular fits (is_singular = TRUE) indicate the random-slope
#   structure for that CpG is unnecessary or variance is essentially zero.


mod_function_lmer <- function(
    probe,
    matrix,
    pheno,
    sample_var,
    id_var,
    age_var,
    covs = NULL,
    return_blups = TRUE,
    use_lmerTest = TRUE, # if TRUE -> Satterthwaite p-values
    REML = TRUE,
    control = NULL,
    ## random_slope = 1 is intercept only model
    ## replace with age_var for random int/slop model
    random_slope = 1) 
  {
  
  # ---- basic checks ----
  if (!is.character(probe) || length(probe) != 1) {
    return(data.frame(CpG = NA_character_, fit_status = "error",
                      warn_message = NA_character_,
                      error_message = "Probe must be a single string",
                      stringsAsFactors = FALSE))
  }
  if (!(probe %in% rownames(matrix))) {
    return(data.frame(CpG = probe, fit_status = "error",
                      warn_message = NA_character_,
                      error_message = "Probe not found in matrix",
                      stringsAsFactors = FALSE))
  }
  
  # ---- assemble analysis data ----
  cpg <- as.numeric(matrix[probe, ])
  names(cpg) <- colnames(matrix)
  
  sid <- pheno[[sample_var]]
  pheno$CpG <- cpg[as.character(sid)]
  
  sid_chr <- as.character(pheno[[sample_var]])
  missing_sids <- setdiff(unique(sid_chr), colnames(matrix))
  if (length(missing_sids) > 0) {
    return(data.frame(
      CpG = probe, fit_status = "error",
      warn_message = NA_character_,
      error_message = paste0("Sample IDs in pheno not found in matrix: ", 
                             length(missing_sids)),
      stringsAsFactors = FALSE
    ))
  }
  
  # keep complete cases for CpG, id, age (and covariates if provided)
  needed <- c("CpG", id_var, age_var, covs)
  needed <- needed[!is.na(needed)]
  forAnalysis <- pheno[stats::complete.cases(pheno[, needed, drop = FALSE]),
                       , drop = FALSE]
  
  if (nrow(forAnalysis) < 3) {
    return(data.frame(CpG = probe, fit_status = "error",
                      warn_message = NA_character_,
                      error_message = "Too few complete observations after filtering",
                      stringsAsFactors = FALSE))
  }
  
  # Ensure grouping factor is a factor (helps lmer)
  forAnalysis[[id_var]] <- as.factor(forAnalysis[[id_var]])
  
  # ---- build formulas ----
  form_fixed <- paste(c(age_var, covs), collapse = " + ")
  # random slope + intercept:
  form <- stats::as.formula(
    # substitute age_var w/ 1 if we want intercept only
    paste0("CpG ~ ", form_fixed, " + (", random_slope , " | ", id_var, ")"))
  
  # ---- fit with clean logging ----
  fit_status <- "ok"
  warn_message <- NA_character_
  error_message <- NA_character_
  model <- NULL
  
  ################################################
  #----------- Function Call ---------------------
  ################################################
  fit_fun <- function() {
    if (use_lmerTest && requireNamespace("lmerTest", quietly = TRUE)) {
      lmerTest::lmer(form, data = forAnalysis, REML = REML, control=control)
    } else {
      lme4::lmer(form, data = forAnalysis, REML = REML, control=control)
    }
  }
  
  model <- tryCatch(
    withCallingHandlers(
      fit_fun(),
      warning = function(w) {
        if (identical(fit_status, "ok")) {
          fit_status <<- "warning"
          warn_message <<- conditionMessage(w)
        }
        invokeRestart("muffleWarning")
      }
    ),
    error = function(e) {
      fit_status <<- "error"
      error_message <<- conditionMessage(e)
      return(NULL)
    }
  )
  
  if (is.null(model)) {
    return(data.frame(CpG = probe, fit_status = fit_status,
                      warn_message = warn_message,
                      error_message = error_message,
                      stringsAsFactors = FALSE))
  }
  
  # ---- fixed effects (Wald tests) ----
  # If lmerTest is used, summary() includes df + Pr(>|t|)
  sm <- summary(model)
  fe_tab <- as.data.frame(sm$coefficients)  # rows = terms
  
  # match column names across lme4 vs lmerTest
  # lme4: Estimate, Std. Error, t value
  # lmerTest: Estimate, Std. Error, df, t value, Pr(>|t|)
  get_col <- function(possible) intersect(possible, colnames(fe_tab))[1]
  
  est_col <- get_col(c("Estimate"))
  se_col <- get_col(c("Std. Error", "Std.Error"))
  t_col <- get_col(c("t value", "t-value", "t.value"))
  df_col <- get_col(c("df", "DF"))
  p_col <- get_col(c("Pr(>|t|)", "Pr(>|z|)", "p.value", "p-value"))
  
  # If no p-values (lme4), compute asymptotic normal p from t-stat
  if (is.na(p_col)) {
    p_col <- "p_asymptotic"
    fe_tab[[p_col]] <- 2 * stats::pnorm(abs(fe_tab[[t_col]]), lower.tail = FALSE)
  }
  
  # pull intercept + age row 
  int_row <- if ("(Intercept)" %in% rownames(fe_tab)) "(Intercept)" else NA_character_
  age_row <- if (age_var %in% rownames(fe_tab)) age_var else NA_character_
  if (is.na(age_row) && tolower(age_var) %in% rownames(fe_tab)) age_row <- tolower(age_var)
  
  fixed_intercept <- if (!is.na(int_row)) fe_tab[int_row, est_col] else NA_real_
  se_intercept <- if (!is.na(int_row)) fe_tab[int_row, se_col] else NA_real_
  t_intercept <- if (!is.na(int_row)) fe_tab[int_row, t_col] else NA_real_
  df_intercept <- if (!is.na(int_row) && !is.na(df_col)) fe_tab[int_row, df_col] else NA_real_
  p_intercept <- if (!is.na(int_row)) fe_tab[int_row, p_col] else NA_real_
  
  fixed_slope <- if (!is.na(age_row)) fe_tab[age_row, est_col] else NA_real_
  se_slope <- if (!is.na(age_row)) fe_tab[age_row, se_col] else NA_real_
  t_slope <- if (!is.na(age_row)) fe_tab[age_row, t_col] else NA_real_
  df_slope <- if (!is.na(age_row) && !is.na(df_col)) fe_tab[age_row, df_col] else NA_real_
  p_slope <- if (!is.na(age_row)) fe_tab[age_row, p_col] else NA_real_
  
  # ---- random effects / BLUPs ----
  subj_int <- subj_slope <- NULL
  if (isTRUE(return_blups)) {
    re <- lme4::ranef(model)[[id_var]]  # data.frame with (Intercept) and age_var columns typically
    
    # column matching
    int_col <- intersect(colnames(re), c("(Intercept)", "Intercept"))[1]
    slope_col <- intersect(colnames(re), c(age_var, tolower(age_var)))[1]
    
    if (is.na(int_col) || is.na(slope_col)) {
      if (identical(fit_status, "ok")) fit_status <- "warning"
      warn_message <- paste0(
        ifelse(is.na(warn_message), "", paste0(warn_message, " | ")),
        "Could not match ranef column names; used column positions."
      )
      int_col <- colnames(re)[1]
      slope_col <- colnames(re)[2]
    }
    
    subject_ids <- rownames(re)
    subject_intercept <- fixed_intercept + re[[int_col]]
    subject_slope <- fixed_slope + re[[slope_col]]
    
    subj_int <- t(stats::setNames(as.numeric(subject_intercept),
                                  paste0("intercept_", subject_ids)))
    subj_slope <- t(stats::setNames(as.numeric(subject_slope),
                                    paste0("slope_", subject_ids)))
  }
  
  # ---- diagnostics ----
  nObs <- stats::nobs(model)
  aic <- stats::AIC(model)
  bic <- stats::BIC(model)
  # neg2LL is typically used for model comparison - may not be applicable here
  neg2LL <- -2 * as.numeric(stats::logLik(model)) 
  
  # lme4 doesn't have apVar like nlme; keep NA for compatibility
  apVar_status <- NA_character_
  
  # variance components
  # VarCorr returns SDs - convert to variances
  vc <- lme4::VarCorr(model)
  sd_int <- attr(vc[[id_var]], "stddev")[1]
  sd_slope <- attr(vc[[id_var]], "stddev")[2]
  var.rand.int <- sd_int^2
  var.rand.slope <- sd_slope^2
  var.resid <- sigma(model)^2
  
  # eigenvalues of random-effects var-cov matrix (2x2)
  ev1 <- ev2 <- NA_real_
  try({
    G <- as.matrix(vc[[id_var]])   # variance-covariance matrix of random effects
    ev <- eigen(G, symmetric = TRUE, only.values = TRUE)$values
    ev1 <- ev[1]; ev2 <- ev[2]
  }, silent = TRUE)
  
  # create a singular fit flag to identify questionable fits
  singular <- lme4::isSingular(model, tol = 1e-5)
  
  out <- data.frame(
    CpG = probe,
    fit_status = fit_status,
    warn_message = warn_message,
    error_message = error_message,
    model_call = paste(deparse(stats::formula(model)), collapse = " "),
    nObs = nObs,
    aic = aic,
    bic = bic,
    neg2LL = neg2LL,
    apVar = apVar_status,
    
    fixed_intercept = fixed_intercept,
    fixed_slope = fixed_slope,
    se_intercept = se_intercept,
    p_intercept = p_intercept,
    se_slope = se_slope,
    p_slope = p_slope,
    t_slope = t_slope,
    df_slope = df_slope,
    
    var.rand.int = var.rand.int,
    var.rand.slope = var.rand.slope,
    var.resid = var.resid,
    varcov.rand.eigen.1 = ev1,
    varcov.rand.eigen.2 = ev2,
    is_singular = singular,
    stringsAsFactors = FALSE
  )
  
  if (isTRUE(return_blups)) {
    out <- cbind(out, subj_int, subj_slope)
  }
  
  out
}
