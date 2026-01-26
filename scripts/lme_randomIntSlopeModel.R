# TMS 20260124

# Modified script created by Lauren & Sarah. This version was edited to
# to increase computational speed and catch specific error messages

# Arguments:
  # probe = single CpG getting summary measures for
  # matrix = batch-corrected DNAm M matrix
  # pheno = longitudinal phenotype file 
  # sample_var = name of column in phenotype file with sample names
    # sample name format needs to match that in M matrix
  # id_var = name of column in phenotype file with subject names
  # age_var = name of column in phenotype file with age at sample collection
  # covs (optional) = list of desired covariates in format c("cov1", "cov2")

# Model: CpG ~ age_var + covs, random = ~age_var | id_var

# Updated 20260124:
# - Recoded for speed and robustness 
# - Wrapped nlme::lme() in tryCatch + withCallingHandlers to capture both 
#   warnings and errors:
#     * fit_status = "ok" / "warning" / "error"
#     * warn_message stores the first warning message (additional warnings are muffled)
#     * error_message stores the error message when fitting fails
# - Removed returning the model object/class; instead return a compact per-CpG summary row.
# - Added diagnostics to aid troubleshooting:
#     * model_call, nObs, AIC, BIC, -2*logLik
#     * apVar status (e.g., "Non-positive definite approximate variance-covariance")
#     * variance components (random intercept, random slope, residual)
#     * eigenvalues of the random-effects var–cov matrix (to flag near-singularity)
# - Optional: return subject-level BLUPs (intercepts/slopes) by adding random 
#   effects to fixed effects.




# faster cleaner version
mod_function <- function(
    probe,
    matrix,
    pheno,
    sample_var,
    id_var,
    age_var,
    covs = NULL,
    return_blups = TRUE
    ) {
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
  
  # ---- assemble analysis data (fast + stable) ----
  # Build CpG vector with sample ids as names
  cpg <- as.numeric(matrix[probe, ])
  names(cpg) <- colnames(matrix)
  
  # Add CpG to pheno by matching sample IDs (preserves order)
  sid <- pheno[[sample_var]]
  pheno$CpG <- cpg[as.character(sid)]
  
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
  
  # ---- build formulas ----
  rhs <- c(age_var, covs)
  form <- stats::as.formula(paste("CpG ~", paste(rhs, collapse = " + ")))
  form_rand <- stats::as.formula(paste("~", age_var, "|", id_var))
  
  # ---- fit with clean logging ----
  fit_status <- "ok"
  warn_message <- NA_character_
  error_message <- NA_character_
  model <- NULL
  
  model <- tryCatch(
    withCallingHandlers(
      nlme::lme(form, random = form_rand, data = forAnalysis,
                na.action = na.omit),
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
  
  # ---- fixed effects (fast) ----
  tt <- summary(model)$tTable
  # tt has rows = terms, cols = Value, Std.Error, DF, t-value, p-value
  get_row <- function(term) {
    if (term %in% rownames(tt)) term else NA_character_
  }
  
  # Intercept might be "(Intercept)" or "Intercept" depending on method
  int_row <- get_row("(Intercept)")
  if (is.na(int_row)) int_row <- get_row("Intercept")
  
  age_row <- get_row(age_var)
  # Some pipelines lower-case age var or rename; fallback tries
  if (is.na(age_row)) age_row <- get_row(tolower(age_var))
  
  fixed_intercept <- if (!is.na(int_row)) tt[int_row, "Value"] else NA_real_
  fixed_slope <- if (!is.na(age_row)) tt[age_row, "Value"] else NA_real_
  
  # Optional: store SE/p for intercept and slope
  se_intercept <- if (!is.na(int_row)) tt[int_row, "Std.Error"] else NA_real_
  p_intercept  <- if (!is.na(int_row)) tt[int_row, "p-value"] else NA_real_
  se_slope      <- if (!is.na(age_row)) tt[age_row, "Std.Error"] else NA_real_
  p_slope       <- if (!is.na(age_row)) tt[age_row, "p-value"] else NA_real_
  
  # ---- random effects / BLUPs ----
  subj_int <- subj_slope <- NULL
  if (isTRUE(return_blups)) {
    re <- nlme::ranef(model)
    
    # robust column matching
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
  neg2LL <- -2 * as.numeric(stats::logLik(model))
  
  # apVar status
  apVar_status <- if (is.character(model$apVar)) model$apVar else NA_character_
  
  # variance components
  vc <- nlme::VarCorr(model)
  var.rand.int <- suppressWarnings(as.numeric(vc[1, "Variance"]))
  var.rand.slope <- suppressWarnings(as.numeric(vc[2, "Variance"]))
  var.resid <- suppressWarnings(as.numeric(vc[nrow(vc), "Variance"]))
  
  # eigenvalues of random-effects var-cov 
  ev1 <- ev2 <- NA_real_
  try({
    G <- as.matrix(nlme::getVarCov(model))
    ev <- eigen(G, symmetric = TRUE, only.values = TRUE)$values
    ev1 <- ev[1]; ev2 <- ev[2]
  }, silent = TRUE)
  
  out <- data.frame(
    CpG = probe,
    fit_status = fit_status,
    warn_message = warn_message,
    error_message = error_message,
    model_call = paste(deparse(model$call), collapse = " "),
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
    var.rand.int = var.rand.int,
    var.rand.slope = var.rand.slope,
    var.resid = var.resid,
    varcov.rand.eigen.1 = ev1,
    varcov.rand.eigen.2 = ev2,
    stringsAsFactors = FALSE
  )
  
  if (isTRUE(return_blups)) {
    out <- cbind(out, subj_int, subj_slope)
  }
  
  out
}
