# TMS 20260322
#
# Modified script created by Lauren & Sarah.
# This version replaces nlme::lme() with lme4::lmer()
# to improve computational speed, numerical stability,
# and scalability for high-throughput EWAS (~700k CpGs).
#
# Updated for spline-based trajectory modeling:
# This version models age using a natural cubic spline with df = 3
# to capture nonlinear early-life CpG trajectories, and returns
# predicted values on a user-specified age grid for downstream
# clustering of CpG sites by trajectory shape.
#
# Arguments:
# probe = single CpG getting summary measures for
# matrix = batch-corrected DNAm M matrix
# pheno = longitudinal phenotype file
# sample_var = name of column in phenotype file with sample names
# sample name format needs to match that in M matrix
# id_var = name of column in phenotype file with subject names
# age_var = name of column in phenotype file with age at sample collection
# covs (optional) = list of desired covariates in format c("cov1", "cov2")
# age_grid (optional) = numeric vector of ages at which to generate
#   population-level predicted CpG values for clustering
# return_predictions = logical; if TRUE, returns predicted values at age_grid
# return_blups = logical; if TRUE, returns subject-level random effects
#   when compatible with chosen random-effects structure
# use_lmerTest = logical; if TRUE, uses lmerTest for Satterthwaite df + p-values
# REML = logical; TRUE (default) for final model estimation,
#   FALSE for ML comparisons of fixed-effects structures
#
# Model:
#   CpG ~ splines::ns(age_var, df = 3) + covs + random effects
#
# Default random-effects structure for clustering applications:
#   CpG ~ splines::ns(age_var, df = 3) + covs + (1 | id_var)
#
# Random-slope model:
#   CpG ~ splines::ns(age_var, df = 3) + covs + (age_var | id_var)
#
# Updated 20260211:
# - Replaced nlme::lme() with lme4::lmer():
#     * More robust handling of near-singular random-effects structures
#     * Employs sparse matrices that improve computational efficiency
#       and numerical stability
#
# - Replaced the single linear age term with a natural spline basis:
#     * Uses splines::ns(age_var, df = 3)
#     * Allows nonlinear modeling of early-life methylation trajectories
#
# - Shifted output from single-slope inference toward trajectory extraction:
#     * Predicted CpG values are generated on a common age grid
#     * These predicted trajectories are intended for downstream clustering
#       of CpGs by similar temporal patterns
#
# - Optional support for lmerTest:
#     * Provides Satterthwaite degrees of freedom
#     * Returns Wald p-values for fixed effects when requested
#     * If lmerTest is not available, asymptotic normal p-values are computed
#
# - Wrapped lmer() in tryCatch + withCallingHandlers to capture both
#   warnings and errors:
#     * fit_status = "ok" / "warning" / "error"
#     * warn_message stores the first warning message (additional warnings muffled)
#     * error_message stores the error message when fitting fails
#
# - Returns compact per-CpG summary row (no model object stored)
#
# - Diagnostics included for model evaluation:
#     * model_call
#     * nObs
#     * AIC, BIC
#     * -2*logLik (neg2LL)
#     * is_singular flag (via lme4::isSingular)
#     * variance components:
#           - random intercept variance
#           - random slope variance (if included)
#           - residual variance
#     * eigenvalues of random-effects variance-covariance matrix
#           (used to detect near-singularity / essentially zero variance)
#
# - Fixed effects output:
#     * Intercept estimate and related statistics
#     * Spline basis coefficient estimates and related statistics
#     * Note: in the spline model there is no single interpretable
#       "age slope" parameter analogous to the linear model
#
# - Predicted trajectory output:
#     * Population-level predictions at each value in age_grid
#     * Returned as columns such as pred_age_0.25, pred_age_0.50, ...
#     * These columns are the primary features for downstream clustering
#
# - Optional: return subject-level BLUPs
#     * Subject-specific random effects can be returned when requested
#     * These are mainly for sensitivity analyses, not the primary
#       clustering workflow
#
# Notes:
# - REML=TRUE is recommended for final genome-wide trajectory estimation
# - Use REML=FALSE only when performing likelihood ratio tests
#   comparing different fixed-effects structures
# - Singular fits (is_singular = TRUE) indicate the random-effects structure
#   may be more complex than needed for that CpG

mod_function_lmer_ns3 <- function(
    probe,
    matrix,
    pheno,
    sample_var,
    id_var,
    age_var,
    covs = NULL,
    age_grid = NULL,
    return_predictions = TRUE,
    use_lmerTest = TRUE,
    REML = TRUE,
    control = NULL,
    random_slope = FALSE) 
{
  
  # ---- basic checks ----
  if (!is.character(probe) || length(probe) != 1) {
    return(data.frame(
      CpG = NA_character_,
      fit_status = "error",
      warn_message = NA_character_,
      error_message = "Probe must be a single string",
      stringsAsFactors = FALSE))
  }
  
  if (!(probe %in% rownames(matrix))) {
    return(data.frame(
      CpG = probe,
      fit_status = "error",
      warn_message = NA_character_,
      error_message = "Probe not found in matrix",
      stringsAsFactors = FALSE))
  }
  
  # ---- assemble analysis data ----
  cpg <- as.numeric(matrix[probe, ])
  names(cpg) <- colnames(matrix)
  
  sid_chr <- as.character(pheno[[sample_var]])
  missing_sids <- setdiff(unique(sid_chr), colnames(matrix))
  if (length(missing_sids) > 0) {
    return(data.frame(
      CpG = probe,
      fit_status = "error",
      warn_message = NA_character_,
      error_message = paste0("Sample IDs in pheno not found in matrix: ", 
                             length(missing_sids)),
      stringsAsFactors = FALSE))
  }
  
  pheno$CpG <- cpg[sid_chr]
  
  needed <- c("CpG", id_var, age_var, covs)
  needed <- needed[!is.na(needed)]
  forAnalysis <- pheno[
    stats::complete.cases(pheno[, needed, drop = FALSE]),
    , drop = FALSE
  ]
  
  if (nrow(forAnalysis) < 5) {
    return(data.frame(
      CpG = probe,
      fit_status = "error",
      warn_message = NA_character_,
      error_message = "Too few complete observations after filtering",
      stringsAsFactors = FALSE))
  }
  
  forAnalysis[[id_var]] <- as.factor(forAnalysis[[id_var]])
  
  # ---- default age grid ----
  if (is.null(age_grid)) {
    age_grid <- quantile(
      forAnalysis[[age_var]],
      probs = seq(0, 1, length.out = 12),
      na.rm = TRUE)
    age_grid <- as.numeric(age_grid)
  }
  
  # ---- build formula ----
  age_term <- paste0("splines::ns(", age_var, ", df = 3)")
  fixed_terms <- c(age_term, covs)
  fixed_rhs <- paste(fixed_terms, collapse = " + ")
  
  rand_term <- if (isTRUE(random_slope)) {
    paste0("(", age_var, " | ", id_var, ")")
  } else {
    paste0("(1 | ", id_var, ")")
  }
  
  form <- stats::as.formula(
    paste0("CpG ~ ", fixed_rhs, " + ", rand_term))
  
  # ---- fit with warning/error capture ----
  fit_status <- "ok"
  warn_message <- NA_character_
  error_message <- NA_character_
  
  ################################################
  #----------- Function Call ---------------------
  ################################################
  fit_fun <- function() {
    if (use_lmerTest && requireNamespace("lmerTest", quietly = TRUE)) {
      lmerTest::lmer(form, data = forAnalysis, REML = REML, control = control)
    } else {
      lme4::lmer(form, data = forAnalysis, REML = REML, control = control)
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
      }),
    error = function(e) {
      fit_status <<- "error"
      error_message <<- conditionMessage(e)
      return(NULL)
    })
  
  if (is.null(model)) {
    return(data.frame(
      CpG = probe,
      fit_status = fit_status,
      warn_message = warn_message,
      error_message = error_message,
      stringsAsFactors = FALSE))
  }
  
  # ---- fixed effects table ----
  sm <- summary(model)
  fe_tab <- as.data.frame(sm$coefficients)
  
  get_col <- function(possible) {
    hit <- intersect(possible, colnames(fe_tab))
    if (length(hit) == 0) return(NA_character_)
    hit[1]
  }
  
  est_col <- get_col(c("Estimate"))
  se_col  <- get_col(c("Std. Error", "Std.Error"))
  t_col   <- get_col(c("t value", "t-value", "t.value"))
  df_col  <- get_col(c("df", "DF"))
  p_col   <- get_col(c("Pr(>|t|)", "Pr(>|z|)", "p.value", "p-value"))
  
  if (is.na(p_col) && !is.na(t_col)) {
    p_col <- "p_asymptotic"
    fe_tab[[p_col]] <- 2 * stats::pnorm(abs(fe_tab[[t_col]]), lower.tail = FALSE)
  }
  
  # keep intercept if present
  int_row <- if ("(Intercept)" %in% rownames(fe_tab)) "(Intercept)" else NA_character_
  
  fixed_intercept <- if (!is.na(int_row)) fe_tab[int_row, est_col] else NA_real_
  se_intercept    <- if (!is.na(int_row)) fe_tab[int_row, se_col] else NA_real_
  t_intercept     <- if (!is.na(int_row) && !is.na(t_col)) fe_tab[int_row, t_col] else NA_real_
  df_intercept    <- if (!is.na(int_row) && !is.na(df_col)) fe_tab[int_row, df_col] else NA_real_
  p_intercept     <- if (!is.na(int_row) && !is.na(p_col)) fe_tab[int_row, p_col] else NA_real_
  
  # optional: keep spline basis coefficients
  ns_rows <- grep(paste0("^splines::ns\\(", age_var, ", df = 3\\)"), rownames(fe_tab))
  
  ns_est <- ns_se <- ns_t <- ns_df <- ns_p <- rep(NA_real_, 3)
  if (length(ns_rows) > 0) {
    take <- seq_len(min(3, length(ns_rows)))
    ns_est[take] <- fe_tab[ns_rows[take], est_col]
    ns_se[take]  <- fe_tab[ns_rows[take], se_col]
    if (!is.na(t_col)) ns_t[take]  <- fe_tab[ns_rows[take], t_col]
    if (!is.na(df_col)) ns_df[take] <- fe_tab[ns_rows[take], df_col]
    if (!is.na(p_col)) ns_p[take]  <- fe_tab[ns_rows[take], p_col]
  }
  
  # ---- diagnostics ----
  nObs   <- stats::nobs(model)
  aic    <- stats::AIC(model)
  bic    <- stats::BIC(model)
  logLik <- as.numeric(stats::logLik(model))
  neg2LL <- -2 * logLik
  
  vc <- lme4::VarCorr(model)
  var.resid <- sigma(model)^2
  singular <- lme4::isSingular(model, tol = 1e-5)
  
  var.rand.int <- NA_real_
  var.rand.slope <- NA_real_
  ev1 <- ev2 <- NA_real_
  
  if (id_var %in% names(vc)) {
    re_sd <- attr(vc[[id_var]], "stddev")
    if (length(re_sd) >= 1) var.rand.int <- re_sd[1]^2
    if (length(re_sd) >= 2) var.rand.slope <- re_sd[2]^2
    
    try({
      G <- as.matrix(vc[[id_var]])
      ev <- eigen(G, symmetric = TRUE, only.values = TRUE)$values
      if (length(ev) >= 1) ev1 <- ev[1]
      if (length(ev) >= 2) ev2 <- ev[2]
    }, silent = TRUE)
  }
  
  out <- data.frame(
    CpG = probe,
    fit_status = fit_status,
    warn_message = warn_message,
    error_message = error_message,
    model_call = paste(deparse(stats::formula(model)), collapse = " "),
    nObs = nObs,
    aic = aic,
    bic = bic,
    logLik = logLik,
    neg2LL = neg2LL,
    
    fixed_intercept = fixed_intercept,
    se_intercept = se_intercept,
    t_intercept = t_intercept,
    df_intercept = df_intercept,
    p_intercept = p_intercept,
    
    ns1_est = ns_est[1],
    ns2_est = ns_est[2],
    ns3_est = ns_est[3],
    ns1_se  = ns_se[1],
    ns2_se  = ns_se[2],
    ns3_se  = ns_se[3],
    ns1_t   = ns_t[1],
    ns2_t   = ns_t[2],
    ns3_t   = ns_t[3],
    ns1_df  = ns_df[1],
    ns2_df  = ns_df[2],
    ns3_df  = ns_df[3],
    ns1_p   = ns_p[1],
    ns2_p   = ns_p[2],
    ns3_p   = ns_p[3],
    
    var.rand.int = var.rand.int,
    var.rand.slope = var.rand.slope,
    var.resid = var.resid,
    varcov.rand.eigen.1 = ev1,
    varcov.rand.eigen.2 = ev2,
    is_singular = singular,
    stringsAsFactors = FALSE
  )
  
  # ---- predicted trajectory on age_grid ----
  if (isTRUE(return_predictions)) {
    pred_df <- data.frame(age_tmp = age_grid)
    names(pred_df)[1] <- age_var
    
    for (cv in covs) {
      x <- forAnalysis[[cv]]
      if (is.factor(x)) {
        pred_df[[cv]] <- factor(levels(x)[1], levels = levels(x))
      } else {
        pred_df[[cv]] <- stats::median(x, na.rm = TRUE)
      }
    }
    
    pred_df[[id_var]] <- forAnalysis[[id_var]][1]
    
    preds <- tryCatch(
      stats::predict(model, newdata = pred_df, re.form = NA, allow.new.levels = TRUE),
      error = function(e) rep(NA_real_, length(age_grid)))
    
    pred_names <- paste0("pred_age_", format(round(age_grid, 3), trim = TRUE))
    pred_row <- as.data.frame(t(preds))
    names(pred_row) <- pred_names
    
    out <- cbind(out, pred_row)
  }
  
  out
}