# SDS 20240830

# Modified from script created by Lauren. This version made flexible, so that
# nothing hardcoded in model code.

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

# Updated 20241017 to tryCatch so can catch warnings as well. Removing model
# class from output and instead catching flag based on success, warning, or
# error. Also saving "apVar" which either calculates an approximate
# variance-covariance matrix for all parameters except for fixed effects or
# is set to "Non-positive definite approximate variance-covariance". Finally,
# saving variance of random terms and residual variance to calculate ICC.

# For DAISY testing
# probe <- "cg21870274"
# pheno <- pheno_filt_ct
# sample_var <- "array"
# id_var <- "ID"
# age_var <- "clinage"
# covs <- cov_list

# For TEDDY testing
# probe <- "cg24875017"
# pheno <- pheno_filt_ct
# sample_var <- "rgName"
# id_var <- "maskid"
# age_var <- "sample_agedys"
# covs <- cov_list
# # covs <- NULL

mod.function = function(
    probe, matrix, pheno, sample_var, id_var, age_var, covs=NULL){
  
  #create a temporary dataset of phenotype and methylation with rename probe you
  #want "CpG" moving forward
  matrix <- as.matrix(matrix)
  tmp = as.data.frame(matrix[which(rownames(matrix)==probe),])
  colnames(tmp) = "CpG"
  forAnalysis = merge(pheno, tmp, by.x=sample_var, by.y=0)
  
  # Set fixed portion of model formula using covariates if passed in
  if (!is.null(covs)) {
    form <- as.formula(
      paste(c(paste0("CpG ~ ", age_var), covs), collapse = " + "))
  } else {
    form <- as.formula(paste0("CpG ~ ", age_var))
  }
  
  # Set random portion of model formula
  form_rand <- as.formula(paste("~", age_var, "|", id_var))
  
  #perform statistical model
  # TO NOTE: setting and updating these global variables is not best practice,
  # but is our current best approach to avoid running each model multiple times
  result <- NULL
  flag <- 0
  message <- NA
  tryCatch({
    model = lme(form, random = form_rand, dat=forAnalysis)

    # extract all fixed effects
    model_ext <- broom.mixed::tidy(model) %>%
      dplyr::mutate(
        term = gsub("[()]", "", tolower(term)),
        group = tolower(group),
        effect_group_term = ifelse(
          is.na(group),
          paste0(effect, ".", term),
          paste0(effect, ".", group, ".", term)))
    model_ext_wd <- model_ext %>%
      dplyr::select(-effect, -group, -term) %>%
      pivot_wider(
        names_from = effect_group_term,
        values_from = c(estimate, std.error, df, statistic, p.value),
        names_glue = "{effect_group_term}.{.value}")
    
    # update intercept and slope names to be consistent with prior versions
    model_ext_wd <- model_ext_wd %>%
      dplyr::rename(
        fixed.intercept = fixed.intercept.estimate,
        fixed.slope = paste0("fixed.", age_var, ".estimate"))

    # extract individual random effects
    random.intercept = ranef(model)[1]
    random.slope = ranef(model)[2]

    # get subject specific intercepts
    subject.intercept = random.intercept + model_ext_wd$fixed.intercept
    rownames(subject.intercept) = paste0("intercept_", rownames(subject.intercept))
    subject.intercept = t(subject.intercept)
  
    # get subject specific slopes
    subject.slope = random.slope + model_ext_wd$fixed.slope
    rownames(subject.slope) = paste0("slope_", rownames(subject.slope))
    subject.slope = t(subject.slope)
    
    # other info - could also use glance
    nObs = nobs(model)
    aic = summary(model)$AIC
    bic = summary(model)$BIC
    neg2LL = -2*summary(model)$logLik
    
    # Adding check on variance covariance (Hessian) matrix
    apVar = ifelse(is.character(model$apVar), model$apVar, NA)
    # Also get values from the variance-covariance matrix for random effects
    # since apVar inconsistent between OS. For ICC, also need to extract
    # residual variance
      # TO NOTE: could use getVarCov() to get more significant digits for
      # variance of random intercept and random slope, but not residual variance
    var.rand.int <- as.numeric(VarCorr(model)[1])
    var.rand.slope <- as.numeric(VarCorr(model)[2])
    var.resid <- as.numeric(VarCorr(model)[3])
    varcov.rand.eigen.1 <- eigen(getVarCov(model))$values[1]
    varcov.rand.eigen.2 <- eigen(getVarCov(model))$values[2]

    # TO NOTE: can't use the simpler as.character(summary(model)$call)[2]
    # because just returns the variable name "form"
    model = paste0(
      as.character(summary(model)$call)[1],
      "(", as.character(summary(model)$terms)[2], " ",
      as.character(summary(model)$terms)[1], " ",
      as.character(summary(model)$terms)[3],
      ", random = ", paste(get(model$call$random), collapse = ""),
      ", data = ", as.character(summary(model)$call)[3], ")")
    
    # Update result
    result <- cbind(
      model_ext_wd,
      data.frame(
        model, nObs, aic, bic, neg2LL, apVar, var.rand.int, var.rand.slope,
        var.resid, varcov.rand.eigen.1, varcov.rand.eigen.2,
        subject.intercept, subject.slope))

  }, error = function(e) {
    flag <<- 2
    message <<- e$message
  }, warning = function(w) {
    flag <<- 1
    message <<- w$message
  })
  
  result$CpG <- probe
  result$flag <- flag
  result$message <- message
  result <- as.data.frame(result) %>% dplyr::relocate(CpG, flag, message)
  return(result)
}


# faster cleaner version
mod_function_fast <- function(
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
