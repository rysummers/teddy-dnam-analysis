mod_function_lmer_traj <- function(
    probe,
    matrix,
    pheno,
    sample_var,
    id_var,
    age_var,
    covs = NULL,
    age_mode = c("linear", "ns"),
    age_df = 3,
    return_predictions = T,
    pred_age_grid = NULL,
    use_lmerTest = F,
    REML = F,
    control = NULL,
    random_slope = F) 
{
  
  age_mode <- match.arg(age_mode)
  
  # ---- checks ----
  if (!is.character(probe) || length(probe) != 1) {
    return(data.frame(
      CpG = NA_character_,
      fit_status = "error",
      error_message = "Probe must be a single string",
      stringsAsFactors = F))
  }
  
  if (!(probe %in% rownames(matrix))) {
    return(data.frame(
      CpG = probe,
      fit_status = "error",
      error_message = "Probe not found in matrix",
      stringsAsFactors = F))
  }
  
  # ---- assemble data ----
  cpg <- as.numeric(matrix[probe, ])
  names(cpg) <- colnames(matrix)
  
  sid_chr <- as.character(pheno[[sample_var]])
  missing_sids <- setdiff(unique(sid_chr), colnames(matrix))
  if (length(missing_sids) > 0) {
    return(data.frame(
      CpG = probe,
      fit_status = "error",
      error_message = paste0("Sample IDs in pheno not found in matrix: ", 
                             length(missing_sids)), stringsAsFactors = F))
  }
  
  pheno$CpG <- cpg[sid_chr]
  
  needed <- c("CpG", id_var, age_var, covs)
  needed <- needed[!is.na(needed)]
  forAnalysis <- pheno[stats::complete.cases(pheno[, needed, drop = F]), ,
                       drop = F]
  
  if (nrow(forAnalysis) < 5) {
    return(data.frame(
      CpG = probe,
      fit_status = "error",
      error_message = "Too few complete observations after filtering",
      stringsAsFactors = F))
  }
  
  forAnalysis[[id_var]] <- as.factor(forAnalysis[[id_var]])
  
  # ---- build age term ----
  age_term <- switch(
    age_mode,
    linear = age_var,
    ns = paste0("splines::ns(", age_var, ", df = ", age_df, ")"))
  
  fixed_terms <- c(age_term, covs)
  fixed_terms <- fixed_terms[nzchar(fixed_terms)]
  fixed_rhs <- paste(fixed_terms, collapse = " + ")
  
  rand_term <- if (isTRUE(random_slope)) {
    paste0("(", age_var, " | ", id_var, ")")
  } else {
    paste0("(1 | ", id_var, ")")
  }
  
  form <- stats::as.formula(
    paste0("CpG ~ ", fixed_rhs, " + ", rand_term))
  
  fit_status <- "ok"
  warn_message <- NA_character_
  error_message <- NA_character_
  
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
      stringsAsFactors = F))
  }
  
  # ---- diagnostics ----
  out <- data.frame(
    CpG = probe,
    fit_status = fit_status,
    warn_message = warn_message,
    error_message = error_message,
    age_mode = age_mode,
    age_df = ifelse(age_mode == "ns", age_df, NA),
    random_slope = random_slope,
    model_call = paste(deparse(stats::formula(model)), collapse = " "),
    nObs = stats::nobs(model),
    AIC = stats::AIC(model),
    BIC = stats::BIC(model),
    logLik = as.numeric(stats::logLik(model)),
    neg2LL = -2 * as.numeric(stats::logLik(model)),
    is_singular = lme4::isSingular(model, tol = 1e-5),
    stringsAsFactors = F)
  
  # ---- optional trajectory predictions ----
  if (isTRUE(return_predictions)) {
    
    if (is.null(pred_age_grid)) {
      pred_age_grid <- seq(
        min(forAnalysis[[age_var]], na.rm = T),
        max(forAnalysis[[age_var]], na.rm = T),
        length.out = 10)
    }
    
    # hold nuisance covariates fixed for population trajectory
    pred_df <- data.frame(age_tmp = pred_age_grid)
    names(pred_df)[1] <- age_var
    
    for (cv in covs) {
      x <- forAnalysis[[cv]]
      if (is.factor(x)) {
        pred_df[[cv]] <- factor(levels(x)[1], levels = levels(x))
      } else {
        pred_df[[cv]] <- stats::median(x, na.rm = TRUE)
      }
    }
    
    # include grouping column for predict; allow population-level prediction
    pred_df[[id_var]] <- forAnalysis[[id_var]][1]
    
    preds <- tryCatch(
      stats::predict(model, newdata = pred_df, re.form = NA, allow.new.levels = TRUE),
      error = function(e) rep(NA_real_, length(pred_age_grid)))
    
    pred_names <- paste0("pred_age_", format(round(pred_age_grid, 3), trim = TRUE))
    pred_row <- as.data.frame(t(preds))
    names(pred_row) <- pred_names
    
    out <- cbind(out, pred_row)
  }
  
  out
}