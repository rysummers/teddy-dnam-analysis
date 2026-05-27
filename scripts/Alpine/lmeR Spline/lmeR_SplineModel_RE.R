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
    return_model = FALSE,
    use_lmerTest = TRUE,
    REML = TRUE,
    control = NULL,
    random_slope = FALSE) 
{
  
  make_error_return <- function(probe_val, msg) {
    summary_row <- data.frame(
      CpG = probe_val,
      fit_status = "error",
      warn_message = NA_character_,
      error_message = msg,
      model_call = NA_character_,
      nObs = NA_integer_,
      aic = NA_real_,
      bic = NA_real_,
      logLik = NA_real_,
      neg2LL = NA_real_,
      fixed_intercept = NA_real_,
      se_intercept = NA_real_,
      t_intercept = NA_real_,
      df_intercept = NA_real_,
      p_intercept = NA_real_,
      ns1_est = NA_real_,
      ns2_est = NA_real_,
      ns3_est = NA_real_,
      ns1_se = NA_real_,
      ns2_se = NA_real_,
      ns3_se = NA_real_,
      ns1_t = NA_real_,
      ns2_t = NA_real_,
      ns3_t = NA_real_,
      ns1_df = NA_real_,
      ns2_df = NA_real_,
      ns3_df = NA_real_,
      ns1_p = NA_real_,
      ns2_p = NA_real_,
      ns3_p = NA_real_,
      var.resid = NA_real_,
      is_singular = NA,
      stringsAsFactors = FALSE
    )
    
    blup_row <- data.frame(CpG = probe_val, stringsAsFactors = FALSE)
    
    return(list(
      summary_row = summary_row,
      blup_row = blup_row,
      blup_long = data.frame(),
      model = NULL
    ))
  }
  
  if (!is.character(probe) || length(probe) != 1) {
    return(make_error_return(NA_character_, "Probe must be a single string"))
  }
  
  if (!(probe %in% rownames(matrix))) {
    return(make_error_return(probe, "Probe not found in matrix"))
  }
  
  cpg <- as.numeric(matrix[probe, ])
  names(cpg) <- colnames(matrix)
  
  sid_chr <- as.character(pheno[[sample_var]])
  missing_sids <- setdiff(unique(sid_chr), colnames(matrix))
  if (length(missing_sids) > 0) {
    return(make_error_return(
      probe,
      paste0("Sample IDs in pheno not found in matrix: ", length(missing_sids))
    ))
  }
  
  pheno$CpG <- cpg[sid_chr]
  
  needed <- c("CpG", id_var, age_var, covs)
  needed <- needed[!is.na(needed)]
  forAnalysis <- pheno[
    stats::complete.cases(pheno[, needed, drop = FALSE]),
    , drop = FALSE
  ]
  
  if (nrow(forAnalysis) < 5) {
    return(make_error_return(probe, "Too few complete observations after filtering"))
  }
  
  forAnalysis[[id_var]] <- as.factor(forAnalysis[[id_var]])
  
  if (is.null(age_grid)) {
    age_grid <- quantile(
      forAnalysis[[age_var]],
      probs = seq(0, 1, length.out = 12),
      na.rm = TRUE
    )
    age_grid <- as.numeric(age_grid)
  }
  
  age_term <- paste0("splines::ns(", age_var, ", df = 3)")
  fixed_terms <- c(age_term, covs)
  fixed_rhs <- paste(fixed_terms, collapse = " + ")
  
  rand_term <- if (isTRUE(random_slope)) {
    paste0("(", age_term, " | ", id_var, ")")
  } else {
    paste0("(1 | ", id_var, ")")
  }
  
  form_txt <- paste0("CpG ~ ", fixed_rhs, " + ", rand_term)
  form <- stats::as.formula(form_txt)
  
  fit_status <- "ok"
  warn_message <- NA_character_
  error_message <- NA_character_
  
  fit_fun <- function() {
    suppressWarnings({
      if (use_lmerTest && requireNamespace("lmerTest", quietly = TRUE)) {
        lmerTest::lmer(form, data = forAnalysis, REML = REML, control = control)
      } else {
        lme4::lmer(form, data = forAnalysis, REML = REML, control = control)
      }
    })
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
      NULL
    }
  )
  
  if (is.null(model)) {
    out <- make_error_return(probe, error_message)
    out$summary_row$fit_status <- fit_status
    out$summary_row$warn_message <- warn_message
    out$summary_row$error_message <- error_message
    return(out)
  }
  
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
  
  int_row <- if ("(Intercept)" %in% rownames(fe_tab)) "(Intercept)" else NA_character_
  
  fixed_intercept <- if (!is.na(int_row)) fe_tab[int_row, est_col] else NA_real_
  se_intercept    <- if (!is.na(int_row)) fe_tab[int_row, se_col] else NA_real_
  t_intercept     <- if (!is.na(int_row) && !is.na(t_col)) fe_tab[int_row, t_col] else NA_real_
  df_intercept    <- if (!is.na(int_row) && !is.na(df_col)) fe_tab[int_row, df_col] else NA_real_
  p_intercept     <- if (!is.na(int_row) && !is.na(p_col)) fe_tab[int_row, p_col] else NA_real_
  
  ns_pattern <- paste0("^splines::ns\\(", age_var, ", df = 3\\)")
  ns_rows <- grep(ns_pattern, rownames(fe_tab))
  
  ns_est <- ns_se <- ns_t <- ns_df <- ns_p <- rep(NA_real_, 3)
  if (length(ns_rows) > 0) {
    take <- seq_len(min(3, length(ns_rows)))
    ns_est[take] <- fe_tab[ns_rows[take], est_col]
    ns_se[take]  <- fe_tab[ns_rows[take], se_col]
    if (!is.na(t_col))  ns_t[take]  <- fe_tab[ns_rows[take], t_col]
    if (!is.na(df_col)) ns_df[take] <- fe_tab[ns_rows[take], df_col]
    if (!is.na(p_col))  ns_p[take]  <- fe_tab[ns_rows[take], p_col]
  }
  
  nObs   <- stats::nobs(model)
  aic    <- stats::AIC(model)
  bic    <- stats::BIC(model)
  logLik <- as.numeric(stats::logLik(model))
  neg2LL <- -2 * logLik
  singular <- lme4::isSingular(model, tol = 1e-5)
  var.resid <- sigma(model)^2
  
  vc <- lme4::VarCorr(model)
  re_var_cols <- list()
  
  if (id_var %in% names(vc)) {
    G <- as.matrix(vc[[id_var]])
    re_terms <- colnames(G)
    
    clean_re_name <- function(x) {
      x <- gsub("^\\(Intercept\\)$", "intercept", x)
      x <- gsub(ns_pattern, "ns", x)
      x <- gsub("[^A-Za-z0-9]+", "", x)
      x
    }
    
    clean_terms <- vapply(re_terms, clean_re_name, character(1))
    
    for (j in seq_along(clean_terms)) {
      re_var_cols[[paste0("var.rand.", clean_terms[j])]] <- G[j, j]
    }
    
    ev <- tryCatch(
      eigen(G, symmetric = TRUE, only.values = TRUE)$values,
      error = function(e) numeric(0)
    )
    if (length(ev) > 0) {
      for (j in seq_along(ev)) {
        re_var_cols[[paste0("varcov.rand.eigen.", j)]] <- ev[j]
      }
    }
  }
  
  summary_row <- data.frame(
    CpG = probe,
    fit_status = fit_status,
    warn_message = warn_message,
    error_message = error_message,
    model_call = form_txt,
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
    var.resid = var.resid,
    is_singular = singular,
    stringsAsFactors = FALSE
  )
  
  if (length(re_var_cols) > 0) {
    summary_row <- cbind(summary_row, as.data.frame(re_var_cols, check.names = FALSE))
  }
  
  if (isTRUE(return_predictions)) {
    pred_vals <- vapply(age_grid, function(a) {
      nd <- forAnalysis
      nd[[age_var]] <- a
      p <- predict(model, newdata = nd, re.form = NA, allow.new.levels = TRUE)
      mean(p, na.rm = TRUE)
    }, numeric(1))
    
    pred_names <- paste0("pred_age_", format(round(age_grid, 3), trim = TRUE))
    pred_row <- as.data.frame(as.list(stats::setNames(pred_vals, pred_names)))
    summary_row <- cbind(summary_row, pred_row)
  }
  
  re_subj <- tryCatch(ranef(model)[[id_var]], error = function(e) NULL)
  
  if (is.null(re_subj)) {
    blup_long <- data.frame()
    blup_row <- data.frame(CpG = probe, stringsAsFactors = FALSE)
  } else {
    re_subj[[id_var]] <- rownames(re_subj)
    
    clean_term <- function(x) {
      if (x == "(Intercept)") return("intercept")
      if (grepl(ns_pattern, x)) {
        idx <- sub(ns_pattern, "", x)
        idx <- gsub("[^0-9]", "", idx)
        if (nzchar(idx)) return(paste0("ns", idx))
        return("ns")
      }
      gsub("[^A-Za-z0-9]+", "_", x)
    }
    
    original_terms <- setdiff(colnames(re_subj), id_var)
    clean_terms <- vapply(original_terms, clean_term, character(1))
    colnames(re_subj)[match(original_terms, colnames(re_subj))] <- clean_terms
    
    blup_long_list <- lapply(clean_terms, function(tt) {
      data.frame(
        CpG = probe,
        subject_id = as.character(re_subj[[id_var]]),
        term = tt,
        blup = re_subj[[tt]],
        stringsAsFactors = FALSE
      )
    })
    blup_long <- do.call(rbind, blup_long_list)
    
    # one CpG per row
    blup_vals <- list(CpG = probe)
    for (tt in clean_terms) {
      for (sid in as.character(re_subj[[id_var]])) {
        col_nm <- paste0(tt, "_", sid)
        blup_vals[[col_nm]] <- re_subj[re_subj[[id_var]] == sid, tt][1]
      }
    }
    blup_row <- as.data.frame(blup_vals, check.names = FALSE, stringsAsFactors = FALSE)
  }
  
  out <- list(
    summary_row = summary_row,
    blup_row = blup_row,
    blup_long = blup_long
  )
  
  if (isTRUE(return_model)) {
    out$model <- model
  }
  
  out
}