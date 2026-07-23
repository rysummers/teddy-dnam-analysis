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
    return_subject_predictions = TRUE,
    subject_pred_age_grid = NULL,
    return_blups = FALSE,
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
      ns1_p = NA_real_,
      ns2_p = NA_real_,
      ns3_p = NA_real_,
      
      age_omnibus_F = NA_real_,
      age_omnibus_numdf = NA_real_,
      age_omnibus_dendf = NA_real_,
      age_omnibus_p = NA_real_,
      
      var.resid = NA_real_,
      is_singular = NA,
      stringsAsFactors = FALSE
    )
    
    blup_row <- data.frame(CpG = probe_val, stringsAsFactors = FALSE)
    # returns two separate tables. 
    # one with typical model stats and one with BLUP estimates
    # another two tables will store the confidence and prediction SEs
    return(list(
      summary_row = summary_row,
      prediction_long = data.frame(),
      blup_row = blup_row,
      blup_long = data.frame(),
      subject_pred_row = data.frame(),
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
  
  # ns_pattern <- paste0("^splines::ns\\(", age_var, ", df = 3\\)")
  # ns_rows <- grep(ns_pattern, rownames(fe_tab))
  # 
  # ns_est <- ns_se <- ns_t <- ns_df <- ns_p <- rep(NA_real_, 3)
  # if (length(ns_rows) > 0) {
  #   take <- seq_len(min(3, length(ns_rows)))
  #   ns_est[take] <- fe_tab[ns_rows[take], est_col]
  #   ns_se[take]  <- fe_tab[ns_rows[take], se_col]
  #   if (!is.na(t_col))  ns_t[take]  <- fe_tab[ns_rows[take], t_col]
  #   if (!is.na(df_col)) ns_df[take] <- fe_tab[ns_rows[take], df_col]
  #   if (!is.na(p_col))  ns_p[take]  <- fe_tab[ns_rows[take], p_col]
  # }
  
  # Match spline coefficient names with or without "splines::"
  ns_pattern <- paste0(
    "^(splines::)?ns\\(",
    age_var,
    ",[[:space:]]*df[[:space:]]*=[[:space:]]*3\\)[1-3]$"
  )
  
  ns_rows <- grep(
    ns_pattern,
    rownames(fe_tab),
    perl = TRUE
  )
  
  # Put coefficients in basis-function order
  if (length(ns_rows) > 0) {
    ns_rows <- ns_rows[
      order(as.integer(
        sub(".*\\)([1-3])$", "\\1", rownames(fe_tab)[ns_rows])
      ))
    ]
  }
  
  ns_est <- ns_se <- ns_t <- ns_df <- ns_p <- rep(NA_real_, 3)
  
  if (length(ns_rows) > 0) {
    take <- seq_len(min(3L, length(ns_rows)))
    
    ns_est[take] <- fe_tab[ns_rows[take], est_col]
    ns_se[take]  <- fe_tab[ns_rows[take], se_col]
    
    if (!is.na(t_col)) {
      ns_t[take] <- fe_tab[ns_rows[take], t_col]
    }
    
    if (!is.na(df_col)) {
      ns_df[take] <- fe_tab[ns_rows[take], df_col]
    }
    
    if (!is.na(p_col)) {
      ns_p[take] <- fe_tab[ns_rows[take], p_col]
    }
  }
  
  ############################################################
  # Omnibus Wald F-test of all spline coefficients
  #
  # H0: ns1 = ns2 = ns3 = 0
  ############################################################
  
  age_omnibus_F <- NA_real_
  age_omnibus_numdf <- NA_real_
  age_omnibus_dendf <- NA_real_
  age_omnibus_p <- NA_real_
  
  if (
    requireNamespace("lmerTest", quietly = T) &&
    length(ns_rows) == 3L
  ) {
    
    coef_names <- names(lme4::fixef(model))
    
    # Match spline coefficients in fixef(model)
    ns_idx <- grep(ns_pattern, coef_names, perl = T)
    
    if (length(ns_idx) == 3L) {
      
      # Put coefficients in basis-function order
      ns_idx <- ns_idx[
        order(as.integer(
          sub(".*\\)([1-3])$", "\\1", coef_names[ns_idx])
        ))
      ]
      
      # Construct a 3 x p joint-contrast matrix
      L <- matrix(
        0,
        nrow = 3L,
        ncol = length(coef_names),
        dimnames = list(
          coef_names[ns_idx],
          coef_names
        )
      )
      
      L[cbind(seq_len(3L), ns_idx)] <- 1
      
      # contestMD requires an lmerTest model
      model_for_test <- if (
        inherits(model, "lmerModLmerTest")
      ) {
        model
      } else {
        tryCatch(
          lmerTest::as_lmerModLmerTest(model),
          error = function(e) NULL
        )
      }
      
      if (!is.null(model_for_test)) {
        age_test <- tryCatch(
          lmerTest::contestMD(
            model_for_test,
            L = L,
            rhs = rep(0, 3L),
            ddf = "Satterthwaite"
          ),
          error = function(e) NULL)
        
        if (!is.null(age_test)) {
          age_omnibus_F <- as.numeric(age_test[["F value"]][1])
          
          age_omnibus_numdf <- as.numeric(age_test[["NumDF"]][1])
          
          age_omnibus_dendf <- as.numeric(age_test[["DenDF"]][1])
          
          age_omnibus_p <- as.numeric(age_test[["Pr(>F)"]][1])
        }
      }
    }
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
    
    age_omnibus_F = age_omnibus_F,
    age_omnibus_numdf = age_omnibus_numdf,
    age_omnibus_dendf = age_omnibus_dendf,
    age_omnibus_p = age_omnibus_p,
    
    var.resid = var.resid,
    is_singular = singular,
    stringsAsFactors = FALSE
  )
  
  if (length(re_var_cols) > 0) {
    summary_row <- cbind(summary_row, as.data.frame(re_var_cols, check.names = FALSE))
  }
  # ############  population level predictions ########### 
  # if (isTRUE(return_predictions)) {
  #   pred_vals <- vapply(age_grid, function(a) {
  #     nd <- forAnalysis
  #     nd[[age_var]] <- a
  #     p <- predict(model, newdata = nd, re.form = NA, allow.new.levels = TRUE)
  #     mean(p, na.rm = TRUE)
  #   }, numeric(1))
  #   
  #   pred_names <- paste0("pred_age_", format(round(age_grid, 3), trim = TRUE))
  #   pred_row <- as.data.frame(as.list(stats::setNames(pred_vals, pred_names)))
  #   summary_row <- cbind(summary_row, pred_row)
  # }
  
  ############################################################
  # Population-level predictions and uncertainty
  ############################################################
  
  prediction_long <- data.frame()
  
  if (isTRUE(return_predictions)) {
    
    beta <- lme4::fixef(model)
    Vbeta <- as.matrix(stats::vcov(model))
    
    # These stored terms preserve the knots used to fit ns(age, df = 3)
    fixed_terms <- stats::delete.response(
      stats::terms(model, fixed.only = TRUE)
    )
    
    original_X <- lme4::getME(model, "X")
    original_contrasts <- attr(original_X, "contrasts")
    
    # Random-effect variance-covariance matrix
    G <- tryCatch(
      as.matrix(lme4::VarCorr(model)[[id_var]]),
      error = function(e) NULL
    )
    
    residual_var <- stats::sigma(model)^2
    
    prediction_list <- lapply(age_grid, function(a) {
      
      nd <- forAnalysis
      nd[[age_var]] <- a
      
      # Population-level predictions, excluding subject BLUPs
      p <- stats::predict(
        model,
        newdata = nd,
        re.form = NA,
        allow.new.levels = TRUE
      )
      
      estimate <- mean(p, na.rm = TRUE)
      
      # Fixed-effect model matrix using the original spline knots
      X <- stats::model.matrix(
        fixed_terms,
        data = nd,
        contrasts.arg = original_contrasts
      )
      
      # Account for columns potentially dropped because of rank deficiency
      missing_beta_cols <- setdiff(names(beta), colnames(X))
      
      if (length(missing_beta_cols) > 0L) {
        stop(
          "Prediction matrix is missing fixed-effect columns: ",
          paste(missing_beta_cols, collapse = ", ")
        )
      }
      
      X <- X[, names(beta), drop = FALSE]
      
      # Average over the observed covariate distribution
      xbar <- colMeans(X)
      
      fixed_var <- as.numeric(
        xbar %*% Vbeta %*% xbar
      )
      
      fixed_var <- max(fixed_var, 0)
      confidence_se <- sqrt(fixed_var)
      
      ##########################################################
      # Random-effect variance for a new subject
      ##########################################################
      
      random_var <- NA_real_
      
      if (!is.null(G)) {
        
        if (!isTRUE(random_slope)) {
          
          # Random-intercept model
          random_var <- as.numeric(G[1, 1])
          
        } else {
          
          # Random intercept plus spline slopes
          z <- rep(0, ncol(G))
          names(z) <- colnames(G)
          
          if ("(Intercept)" %in% names(z)) {
            z["(Intercept)"] <- 1
          }
          
          random_slope_names <- setdiff(
            names(z),
            "(Intercept)"
          )
          
          matched_names <- intersect(
            random_slope_names,
            names(xbar)
          )
          
          z[matched_names] <- xbar[matched_names]
          
          random_var <- as.numeric(
            t(z) %*% G %*% z
          )
        }
      }
      
      random_var <- if (
        is.na(random_var)
      ) {
        NA_real_
      } else {
        max(random_var, 0)
      }
      
      # Prediction interval for one new observation from a new subject
      prediction_var <- fixed_var + random_var + residual_var
      prediction_se <- sqrt(prediction_var)
      
      data.frame(
        CpG = probe,
        age = as.numeric(a),
        
        predicted_CpG = estimate,
        
        confidence_se = confidence_se,
        confidence_lower = estimate - 1.96 * confidence_se,
        confidence_upper = estimate + 1.96 * confidence_se,
        
        random_effect_var = random_var,
        residual_var = residual_var,
        
        prediction_se = prediction_se,
        prediction_lower = estimate - 1.96 * prediction_se,
        prediction_upper = estimate + 1.96 * prediction_se,
        
        stringsAsFactors = FALSE
      )
    })
    
    prediction_long <- data.table::rbindlist(
      prediction_list,
      fill = TRUE
    )
    
    # Optional: retain the original wide predictions in summary_row
    pred_vals <- prediction_long$predicted_CpG
    
    pred_names <- paste0(
      "pred_age_",
      format(round(age_grid, 3), trim = TRUE)
    )
    
    pred_row <- as.data.frame(
      as.list(stats::setNames(pred_vals, pred_names)),
      check.names = FALSE
    )
    
    summary_row <- cbind(summary_row, pred_row)
  }
  ########### subject-level predictions ########### 
  subject_pred_row <- data.frame(CpG = probe, stringsAsFactors = FALSE)
  
  if (isTRUE(return_subject_predictions)) {
    
    if (is.null(subject_pred_age_grid)) {
      subject_pred_age_grid <- age_grid
    }
    
    # one row per subject
    subj_ids <- unique(as.character(forAnalysis[[id_var]]))
    
    for (a in subject_pred_age_grid) {
      
      nd <- forAnalysis[match(subj_ids, as.character(forAnalysis[[id_var]])), , drop = FALSE]
      nd[[age_var]] <- a
      
      p_subj <- predict(model,newdata = nd,re.form = NULL,allow.new.levels = FALSE)
      
      pred_nm <- paste0("subjpred_age_",format(round(a, 3), trim = TRUE),"_")
      
      vals <- as.list(as.numeric(p_subj))
      names(vals) <- paste0(pred_nm, subj_ids)
      
      subject_pred_row <- cbind(
        subject_pred_row,
        as.data.frame(vals, check.names = FALSE))
    }
  }
  
  blup_long <- data.frame()
  blup_row <- data.frame(CpG = probe, stringsAsFactors = FALSE)
  
  if (isTRUE(return_blups)) {
    re_subj <- tryCatch(ranef(model)[[id_var]], error = function(e) NULL)
    
    if (!is.null(re_subj)) {
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
          stringsAsFactors = FALSE)
      })
      blup_long <- do.call(rbind, blup_long_list)
      blup_vals <- list(CpG = probe)
      for (tt in clean_terms) {
        for (sid in as.character(re_subj[[id_var]])) {
          col_nm <- paste0(tt, "_", sid)
          blup_vals[[col_nm]] <- re_subj[re_subj[[id_var]] == sid, tt][1]
        }
      }
      
      blup_row <- as.data.frame(blup_vals, check.names = FALSE, stringsAsFactors = FALSE)
    }
  }

  out <- list(
    summary_row = summary_row,
    prediction_long = prediction_long,
    blup_row = blup_row,
    blup_long = blup_long,
    subject_pred_row = subject_pred_row
  )
  
  if (isTRUE(return_model)) {
    out$model <- model
  }
  
  out
}