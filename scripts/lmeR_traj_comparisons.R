compare_age_models_one_cpg <- function(
    probe,
    matrix,
    pheno,
    sample_var,
    id_var,
    age_var,
    covs = NULL,
    control = NULL,
    pred_age_grid = NULL,
    random_slope = FALSE) 
{
  
  fit_model <- function(age_mode, age_df = NULL) {
    mod_function_lmer_traj(
      probe = probe,
      matrix = matrix,
      pheno = pheno,
      sample_var = sample_var,
      id_var = id_var,
      age_var = age_var,
      covs = covs,
      age_mode = age_mode,
      age_df = ifelse(is.null(age_df), 3, age_df),
      return_predictions = TRUE,
      pred_age_grid = pred_age_grid,
      REML = FALSE,
      control = control,
      random_slope = random_slope)
  }
  
  # fit summary rows
  res_linear <- fit_model("linear")
  res_ns3    <- fit_model("ns", 3)
  res_ns4    <- fit_model("ns", 4)
  
  # refit actual model objects for LRT
  build_data <- function() {
    cpg <- as.numeric(matrix[probe, ])
    names(cpg) <- colnames(matrix)
    ph <- pheno
    ph$CpG <- cpg[as.character(ph[[sample_var]])]
    needed <- c("CpG", id_var, age_var, covs)
    needed <- needed[!is.na(needed)]
    ph <- ph[stats::complete.cases(ph[, needed, drop = FALSE]), , drop = FALSE]
    ph[[id_var]] <- as.factor(ph[[id_var]])
    ph
  }
  
  dat <- build_data()
  
  rand_term <- if (random_slope) {
    paste0("(", age_var, " | ", id_var, ")")
  } else {
    paste0("(1 | ", id_var, ")")
  }
  
  cov_txt <- if (length(covs)) paste("+", paste(covs, collapse = " + ")) else ""
  
  f_linear <- stats::as.formula(
    paste0("CpG ~ ", age_var, " ", cov_txt, " + ", rand_term))
  
  f_ns3 <- stats::as.formula(
    paste0("CpG ~ splines::ns(", age_var, ", df = 3) ", cov_txt, " + ", rand_term))
  
  f_ns4 <- stats::as.formula(
    paste0("CpG ~ splines::ns(", age_var, ", df = 4) ", cov_txt, " + ", rand_term))
  
  m_linear <- tryCatch(lme4::lmer(f_linear, data = dat, REML = FALSE, 
                                  control = control), error = function(e) NULL)
  m_ns3    <- tryCatch(lme4::lmer(f_ns3,    data = dat, REML = FALSE, 
                                  control = control), error = function(e) NULL)
  m_ns4    <- tryCatch(lme4::lmer(f_ns4,    data = dat, REML = FALSE, 
                                  control = control), error = function(e) NULL)
  
  lrt_linear_vs_ns3 <- NA
  lrt_ns3_vs_ns4 <- NA
  
  if (!is.null(m_linear) && !is.null(m_ns3)) {
    lrt_linear_vs_ns3 <- tryCatch(anova(m_linear, m_ns3)$`Pr(>Chisq)`[2],
                                  error = function(e) NA_real_)
  }
  
  if (!is.null(m_ns3) && !is.null(m_ns4)) {
    lrt_ns3_vs_ns4 <- tryCatch(anova(m_ns3, m_ns4)$`Pr(>Chisq)`[2], 
                               error = function(e) NA_real_)
  }
  
  bind_rows <- function(...) {
    data.table::rbindlist(list(...), fill = TRUE)
  }
  
  out <- bind_rows(res_linear, res_ns3, res_ns4)
  out$lrt_linear_vs_ns3 <- lrt_linear_vs_ns3
  out$lrt_ns3_vs_ns4 <- lrt_ns3_vs_ns4
  out
}