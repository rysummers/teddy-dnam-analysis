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
