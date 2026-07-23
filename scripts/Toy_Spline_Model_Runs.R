
suppressPackageStartupMessages({
  library(lme4)
  library(lmerTest)
  library(splines)
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(ggplot2)
  library(readr)
  library(future)
  library(future.apply)
  library(progressr)
  library(qs)
})

set.seed(2522)

############################################################
# User settings
############################################################

output_dir <- "/Users/ryan_summers/GitHub/teddy-dnam-analysis/results"

# load required objects
pheno_complete <- suppressMessages(
  read_csv(
    "/Users/ryan_summers/Library/Mobile\ Documents/com~apple~CloudDocs/GitHub/teddy-dnam-analysis/data/processed/pheno_complete.csv") 
  %>% select(-c(1:2)))

matrix.clean <- qread(
  "/Users/ryan_summers/GitHub/teddy-dnam-analysis/data/processed/matrix.ctrls.qs",
  nthreads = 12)

pheno_complete <- pheno_complete |> 
  mutate(across(c(maskid, cc, gender, CaseControl), as.factor)) |> 
  dplyr::filter(rgName %in% colnames(matrix.clean)) |>  # match samples w/Matrix****
  mutate(across(c(maskid, cc, gender, CaseControl), droplevels))


set.seed(2522)

n_toy_cpgs <- 25000L # ≈5mins for comparisons but ≈30+min for predicted curves

required_probes <- "cg26928153"

available_probes <- setdiff(
  rownames(matrix.clean),
  required_probes)

selected_probes <- unique(c(
  required_probes,
  sample(
    available_probes,
    size = n_toy_cpgs - length(required_probes),
    replace = FALSE
  )
))

selected_probes <- intersect(selected_probes,rownames(matrix.clean))

# Only retain phenotype samples that are needed
toy_sample_ids <- intersect(
  as.character(pheno_complete$rgName),
  colnames(matrix.clean)
)

# Critical: do not send the entire methylation matrix to each worker
matrix.toy <- matrix.clean[selected_probes,toy_sample_ids,drop = FALSE]

message("Toy CpGs: ", nrow(matrix.toy))
message("Toy samples: ", ncol(matrix.toy))


# Spline specifications to compare
spline_df_values <- c(3L, 4L, 5L)

sample_var <- "rgName"
id_var <- "maskid"
age_var <- "age_yrs"

covariates <- c(
  "gender",
  "cc",
  "new_Bcell",
  "new_CD4T",
  "new_CD8T",
  "new_Mono",
  "new_NK"
)

ctrl <- lmerControl(optimizer = "bobyqa",optCtrl = list(maxfun = 2e5))

############################################################
# Prepare phenotype
############################################################

pheno_toy <- pheno_complete %>%
  mutate(
    maskid = as.factor(maskid),
    gender = as.factor(gender),
    cc = as.factor(cc)
  )

# Check that requested probes and samples exist
selected_probes <- intersect(
  selected_probes,
  rownames(matrix.clean)
)

if (length(selected_probes) == 0L) {
  stop("None of the selected CpGs were found in matrix.clean")
}

missing_samples <- setdiff(
  as.character(pheno_toy[[sample_var]]),
  colnames(matrix.clean)
)

if (length(missing_samples) > 0L) {
  stop(
    length(missing_samples),
    " phenotype sample IDs were not found in matrix.clean"
  )
}


fit_one_cpg <- function(probe, spline_dfs = c(3L, 4L)) {
  
  message("Fitting ", probe)
  
  # Attach methylation values by sample ID
  cpg <- as.numeric(matrix.toy[probe, ])
  names(cpg) <- colnames(matrix.toy)
  
  dat <- pheno_toy
  dat$CpG <- cpg[as.character(dat[[sample_var]])]
  
  vars_needed <- c(
    "CpG",
    sample_var,
    id_var,
    age_var,
    covariates
  )
  
  dat <- dat[
    complete.cases(dat[, vars_needed, drop = FALSE]),
    vars_needed,
    drop = FALSE
  ]
  
  dat[[id_var]] <- droplevels(as.factor(dat[[id_var]]))
  
  if (nrow(dat) < 10L) {
    stop("Too few complete observations for ", probe)
  }
  
  if (length(unique(dat[[id_var]])) < 3L) {
    stop("Too few subjects for ", probe)
  }
  
  fixed_covariates <- paste(covariates, collapse = " + ")
  
  model_list <- setNames(
    vector("list", length(spline_dfs)),
    paste0("ns_df", spline_dfs)
  )
  
  for (spline_df in spline_dfs) {
    
    formula_text <- paste0(
      "CpG ~ splines::ns(",
      age_var,
      ", df = ",
      spline_df,
      ") + ",
      fixed_covariates,
      " + (1 | ",
      id_var,
      ")"
    )
    
    model_list[[paste0("ns_df", spline_df)]] <- lmerTest::lmer(
      formula = as.formula(formula_text),
      data = dat,
      
      # Use ML when comparing different fixed effects
      REML = FALSE,
      control = ctrl
    )
  }
  
  list(
    CpG = probe,
    models = model_list
  )
}

toy_results <- setNames(
  vector("list", length(selected_probes)),
  selected_probes
)

############################################################
# Parallel configuration
############################################################

Sys.setenv(
  OMP_NUM_THREADS = "1",
  OPENBLAS_NUM_THREADS = "1",
  MKL_NUM_THREADS = "1",
  VECLIB_MAXIMUM_THREADS = "1"
)

available_cores <- future::availableCores()

# Start conservatively on a laptop
workers <- 12

message("Available cores: ", available_cores)
message("Parallel workers: ", workers)
message(
  "Total models expected: ",
  length(selected_probes) * length(spline_df_values)
)

# multisession works in RStudio and regular R sessions
future::plan(
  future::multisession,
  workers = workers
)

# Increase only if future reports that exported globals are too large
options(future.globals.maxSize = 8 * 1024^3)

progressr::handlers(
  progressr::handler_txtprogressbar(style = 3,width = 60))

start_time <- Sys.time()

toy_results <- progressr::with_progress({
  
  progress_counter <- progressr::progressor(
    steps = length(selected_probes)
  )
  
  results <- future.apply::future_lapply(
    selected_probes,
    
    function(probe) {
      
      result <- tryCatch(
        fit_one_cpg(
          probe = probe,
          spline_dfs = spline_df_values
        ),
        error = function(e) {
          list(
            CpG = probe,
            error = conditionMessage(e),
            models = NULL
          )
        }
      )
      
      progress_counter(
        message = probe
      )
      
      result
    },
    
    future.seed = TRUE,
    
    future.packages = c(
      "lme4",
      "lmerTest",
      "splines"
    )
  )
  
  results
})

names(toy_results) <- selected_probes

# Return R to sequential processing
future::plan(future::sequential)


# qs::qsave(
#   toy_results,
#   file =file.path(output_dir, "toy_spline_model_objects.qs"),
#   preset = "balanced",
#   nthreads = 12)

error_tbl <- purrr::map_dfr(
  toy_results,
  function(result) {
    
    if (!is.null(result$models)) {
      return(NULL)
    }
    
    data.frame(
      CpG = result$CpG,
      error_message = result$error,
      stringsAsFactors = FALSE
    )
  }
)


qs::qsave(
  error_tbl,
  file =file.path(output_dir, "error_tbl.qs"),
  preset = "balanced",
  nthreads = 12)

message("Model errors recorded: ", nrow(error_tbl))

###### Model Comparison Table ######

extract_model_stats <- function(model, probe, spline_df) {
  
  convergence_messages <- model@optinfo$conv$lme4$messages
  
  if (is.null(convergence_messages)) {
    convergence_messages <- NA_character_
  } else {
    convergence_messages <- paste(
      convergence_messages,
      collapse = "; "
    )
  }
  
  vc <- as.data.frame(VarCorr(model))
  
  random_intercept_var <- vc$vcov[
    vc$grp == id_var &
      vc$var1 == "(Intercept)" &
      is.na(vc$var2)
  ]
  
  if (length(random_intercept_var) == 0L) {
    random_intercept_var <- NA_real_
  }
  
  data.frame(
    CpG = probe,
    spline_df = spline_df,
    n_obs = nobs(model),
    n_subjects = nlevels(model.frame(model)[[id_var]]),
    n_parameters = attr(logLik(model), "df"),
    logLik = as.numeric(logLik(model)),
    AIC = AIC(model),
    BIC = BIC(model),
    residual_variance = sigma(model)^2,
    random_intercept_variance = random_intercept_var[1],
    singular = isSingular(model, tol = 1e-5),
    convergence_message = convergence_messages,
    stringsAsFactors = FALSE
  )
}

comparison_tbl <- map_dfr(
  toy_results,
  function(result) {
    
    if (is.null(result$models)) {
      return(
        data.frame(
          CpG = result$CpG,
          spline_df = NA_integer_,
          error = result$error,
          stringsAsFactors = FALSE
        )
      )
    }
    
    map2_dfr(
      result$models,
      spline_df_values,
      ~ extract_model_stats(
        model = .x,
        probe = result$CpG,
        spline_df = .y
      )
    )
  }
)

comparison_tbl <- comparison_tbl %>%
  group_by(CpG) %>%
  mutate(
    delta_AIC = AIC - min(AIC, na.rm = TRUE),
    delta_BIC = BIC - min(BIC, na.rm = TRUE),
    preferred_by_AIC = AIC == min(AIC, na.rm = TRUE),
    preferred_by_BIC = BIC == min(BIC, na.rm = TRUE)
  ) %>%
  ungroup() %>%
  arrange(CpG, spline_df)


qs::qsave(
  comparison_tbl,
  file =file.path(output_dir, "comparison_tbl.qs"),
  preset = "balanced",
  nthreads = 12)


##### Compare Fitted Trajectories #####

make_prediction_curves <- function(result, grid_length = 100L) 
  {
  # Skip failed CpGs
  if (is.null(result$models)) {
    return(data.frame())
  }
  
  probe <- result$CpG
  
  if (!(probe %in% rownames(matrix.toy))) {
    warning("Probe not found in matrix.toy: ", probe)
    return(data.frame())
  }
  
  ############################################################
  # Reconstruct the analysis data for this CpG
  ############################################################
  
  cpg <- as.numeric(matrix.toy[probe, ])
  names(cpg) <- colnames(matrix.toy)
  
  dat <- pheno_toy
  
  dat$CpG <- cpg[as.character(dat[[sample_var]])]
  
  vars_needed <- c(
    "CpG",
    sample_var,
    id_var,
    age_var,
    covariates
  )
  
  dat <- dat[
    complete.cases(dat[, vars_needed, drop = FALSE]),
    vars_needed,
    drop = FALSE
  ]
  
  dat[[id_var]] <- droplevels(
    as.factor(dat[[id_var]])
  )
  
  ############################################################
  # Validate age before calling seq()
  ############################################################
  
  observed_ages <- dat[[age_var]]
  
  observed_ages <- observed_ages[
    is.finite(observed_ages)
  ]
  
  if (length(observed_ages) < 2L) {
    warning(
      "Insufficient finite ages for ",
      probe
    )
    return(data.frame())
  }
  
  age_min <- min(observed_ages)
  age_max <- max(observed_ages)
  
  if (!is.finite(age_min) || !is.finite(age_max)) {
    warning(
      "Non-finite age range for ",
      probe
    )
    return(data.frame())
  }
  
  if (age_min == age_max) {
    warning(
      "Age has no variability for ",
      probe
    )
    return(data.frame())
  }
  
  age_grid <- seq(
    from = age_min,
    to = age_max,
    length.out = grid_length
  )
  
  ############################################################
  # Predictions from each fitted spline model
  ############################################################
  
  purrr::imap_dfr(
    result$models,
    function(model, model_name) {
      
      spline_df <- suppressWarnings(
        as.integer(
          sub("^ns_df", "", model_name)
        )
      )
      
      purrr::map_dfr(
        age_grid,
        function(a) {
          
          nd <- dat
          nd[[age_var]] <- a
          
          population_predictions <- predict(
            model,
            newdata = nd,
            re.form = NA,
            allow.new.levels = TRUE
          )
          
          data.frame(
            CpG = probe,
            model = model_name,
            spline_df = spline_df,
            age_yrs = a,
            predicted_CpG = mean(
              population_predictions,
              na.rm = TRUE
            ),
            stringsAsFactors = FALSE
          )
        }
      )
    }
  )
}

curve_tbl <- map_dfr(toy_results,make_prediction_curves)


qs::qsave(
  curve_tbl,
  file =file.path(output_dir, "toy_spline_prediction_curves.qs"),
  preset = "balanced",
  nthreads = 12)

end_time <- Sys.time()

elapsed_seconds <- as.numeric(
  difftime(end_time, start_time, units = "secs"))

elapsed_minutes <- elapsed_seconds / 60

successful_fits <- sum(
  vapply(toy_results, function(x) !is.null(x$models),logical(1)))

failed_fits <- length(toy_results) - successful_fits

message("")
message("Parallel model fitting completed")
message("CpGs attempted: ", length(toy_results))
message("CpGs successful: ", successful_fits)
message("CpGs failed: ", failed_fits)
message("Elapsed time: ", round(elapsed_minutes, 2), " minutes")
message(
  "Average time per CpG: ",
  round(elapsed_seconds / length(toy_results), 2),
  " seconds")

