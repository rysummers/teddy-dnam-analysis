```{r}
benchmark_summary <- benchmark_res %>%
  group_by(age_mode, age_df) %>%
  summarise(
    n = n(),
    mean_AIC = mean(AIC, na.rm = TRUE),
    mean_BIC = mean(BIC, na.rm = TRUE),
    prop_warning = mean(fit_status == "warning", na.rm = TRUE),
    prop_error = mean(fit_status == "error", na.rm = TRUE),
    prop_singular = mean(is_singular, na.rm = TRUE),
    .groups = "drop")

print(benchmark_summary)
```

```{r}
per_cpg_comp <- benchmark_res %>%
  select(CpG, age_mode, age_df, AIC, BIC, lrt_linear_vs_ns3, lrt_ns3_vs_ns4) %>%
  distinct()

# AIC percentage of NS (df=4) being better than linear age fit
aic_pct <- per_cpg_comp %>%
  mutate(age_mode_df = ifelse(age_mode == "linear", 
                              "linear", 
                              paste0(age_mode, "_", age_df))) %>% 
  select(CpG, age_mode_df, AIC) %>%
  distinct() %>%
  pivot_wider(names_from = age_mode_df, values_from = AIC) %>%
  summarise(pct = mean(ns_4 < linear, na.rm = TRUE) * 100) %>%
  pull(pct)

# AIC percentage of NS (df=3) being better than linear age fit
aic_pct2 <- benchmark_res %>%
  mutate(age_mode_df = ifelse(age_mode == "linear", 
                              "linear", 
                              paste0(age_mode, "_", age_df))) %>% 
  select(CpG, age_mode_df, AIC) %>%
  distinct() %>%
  pivot_wider(names_from = age_mode_df, values_from = AIC) %>%
  summarise(pct = mean(ns_3 < linear, na.rm = TRUE) * 100) %>%
  pull(pct)

# LRT - will always be better - more complexity better model fit

aic_pct
aic_pct2
```

```{r}
# results should be in long format with at least:
# CpG, model, AIC
# where model is something like "linear" and "ns3"

aic_change_summary <- benchmark_res %>%
  mutate(age_mode_df = ifelse(age_mode == "linear", 
                              "linear", 
                              paste0(age_mode, "_", age_df))) %>% 
  select(CpG, age_mode_df, AIC) %>%
  distinct() %>%
  pivot_wider(names_from = age_mode_df, values_from = AIC) %>%
  mutate(
    delta_AIC = linear - ns_3,
    ns3_better = delta_AIC > 0,
    ns3_better_2 = delta_AIC > 2,
    ns3_better_10 = delta_AIC > 10) %>%
  summarise(
    n_cpg = n(),
    pct_ns3_better = mean(ns3_better, na.rm = TRUE) * 100,
    pct_delta_AIC_gt_2 = mean(ns3_better_2, na.rm = TRUE) * 100,
    pct_delta_AIC_gt_10 = mean(ns3_better_10, na.rm = TRUE) * 100,
    mean_delta_AIC = mean(delta_AIC, na.rm = TRUE),
    median_delta_AIC = median(delta_AIC, na.rm = TRUE),
    q25_delta_AIC = quantile(delta_AIC, 0.25, na.rm = TRUE),
    q75_delta_AIC = quantile(delta_AIC, 0.75, na.rm = TRUE),
    q90_delta_AIC = quantile(delta_AIC, 0.90, na.rm = TRUE),
    q95_delta_AIC = quantile(delta_AIC, 0.95, na.rm = TRUE),
    max_delta_AIC = max(delta_AIC, na.rm = TRUE))

print(t(aic_change_summary))
```

Based on median, for more than half of CpGs, linear is slightly better or essentially equivalent.
But, a subset of CpGs are strongly nonlinear:
  q90 = 41 
q95 = 87  
max = 1984

~36% show real improvement of AIC greater than 2
~21% show strong improvement of AIC greater than 10

Per CpG AIC change
```{r}
aic_change_tbl <- benchmark_res %>%
  mutate(age_mode_df = ifelse(age_mode == "linear", 
                              "linear", 
                              paste0(age_mode, "_", age_df))) %>% 
  select(CpG, age_mode_df, AIC) %>%
  distinct() %>%
  pivot_wider(names_from = age_mode_df, values_from = AIC) %>%
  mutate(
    delta_AIC = linear - ns_3,
    ns3_better = delta_AIC > 0,
    ns3_better_2 = delta_AIC > 2,
    ns3_better_10 = delta_AIC > 10)

head(aic_change_tbl)
```

```{r}
hist(
  aic_change_tbl$delta_AIC,
  breaks = 100,
  main = "Delta AIC: linear - ns3",
  xlab = "Delta AIC")
abline(v = c(0, 2, 10), lty = 2)
```

```{r}
traj_cols <- grep("^pred_age_", names(benchmark_res), value = TRUE)

traj_mat <- as.matrix(benchmark_res[, ..traj_cols])
rownames(traj_mat) <- benchmark_res$CpG

```

```{r}
# split models
dt_linear <- benchmark_res[age_mode == "linear"]
dt_ns3    <- benchmark_res[age_mode == "ns" & age_df == 3]

# ensure correct order
setkey(dt_linear, CpG)
setkey(dt_ns3, CpG)

all(dt_linear$CpG == dt_ns3$CpG)
```

```{r}
# create matrices
traj_linear <- as.matrix(dt_linear[, ..traj_cols])
traj_ns3    <- as.matrix(dt_ns3[, ..traj_cols])
```

```{r}
cor(as.vector(dist(traj_linear)), as.vector(dist(traj_ns3)))
```

## Alpine Script

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