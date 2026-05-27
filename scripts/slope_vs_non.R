### Is Random Slopes Necessary?

# On spline models

# convert to data.table for efficiency
rs <- as.data.table(lmeR_spline_results)
ri <- as.data.table(lmeR_spline_int_res)

# keep one row per CpG from each result table
rs_sub <- rs[, .(
  CpG,
  aic_slope = aic,
  bic_slope = bic,
  neg2LL_slope = neg2LL,
  fit_status_slope = fit_status,
  is_singular_slope = is_singular,
  var.rand.int_slope = var.rand.int,
  var.rand.slope = var.rand.slope,
  eigen1_slope = varcov.rand.eigen.1,
  eigen2_slope = varcov.rand.eigen.2)]

ri_sub <- ri[, .(
  CpG,
  aic_int = aic,
  bic_int = bic,
  neg2LL_int = neg2LL,
  fit_status_int = fit_status,
  is_singular_int = is_singular,
  var.rand.int_int = var.rand.int,
  eigen1_int = varcov.rand.eigen.1)]

cmp <- merge(rs_sub, ri_sub, by = "CpG", all = FALSE)

# # AIC/BIC differences:
# # positive delta means random slope model is better
# cmp[, delta_AIC := aic_int - aic_slope]
# cmp[, delta_BIC := bic_int - bic_slope]
# 
# # useful flags
# cmp[, slope_model_better_AIC := delta_AIC > 0]
# cmp[, slope_model_better_AIC_2 := delta_AIC > 2]
# cmp[, slope_model_better_AIC_10 := delta_AIC > 10]
# 
# cmp[, slope_model_better_BIC := delta_BIC > 0]
# cmp[, slope_model_better_BIC_2 := delta_BIC > 2]
# cmp[, slope_model_better_BIC_10 := delta_BIC > 10]



summary_tbl <- cmp[, .(
  n_cpg = .N,
  pct_slope_better_AIC = mean(slope_model_better_AIC, na.rm = TRUE) * 100,
  pct_slope_better_AIC_gt2 = mean(slope_model_better_AIC_2, na.rm = TRUE) * 100,
  pct_slope_better_AIC_gt10 = mean(slope_model_better_AIC_10, na.rm = TRUE) * 100,
  
  pct_slope_better_BIC = mean(slope_model_better_BIC, na.rm = TRUE) * 100,
  pct_slope_better_BIC_gt2 = mean(slope_model_better_BIC_2, na.rm = TRUE) * 100,
  pct_slope_better_BIC_gt10 = mean(slope_model_better_BIC_10, na.rm = TRUE) * 100,
  
  mean_delta_AIC = mean(delta_AIC, na.rm = TRUE),
  median_delta_AIC = median(delta_AIC, na.rm = TRUE),
  q25_delta_AIC = quantile(delta_AIC, 0.25, na.rm = TRUE),
  q75_delta_AIC = quantile(delta_AIC, 0.75, na.rm = TRUE),
  
  mean_delta_BIC = mean(delta_BIC, na.rm = TRUE),
  median_delta_BIC = median(delta_BIC, na.rm = TRUE),
  
  pct_singular_slope = mean(is_singular_slope, na.rm = TRUE) * 100,
  pct_eigen2_near0 = mean(eigen2_slope < 1e-6, na.rm = TRUE) * 100)]

round(t(summary_tbl),2)



hist(cmp$delta_AIC,
     breaks = 100,
     main = "Delta AIC: intercept-only minus random-slope",
     xlab = "delta_AIC")
abline(v = c(0, 2, 10), lty = 2)

hist(cmp$delta_BIC,
     breaks = 100,
     main = "Delta BIC: intercept-only minus random-slope",
     xlab = "delta_BIC")

