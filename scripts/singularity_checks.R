### Exploration Convergence Issues

length(lmeR_spline_results$var.rand.slope < 1e-8) # should matrix M-matrix
mean(lmeR_spline_results$var.rand.slope < 1e-8, na.rm=TRUE)
summary(lmeR_spline_results$var.rand.slope)


*	var.rand.int: how much subjects differ in baseline methylation
*	var.rand.slope: how much subjects differ in linear age drift around the spline mean
*	is_singular: whether that random slope is effectively 0
*	varcov.rand.eigen.2 near 0: same idea, the second random-effects dimension is collapsing

~0.22% of CpGs have effectively zero random‐slope variance
* there is no between‐subject heterogeneity in age slope


# check data matrix (X) for redundancy - DOES NOT EVALUATE Z
mv <- ~ age_yrs + gender + cc + CaseControl + new_Bcell + new_CD4T + new_CD8T + new_Mono + new_NK
X <- model.matrix(mv, data = pheno_complete)

qr(X)$rank
ncol(X)
kappa(X)


* Full rank -> no exact collinearity?
* How close the matrix is to being singular
* Large singular value -> strong direction in predictor space
* Small singular value -> nearly redundant / weakly identifiable direction

It does not say anything about:
* random-effects variance structure
* boundary variance estimates
* random slope identifiability
* covariance singularities

This is why kappa looks good, and we still get convergence warnings. They come from the random-effects (Z), not from X


# check that whole-blood samples have variability
covs <- c("new_Bcell","new_CD4T","new_CD8T","new_Mono","new_NK")
sds <- sapply(pheno_filt_ct_unscaled[, covs], sd, na.rm=TRUE)
which(sds == 0 | is.na(sds))
colSums(!is.finite(as.matrix(pheno_filt_ct_unscaled[, covs])))


vars <- c("sample_agedys","new_Bcell","new_CD4T","new_CD8T","new_Mono","new_NK")
sapply(pheno_filt_ct_unscaled[, vars], function(x) sum(!is.finite(x)))

sapply(pheno_filt_ct_unscaled[, c("new_Bcell","new_CD4T","new_CD8T","new_Mono","new_NK")],
       function(x) sd(x, na.rm = TRUE))





### Identify CpGs that did not converge


# identify CpGs not fit in Step1B model
step1b_bad_cpgs <- step1b_og_results %>% filter(!is.na(message))


# Step1B with scaled variables - CpGs that did not fit
step1b_scaled_bad_cpgs <- step1b_scaled_results %>% filter(!is.na(message))


# lme model CpGs that did not fit (my script)
lme_bad_cpgs <- lme_results %>% filter(!is.na(error_message))



# lmeR model CpGs that did not fit (my script)
## the lmeR provides estimates for every CpG but comes w/ convergence warnings
lmeR_bad_cpgs <- lmeR_results %>% filter(!is.na(warn_message))

# lmeR model, using bobyqa optimizer, CpGs that did not fit
## the lmeR provides estimates for every CpG but comes w/ convergence warnings
lmeR_bobyqa_bad_cpgs <- lmeR_bobyqa_results %>% filter(!is.na(warn_message))

# how many of the "bad" CpGs overlap between the models
cat("Step1B, lme, lmeR, and lmeR w/bobyqa had", 
    length(intersect(step1b_bad_cpgs$CpG, lme_bad_cpgs$CpG, 
                     lmeR_bad_cpgs, lmeR_bobyqa_bad_cpgs)), 
    "overlapping CpGs that failed\n")

# unique CpGs to all models
unique_step1b <- setdiff( # finds all rows in x that aren't in y : setdiff(x,y)
  step1b_bad_cpgs$CpG,
  lme_bad_cpgs$CpG)

cat("Step1B had", length(unique_step1b), "unique CpGs that failed\n")

unique_lme <- setdiff(
  lme_bad_cpgs$CpG,
  step1b_bad_cpgs$CpG)

cat("LME had", length(unique_lme), "unique CpGs that failed\n")



Number of singular fits using lmeR

cat("lmeR model w/ bobyqa had",
    nrow(lmeR_MAD_results[lmeR_MAD_results$is_singular==TRUE,]),
    "singular fits\n")

cat("lmeR model w/splines had",
    nrow(lmeR_spline_results[lmeR_spline_results$is_singular==TRUE,]),
    "singular fits\n")

cat("lmeR model w/splines-int had",
    nrow(lmeR_spline_int_res[lmeR_spline_int_res$is_singular==TRUE,]),
    "singular fits\n")

cat(
  sprintf("%.3f%%\n", 
          nrow(lmeR_MAD_results[lmeR_MAD_results$is_singular==TRUE,]) / 
            nrow(lmeR_MAD_results) * 100))

cat(
  sprintf("%.3f%%\n", 
          nrow(lmeR_spline_results[lmeR_spline_results$is_singular==TRUE,]) / 
            nrow(lmeR_spline_results) * 100))

cat(
  sprintf("%.3f%%\n", 
          nrow(lmeR_spline_int_res[lmeR_spline_int_res$is_singular==TRUE,]) /
            nrow(lmeR_spline_int_res) * 100))


Compare number of significant slopes across methods before multiple-testing adjustment

# lmeR w/ bobyqa optimization and outliers removed
cat("lmeR model w/bobyqa optimization & MAD outliers removed had", 
    nrow(subset(lmeR_MAD_results, p_slope < 0.05)), 
    "significant slopes representing", 
    100*round(nrow(subset(lmeR_MAD_results, p_slope < 0.05)) / 
                nrow(lmeR_MAD_results),3), "% of total CpG tested\n")


* Model w/BOBYQA and MAD outliers removed retained 14450 more CpGs

Identify CpGs that had a significant slope but was singular - var(slope)≈0

lmeR_bobyqa_results <- lmeR_bobyqa_results %>% 
  mutate(sig_singular = ifelse(
    is_singular == TRUE & p_slope<0.05, "Trouble", "Good"))

lmeR_results <- lmeR_results %>% 
  mutate(sig_singular = ifelse(
    is_singular == TRUE & p_slope<0.05, "Trouble", "Good"))

lmeR_MAD_results <- lmeR_MAD_results %>% 
  mutate(sig_singular = ifelse(
    is_singular == TRUE & p_slope<0.05, "Trouble", "Good"))


# number of CpGs with a significant slope but was singular
cat("lmeR w/default optimization had",
    nrow(lmeR_results[lmeR_results$sig_singular=="Trouble",]),
    "CpGs with a significant slope but was singular\n", 
    "(subjects are inferred to have same slope)\n\n")

cat("lmeR w/bobyqa optimization had",
    nrow(lmeR_bobyqa_results[lmeR_bobyqa_results$sig_singular=="Trouble",]),
    "CpGs with a significant slope but was singular\n", 
    "(subjects are inferred to have same slope)\n\n")

cat("lmeR w/bobyqa optimization and MAD outliers removed had ",
    nrow(lmeR_bobyqa_results[lmeR_MAD_results$sig_singular=="Trouble",]),
    "CpGs with a significant slope but was singular\n", 
    "(subjects are inferred to have same slope)")


* 9474 more CpGs with significant slopes and singular

Using the bobyqa optimization and increasing maxfun let the optimizer converge more fully, and for many CpGs the REML estimate for the random-slope variance is on the boundary (≈0), so more fits are flagged as singular. This suggests the random slope is not needed (no heterogeneity) for those CpGs, and AIC-based screening should be used with caution and looked at explicitly.


"Model failed to converge with max|grad| = 0.00399344 (tol = 0.002, component 1)":
  
  That warning means lmer returned estimates, but the optimizer stopped with a gradient that’s a bit larger than the default tolerance. In plain terms:
  
  * the algorithm got close to the optimum, but not close enough to declare full convergence.

So the results aren’t automatically “wrong,” but you should treat that CpG’s fit as potentially unreliable (especially SEs / variance components).


Eigenvalue: 
  0.0345 - Substantial random-effects variability
5×10⁻⁹- Essentially zero variability

Interpretation:
  There is strong subject-to-subject variability in one direction (intercept), and almost none in the other (slope-related direction).

* Eigenvalues > 0 -> matrix is positive definite (good)
* Very small eigenvalues -> one direction has almost no variance
*	≈ 0 eigenvalue -> model is near-singular

