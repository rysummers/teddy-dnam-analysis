### Generalized Additive Model (GAM)


gam.lin <- gam(
  CpG ~ age_yrs + gender + cc + CaseControl +
    new_Bcell + new_CD4T + new_CD8T + new_Mono + new_NK,
  data = dat,
  method = "ML")

gam.smooth <- gam(
  CpG ~ s(age_yrs, bs = "cr") + gender + cc + CaseControl +
    new_Bcell + new_CD4T + new_CD8T + new_Mono + new_NK,
  data = dat,
  method = "ML")

gam.re <- gam(
  CpG ~ s(age_yrs, bs = "cr") +
    gender + cc + CaseControl +
    new_Bcell + new_CD4T + new_CD8T + new_Mono + new_NK +
    s(maskid, bs = "re") +
    s(maskid, by = age_yrs, bs = "re"),
  data = dat,
  method = "ML")

summary(gam.lin)
summary(gam.smooth)
summary(gam.re)

anova(gam.lin, gam.smooth, gam.re,
      test = "Chisq")

AIC(gam.lin, gam.smooth,gam.re)



gam.check(gam.smooth)




plot(gam.smooth)
