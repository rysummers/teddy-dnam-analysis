## Mediation

library(mediation)
mediators <- c("new_Bcell", "new_CD4T", "new_CD8T", "new_Mono", "new_NK")

results <- lapply(mediators, function(med) {
  
  med_form <- as.formula(
    sprintf("%s ~ age_yrs + gender + cc + CaseControl", med))
  
  out_form <- as.formula(
    sprintf("CpG ~ age_yrs + %s + gender + cc + CaseControl", med))
  
  m_med <- lm(med_form, data = dat)
  m_out <- lm(out_form, data = dat)
  
  mediate(
    model.m = m_med,
    model.y = m_out,
    treat = "age_yrs",
    mediator = med,
    sims = 500)
})

# put mediation results in a table
summary_table <- lapply(seq_along(results), function(i) {
  res <- summary(results[[i]])
  data.frame(
    mediator = mediators[i],
    ACME = res$d0,
    ACME_p = res$d0.p,
    ADE = res$z0,
    ADE_p = res$z0.p,
    total = res$tau.coef,
    total_p = res$tau.p,
    prop_med = res$n0)
}) %>%
  bind_rows() %>% 
  mutate(
    across(c(ACME, ADE, total), ~ round(.x, 4)),
    across(c(prop_med), ~ round(.x, 3)),
    across(ends_with("_p"), ~ round(.x, 10)))

summary_table


ACME = mediated effect
ADE = what remains after adjusting for mediator
Total effect = w/out mediator : age --> CpG
Prop. Med =  % of the age effect explained by the mediator

In single-mediator analyses, the association between age and CpG methylation was
partially mediated by estimated blood cell composition, with the strongest mediation 
observed for CD4T cells (ACME = 0.0575, p < 0.001; proportion mediated = 0.61) and
B cells (ACME = 0.0433, p < 0.001; proportion mediated = 0.46). CD8T showed a 
smaller but significant mediated effect, whereas Mono and NK showed little or no 
evidence of mediation. Because the direct effect of age remained significant in 
all models, the results support partial rather than complete mediation
