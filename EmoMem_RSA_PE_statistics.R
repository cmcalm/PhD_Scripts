library(readxl)
library(lme4)
library(lmerTest)
library(dplyr)
library(writexl)
library(lavaan)
library(car)
library(ggplot2)
library(broom.mixed)
library(ggeffects)


## RSA
Gamma_df <- read_excel("//cbsu/data/Group/Camcan/sandbox/cc07/Project3_IAPS_RSA/Results/Gammas_RSA.xlsx")


# De-differentiation analysis: fit γ ~ Valence * poly(Age,2) + (1 | Subject)
Gamma_df$ppt <- as.factor(Gamma_df$ppt)
Gamma_df$val <- as.factor(Gamma_df$val)
Gamma_df$thres <- as.factor(Gamma_df$thres)

Gamma_u <- Gamma_df %>% filter(thres == "u")
Gamma_t <- Gamma_df %>% filter(thres == "t")

sum(Gamma_u$gamma < 0)
sum(Gamma_t$gamma < 0)
model_u <- lmer(gamma ~ val * poly(Age, 2) + (1 | ppt), data = Gamma_u)
model_t <- lmer(gamma ~ val * poly(Age, 2) + (1 | ppt), data = Gamma_t)

# assumption checks: linearity, normality of residuals, homoscedasticity of residuals
plot(model_u)
qqnorm(residuals(model_u))
qqline(residuals(model_u))
plot(fitted(model_u), residuals(model_u))

plot(model_t)
qqnorm(residuals(model_t))
qqline(residuals(model_t))
plot(fitted(model_t), residuals(model_t))

summary(model_u)
summary(model_t)

ggplot(Gamma_t, aes(x = Age, y = gamma, colour = val)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  labs(
    x = "Age",
    y = "Gamma",
    colour = "Valence"
  ) +
  theme_minimal() +
  labs(title = "Age-Related Change in Neural Differentiation by Valence") + 
  theme(axis.title = element_text(size = 12), axis.text  = element_text(size = 10),
        legend.title = element_text(size = 12), legend.text  = element_text(size = 10))

ggplot(Gamma_u, aes(x = Age, y = gamma, colour = val)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  labs(
    x = "Age",
    y = "Gamma",
    colour = "Valence"
  ) +
  theme_minimal() +
  labs(title = "Age-Related Change in Neural Differentiation by Valence (Unthresholded)") + 
  theme(axis.title = element_text(size = 12), axis.text  = element_text(size = 10),
        legend.title = element_text(size = 12), legend.text  = element_text(size = 10))

coef_u <- summary(model_u)$coefficients
coef_t <- summary(model_t)$coefficients

confint(model_u, method = "Wald")
confint(model_t, method = "Wald")


coef_u_df <- as.data.frame(coef_u)
coef_u_df$term <- rownames(coef_u_df)
coef_u_df$thres <- "u"

coef_t_df <- as.data.frame(coef_t)
coef_t_df$term <- rownames(coef_t_df)
coef_t_df$thres <- "t"

ci_u <- as.data.frame(confint(model_u, method = "Wald"))
ci_u$term <- rownames(ci_u)
ci_u$thres <- "u"

ci_t <- as.data.frame(confint(model_t, method = "Wald"))
ci_t$term <- rownames(ci_t)
ci_t$thres <- "t"

coef_u_full <- left_join(coef_u_df, ci_u, by = c("term", "thres"))
coef_t_full <- left_join(coef_t_df, ci_t, by = c("term", "thres"))

coef_full <- bind_rows(coef_u_full, coef_t_full)

write_xlsx(coef_full, "//cbsu/data/Group/Camcan/sandbox/cc07/Project3_IAPS_RSA/Results/Dediff_LME_result.xlsx")


## Moderated Mediation
# Fit 4 Models: A simple mediation model without age as a moderator, a model including age as a moderator on the a-path 
# (representational similarity to connectivity), on the b-path (connectivity to memory), and on both parts of the indirect effect
# Compare models with significant moderator terms using a Chi-square test
# conditional indirect effects (moderated mediation) will be quantified by examining the effect at the mean age +/- 1 SD
# the p-value for the coefficient of the interaction will be Holm Bonferroni-corrected (4 model comparisons)

Data_for_ModMed <- read_excel("//cbsu/data/Group/Camcan/sandbox/cc07/Project3_IAPS_RSA/Results/Data_for_ModMed.xlsx")

# plots for assumption checks:
plot(Data_for_ModMed$beta_t, Data_for_ModMed$dCor_Amy_u) # linearity (looks ok)
plot(Data_for_ModMed$beta_t, Data_for_ModMed$dCor_Hip_u)
plot(Data_for_ModMed$beta_t, Data_for_ModMed$dCor_dmPFC_t)
plot(Data_for_ModMed$beta_t, Data_for_ModMed$dCor_IPL_t)
plot(Data_for_ModMed$beta_t_Age, Data_for_ModMed$dCor_Amy_u)
plot(Data_for_ModMed$beta_t_Age, Data_for_ModMed$dCor_Hip_u)
plot(Data_for_ModMed$beta_t_Age, Data_for_ModMed$dCor_dmPFC_t)
plot(Data_for_ModMed$beta_t_Age, Data_for_ModMed$dCor_IPL_t)
plot(Data_for_ModMed$dCor_Amy_u, Data_for_ModMed$Rec_Diff)
plot(Data_for_ModMed$dCor_Hip_u, Data_for_ModMed$Rec_Diff)
plot(Data_for_ModMed$dCor_dmPFC_t, Data_for_ModMed$Rec_Diff)
plot(Data_for_ModMed$dCor_IPL_t, Data_for_ModMed$Rec_Diff)
plot(Data_for_ModMed$dCor_Amy_u_Age, Data_for_ModMed$Rec_Diff)
plot(Data_for_ModMed$dCor_Hip_u_Age, Data_for_ModMed$Rec_Diff)
plot(Data_for_ModMed$dCor_dmPFC_t_Age, Data_for_ModMed$Rec_Diff)
plot(Data_for_ModMed$dCor_IPL_t_Age, Data_for_ModMed$Rec_Diff)

subset_data <- subset(Data_for_ModMed[, 3:ncol(Data_for_ModMed)], select = -c(Age_c, Age_c2)) 
boxplot(subset_data[, sapply(subset_data, is.numeric)]) # a lot of outliers for the interactions


# model specifications:
model0_template <- "
  # Mediator
  {mediator} ~ a1*{beta}

  # Outcome
  Rec_Diff ~ b1*{mediator} + c1*{beta}
  
  ind := a1*b1
"

modela_template <- "
  {mediator} ~ a1*{beta} + a2*{beta}_Age + a3*Age

  Rec_Diff ~ b1*{mediator} + c1*{beta} + c2*{beta}_Age + c3*Age
  
  ind_mean := (a1 + a2*0)*b1
  ind_low  := (a1 + a2*(-1))*b1
  ind_high := (a1 + a2*(1))*b1
"

modelb_template <- "
  {mediator} ~ a1*{beta}

  Rec_Diff ~ b1*{mediator} + b2*{mediator}_Age + c1*{beta} + c2*{beta}_Age + c3*Age
  
  ind_mean := a1*(b1 + b2*0)
  ind_low  := a1*(b1 + b2*(-1))
  ind_high := a1*(b1 + b2*(1))
"

modelab_template <- "
  {mediator} ~ a1*{beta} + a2*{beta}_Age + a3*Age

  Rec_Diff ~ b1*{mediator} + b2*{mediator}_Age + c1*{beta} + c2*{beta}_Age + c3*Age
  

  ind_mean := (a1 + a2*0)*(b1 + b2*0)
  ind_low  := (a1 + a2*(-1))*(b1 + b2*(-1))
  ind_high := (a1 + a2*(1))*(b1 + b2*(1))
  
"

mediators <- list(
  dCor_Amy_u = "beta_u",
  dCor_Hip_u = "beta_u",
  dCor_IPL_t = "beta_t",
  dCor_IPL_u = "beta_u",
  dCor_dmPFC_t = "beta_t",
  dCor_dmPFC_u = "beta_u"
)

build_model <- function(template, mediator, beta) {
  model <- gsub("\\{mediator\\}", mediator, template)
  model <- gsub("\\{beta\\}", beta, model)
  return(model)
}


ModMed_results <- list()

for (med in names(mediators)) {
  beta <- mediators[[med]]
  m0  <- build_model(model0_template, med, beta)
  ma  <- build_model(modela_template, med, beta)
  mb  <- build_model(modelb_template, med, beta)
  mab <- build_model(modelab_template, med, beta)
  fit0  <- sem(m0,  data = Data_for_ModMed, meanstructure = TRUE, se = "bootstrap", bootstrap = 5000)
  fita  <- sem(ma,  data = Data_for_ModMed, meanstructure = TRUE, se = "bootstrap", bootstrap = 5000)
  fitb  <- sem(mb,  data = Data_for_ModMed, meanstructure = TRUE, se = "bootstrap", bootstrap = 5000)
  fitab <- sem(mab, data = Data_for_ModMed, meanstructure = TRUE, se = "bootstrap", bootstrap = 5000)
  ModMed_results[[med]] <- list(fit0 = fit0, fita = fita, fitb = fitb, fitab = fitab)
}


extract_params <- function(fit, model_name, mediator) {
  pe <- parameterEstimates(fit, ci = TRUE, standardized = TRUE)
  pe$Param <- paste0(pe$lhs, "_", pe$op, "_", pe$rhs)
  pe$model <- model_name
  pe$mediator <- mediator
  
  pe <- pe %>%
    dplyr::filter(op %in% c("~", ":=")) %>%
    dplyr::mutate(
      is_interaction = grepl("_Age", Param)
    ) %>%
    dplyr::select(mediator, model, Param, est, std.all, se, z, pvalue, ci.lower, ci.upper, is_interaction)
  
  return(pe)
}

param_results <- dplyr::bind_rows(
  lapply(names(ModMed_results), function(med) {
    
    fits <- ModMed_results[[med]]
    
    dplyr::bind_rows(
      extract_params(fits$fit0,  "model0", med),
      extract_params(fits$fita,  "modela", med),
      extract_params(fits$fitb,  "modelb", med),
      extract_params(fits$fitab, "modelab", med)
    )
    
  })
)

param_results <- param_results %>%
  dplyr::mutate(
    thres = ifelse(grepl("_t$", mediator), "t",
                   ifelse(grepl("_u$", mediator), "u", NA))
  ) %>%
  dplyr::group_by(thres, Param) %>%
  dplyr::mutate(
    p_holm = p.adjust(pvalue, method = "holm")
  ) %>%
  dplyr::ungroup()


extract_fitmeasures <- function(fit, model_name, mediator) {
  fm <- fitMeasures(fit, c("chisq", "df", "pvalue", "cfi", "rmsea", "aic", "bic"))
  data.frame(
    mediator = mediator,
    model = model_name,
    chisq = fm["chisq"],
    df = fm["df"],
    pvalue = fm["pvalue"],
    cfi = fm["cfi"],
    rmsea = fm["rmsea"],
    aic = fm["aic"],
    bic = fm["bic"]
  )
}

fitmeasures_results <- dplyr::bind_rows(
  lapply(names(ModMed_results), function(med) {
    
    fits <- ModMed_results[[med]]
    
    dplyr::bind_rows(
      extract_fitmeasures(fits$fit0,  "model0", med),
      extract_fitmeasures(fits$fita,  "modela", med),
      extract_fitmeasures(fits$fitb,  "modelb", med),
      extract_fitmeasures(fits$fitab, "modelab", med)
    )
    
  })
)

# For assumption checks, inspect fit indices: 
# we want CFI ≥ .95 (good), ≥ .90 (acceptable); RMSEA ≤ .06 (good), ≤ .08 (acceptable); χ²: non-sig.

write_xlsx(fitmeasures_results, "//cbsu/data/Group/Camcan/sandbox/cc07/Project3_IAPS_RSA/Results/ModMed_fit.xlsx")
write_xlsx(param_results, "//cbsu/data/Group/Camcan/sandbox/cc07/Project3_IAPS_RSA/Results/ModMed_params.xlsx")


# Subsequent Memory Analysis
# implement a binomial generalised linear mixed model, predicting the binary variable of whether an item is 
# later correctly recalled (gist or detail memory) based on the trial-level RSA result for positive and negative 
# backgrounds. 

Trial_RSA <- read_excel("//cbsu/data/Group/Camcan/sandbox/cc07/Project3_IAPS_RSA/Results/trialRSA_results.xlsx")
Trial_RSA$Recall <- ifelse(Trial_RSA$Recall == "Remembered", 1, 0)
tRSA_df <- filter(Trial_RSA, Valence %in% c('Positive','Negative'))
tRSA_df$Valence <- factor(tRSA_df$Valence, levels = c("Negative", "Positive"))
tRSA_df$Age_c <- scale(tRSA_df$Age, center = TRUE, scale = FALSE)
tRSA_df$rsa_t_c <- scale(tRSA_df$rsa_t, center = TRUE, scale = FALSE)
tRSA_df$rsa_u_c <- scale(tRSA_df$rsa_u, center = TRUE, scale = FALSE)
tRSA_df$participant <- factor(tRSA_df$participant)
tRSA_df$stimulus <- factor(tRSA_df$stimulus)

sm_mod_t <- glmer(Recall ~ rsa_t_c  * Age_c * Valence + (1 | participant) + (1 | stimulus) , data = tRSA_df, family = binomial)
summary(sm_mod_t)

sm_mod_u <- glmer(Recall ~ rsa_u_c  * Age_c * Valence + (1 | participant) + (1 | stimulus) , data = tRSA_df, family = binomial)
summary(sm_mod_u)

# assumption checks
vif(lm(rsa_t ~ Age_c * Valence, data = tRSA_df)) # some collinearity, but below concerning threshold
vif(lm(rsa_u ~ Age_c * Valence, data = tRSA_df))

ggplot(tRSA_df, aes(rsa_t, Recall)) +
  geom_smooth(method = "glm", method.args = list(family = "binomial")) # looks like linear relationship for both
ggplot(tRSA_df, aes(rsa_u, Recall)) +
  geom_smooth(method = "glm", method.args = list(family = "binomial"))

overdisp_fun <- function(model) {
  rdf <- df.residual(model)
  rp <- residuals(model, type = "pearson")
  sqrt(sum(rp^2) / rdf)
}

overdisp_fun(sm_mod_t) # no overdispersion
overdisp_fun(sm_mod_u)


extract_glmer_results <- function(model, model_name) {
  res <- broom.mixed::tidy(model, effects = "fixed", conf.int = TRUE)
  
  res <- res %>%
    mutate(
      model = model_name,
      odds_ratio = exp(estimate),    # result is in log odds, transform to get odds ratio
      OR_low = exp(conf.low),    
      OR_high = exp(conf.high)
    )
  
  return(res)
}

res_t <- extract_glmer_results(sm_mod_t, "rsa_t") %>%
  mutate(thres = "t")

res_u <- extract_glmer_results(sm_mod_u, "rsa_u") %>%
  mutate(thres = "u")

glmm_results <- bind_rows(res_t, res_u)

write_xlsx(glmm_results,
           "//cbsu/data/Group/Camcan/sandbox/cc07/Project3_IAPS_RSA/Results/SubMem_GLMM_results.xlsx"
)


pred_t <- ggpredict(sm_mod_t, terms = c("rsa_t_c", "Age_c [-1,0,1]", "Valence")) %>%
  as.data.frame() %>%
  mutate(thres = "t")

pred_u <- ggpredict(sm_mod_u, terms = c("rsa_u_c", "Age_c [-1,0,1]", "Valence")) %>%
  as.data.frame() %>%
  mutate(thres = "u")

pred_comb <- bind_rows(pred_t, pred_u)

write_xlsx(as.data.frame(pred_comb),
           "//cbsu/data/Group/Camcan/sandbox/cc07/Project3_IAPS_RSA/Results/SubMem_pred.xlsx"
)


ggplot(pred_t, aes(x = x, y = predicted, colour = group)) +
  geom_line(size = 1.2) +
  facet_wrap(~facet) +
  labs(
    x = "Trial-level RSA (centered)",
    y = "Predicted Recall Probability",
    colour = "Age",
    title = "Subsequent Memory Effect"
  ) +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5),
        strip.text = element_text(size = 12),
        strip.background = element_blank())

ggplot(pred_u, aes(x = x, y = predicted, colour = group)) +
  geom_line(size = 1.2) +
  facet_wrap(~facet) +
  labs(
    x = "Trial-level RSA (centered)",
    y = "Predicted Recall Probability",
    colour = "Age",
    title = "Subsequent Memory Effect (unthresholded)" 
  ) +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5),
        strip.text = element_text(size = 12),
        strip.background = element_blank())


# repeat analyses with arousal-level specific data

## RSA
Gamma_df <- read_excel("//cbsu/data/Group/Camcan/sandbox/cc07/Project3_IAPS_RSA/Results/Gammas_RSA_aro_split.xlsx")


# De-differentiation analysis: fit γ ~ Valence * poly(Age,2) + (1 | Subject)
Gamma_df$ppt <- as.factor(Gamma_df$ppt)
Gamma_df$val <- as.factor(Gamma_df$val)
Gamma_df$thres <- as.factor(Gamma_df$thres)

Gamma_low <- Gamma_df %>% filter(arousal == "low")
Gamma_high <- Gamma_df %>% filter(arousal == "high")
Gamma_low_u <- Gamma_low %>% filter(thres == "u")
Gamma_low_t <- Gamma_low %>% filter(thres == "t")
Gamma_high_u <- Gamma_high %>% filter(thres == "u")
Gamma_high_t <- Gamma_high %>% filter(thres == "t")

sum(Gamma_low_t$gamma < 0)
sum(Gamma_low_u$gamma < 0)
sum(Gamma_high_t$gamma < 0)
sum(Gamma_high_u$gamma < 0) # negative values for gamma are likely to reflect noise, should we clip them to 0?

model_u_low <- lmer(gamma ~ val * poly(Age, 2) + (1 | ppt), data = Gamma_low_u)
model_t_low <- lmer(gamma ~ val * poly(Age, 2) + (1 | ppt), data = Gamma_low_t)
model_u_high <- lmer(gamma ~ val * poly(Age, 2) + (1 | ppt), data = Gamma_high_u)
model_t_high <- lmer(gamma ~ val * poly(Age, 2) + (1 | ppt), data = Gamma_high_t)

# assumption checks: linearity, normality of residuals, homoscedasticity of residuals
plot(model_u_low)
qqnorm(residuals(model_u_low))
qqline(residuals(model_u_low))
plot(fitted(model_u_low), residuals(model_u_low))

plot(model_t_low)
qqnorm(residuals(model_t_low))
qqline(residuals(model_t_low))
plot(fitted(model_t_low), residuals(model_t_low))

plot(model_u_high)
qqnorm(residuals(model_u_high))
qqline(residuals(model_u_high))
plot(fitted(model_u_high), residuals(model_u_high))

plot(model_t_high)
qqnorm(residuals(model_t_high))
qqline(residuals(model_t_high))
plot(fitted(model_t_high), residuals(model_t_high))
# several violations, do we need to address them? 

summary(model_u_low)
summary(model_t_low)
summary(model_u_high)
summary(model_t_high)

ggplot(Gamma_low_t, aes(x = Age, y = gamma, colour = val)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  labs(
    x = "Age",
    y = "Gamma",
    colour = "Valence"
  ) +
  theme_minimal() +
  labs(title = "Age-Related Change in Neural Differentiation by Valence, Low Arousal") + 
  theme(axis.title = element_text(size = 12), axis.text  = element_text(size = 10),
        legend.title = element_text(size = 12), legend.text  = element_text(size = 10))

ggplot(Gamma_high_t, aes(x = Age, y = gamma, colour = val)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  labs(
    x = "Age",
    y = "Gamma",
    colour = "val"
  ) +
  theme_minimal() +
  labs(title = "Age-Related Change in Neural Differentiation by Valence, High Arousal") + 
  theme(axis.title = element_text(size = 12), axis.text  = element_text(size = 10),
        legend.title = element_text(size = 12), legend.text  = element_text(size = 10))

coef_u_low <- summary(model_u_low)$coefficients
coef_t_low <- summary(model_t_low)$coefficients
coef_u_high <- summary(model_u_high)$coefficients
coef_t_high <- summary(model_t_high)$coefficients

coef_u_low_df <- as.data.frame(coef_u_low)
coef_u_low_df$term <- rownames(coef_u_low_df)
coef_u_low_df$thres <- "u"

coef_t_low_df <- as.data.frame(coef_t_low)
coef_t_low_df$term <- rownames(coef_t_low_df)
coef_t_low_df$thres <- "t"

coef_u_high_df <- as.data.frame(coef_u_high)
coef_u_high_df$term <- rownames(coef_u_high_df)
coef_u_high_df$thres <- "u"

coef_t_high_df <- as.data.frame(coef_t_high)
coef_t_high_df$term <- rownames(coef_t_high_df)
coef_t_high_df$thres <- "t"

ci_u_low <- as.data.frame(confint(model_u_low, method = "Wald"))
ci_u_low$term <- rownames(ci_u_low)
ci_u_low$thres <- "u"
ci_u_low$arousal <- "low"

ci_t_low <- as.data.frame(confint(model_t_low, method = "Wald"))
ci_t_low$term <- rownames(ci_t_low)
ci_t_low$thres <- "t"
ci_t_low$arousal <- "low"

ci_u_high <- as.data.frame(confint(model_u_high, method = "Wald"))
ci_u_high$term <- rownames(ci_u_high)
ci_u_high$thres <- "u"
ci_u_high$arousal <- "high"

ci_t_high <- as.data.frame(confint(model_t_high, method = "Wald"))
ci_t_high$term <- rownames(ci_t_high)
ci_t_high$thres <- "t"
ci_t_high$arousal <- "high"

coef_u_full_low <- left_join(coef_u_low_df, ci_u_low, by = c("term", "thres"))
coef_t_full_low <- left_join(coef_t_low_df, ci_t_low, by = c("term", "thres"))
coef_u_full_high <- left_join(coef_u_high_df, ci_u_high, by = c("term", "thres"))
coef_t_full_high <- left_join(coef_t_high_df, ci_t_high, by = c("term", "thres"))

coef_full_low <- bind_rows(coef_u_full_low, coef_t_full_low, coef_u_full_high, coef_t_full_high)

write_xlsx(coef_full_low, "//cbsu/data/Group/Camcan/sandbox/cc07/Project3_IAPS_RSA/Results/Dediff_LME_aro_split_result.xlsx")


## Moderated Mediation
Data_for_ModMed_s <- read_excel("//cbsu/data/Group/Camcan/sandbox/cc07/Project3_IAPS_RSA/Results/Data_for_ModMed_aro_split.xlsx")
Data_for_ModMed_low <- filter(Data_for_ModMed_s, arousal == 'low')
Data_for_ModMed_high <- filter(Data_for_ModMed_s, arousal == 'high')

# plots for assumption checks:
plot(Data_for_ModMed_low$beta_t, Data_for_ModMed_low$dCor_Hip_u)
plot(Data_for_ModMed_low$beta_t, Data_for_ModMed_low$dCor_dmPFC_t)
plot(Data_for_ModMed_low$beta_t, Data_for_ModMed_low$dCor_IPL_t)
plot(Data_for_ModMed_low$beta_t_Age, Data_for_ModMed_low$dCor_Amy_u)
plot(Data_for_ModMed_low$beta_t_Age, Data_for_ModMed_low$dCor_Hip_u)
plot(Data_for_ModMed_low$beta_t_Age, Data_for_ModMed_low$dCor_dmPFC_t)
plot(Data_for_ModMed_low$beta_t_Age, Data_for_ModMed_low$dCor_IPL_t)
plot(Data_for_ModMed_low$dCor_Amy_u, Data_for_ModMed_low$Rec_Diff)
plot(Data_for_ModMed_low$dCor_Hip_u, Data_for_ModMed_low$Rec_Diff)
plot(Data_for_ModMed_low$dCor_dmPFC_t, Data_for_ModMed_low$Rec_Diff)
plot(Data_for_ModMed_low$dCor_IPL_t, Data_for_ModMed_low$Rec_Diff)
plot(Data_for_ModMed_low$dCor_Amy_u_Age, Data_for_ModMed_low$Rec_Diff)
plot(Data_for_ModMed_low$dCor_Hip_u_Age, Data_for_ModMed_low$Rec_Diff)
plot(Data_for_ModMed_low$dCor_dmPFC_t_Age, Data_for_ModMed_low$Rec_Diff)
plot(Data_for_ModMed_low$dCor_IPL_t_Age, Data_for_ModMed_low$Rec_Diff)

subset_data <- subset(Data_for_ModMed_low[, 3:ncol(Data_for_ModMed_low)], select = -c(Age_c, Age_c2)) 
boxplot(subset_data[, sapply(subset_data, is.numeric)])


ModMed_results_low <- list()
for (med in names(mediators)) {
  beta <- mediators[[med]]
  m0  <- build_model(model0_template, med, beta)
  ma  <- build_model(modela_template, med, beta)
  mb  <- build_model(modelb_template, med, beta)
  mab <- build_model(modelab_template, med, beta)
  fit0  <- sem(m0,  data = Data_for_ModMed_low, meanstructure = TRUE, se = "bootstrap", bootstrap = 5000)
  fita  <- sem(ma,  data = Data_for_ModMed_low, meanstructure = TRUE, se = "bootstrap", bootstrap = 5000)
  fitb  <- sem(mb,  data = Data_for_ModMed_low, meanstructure = TRUE, se = "bootstrap", bootstrap = 5000)
  fitab <- sem(mab, data = Data_for_ModMed_low, meanstructure = TRUE, se = "bootstrap", bootstrap = 5000)
  ModMed_results_low[[med]] <- list(fit0 = fit0, fita = fita, fitb = fitb, fitab = fitab)
}

ModMed_results_high <- list()
for (med in names(mediators)) {
  beta <- mediators[[med]]
  m0  <- build_model(model0_template, med, beta)
  ma  <- build_model(modela_template, med, beta)
  mb  <- build_model(modelb_template, med, beta)
  mab <- build_model(modelab_template, med, beta)
  fit0  <- sem(m0,  data = Data_for_ModMed_high, meanstructure = TRUE, se = "bootstrap", bootstrap = 5000)
  fita  <- sem(ma,  data = Data_for_ModMed_high, meanstructure = TRUE, se = "bootstrap", bootstrap = 5000)
  fitb  <- sem(mb,  data = Data_for_ModMed_high, meanstructure = TRUE, se = "bootstrap", bootstrap = 5000)
  fitab <- sem(mab, data = Data_for_ModMed_high, meanstructure = TRUE, se = "bootstrap", bootstrap = 5000)
  ModMed_results_high[[med]] <- list(fit0 = fit0, fita = fita, fitb = fitb, fitab = fitab)
}


param_results_low <- dplyr::bind_rows(
  lapply(names(ModMed_results_low), function(med) {
    
    fits <- ModMed_results_low[[med]]
    
    dplyr::bind_rows(
      extract_params(fits$fit0,  "model0", med),
      extract_params(fits$fita,  "modela", med),
      extract_params(fits$fitb,  "modelb", med),
      extract_params(fits$fitab, "modelab", med)
    )
    
  })
)

param_results_low <- param_results_low %>%
  dplyr::mutate(
    thres = ifelse(grepl("_t$", mediator), "t",
                   ifelse(grepl("_u$", mediator), "u", NA))
  ) %>%
  dplyr::group_by(thres, Param) %>%
  dplyr::mutate(
    p_holm = p.adjust(pvalue, method = "holm")
  ) %>%
  dplyr::ungroup()


param_results_high <- dplyr::bind_rows(
  lapply(names(ModMed_results_high), function(med) {
    
    fits <- ModMed_results_high[[med]]
    
    dplyr::bind_rows(
      extract_params(fits$fit0,  "model0", med),
      extract_params(fits$fita,  "modela", med),
      extract_params(fits$fitb,  "modelb", med),
      extract_params(fits$fitab, "modelab", med)
    )
    
  })
)

param_results_high <- param_results_high %>%
  dplyr::mutate(
    thres = ifelse(grepl("_t$", mediator), "t",
                   ifelse(grepl("_u$", mediator), "u", NA))
  ) %>%
  dplyr::group_by(thres, Param) %>%
  dplyr::mutate(
    p_holm = p.adjust(pvalue, method = "holm")
  ) %>%
  dplyr::ungroup()



fitmeasures_results_low <- dplyr::bind_rows(
  lapply(names(ModMed_results_low), function(med) {
    
    fits <- ModMed_results_low[[med]]
    
    dplyr::bind_rows(
      extract_fitmeasures(fits$fit0,  "model0", med),
      extract_fitmeasures(fits$fita,  "modela", med),
      extract_fitmeasures(fits$fitb,  "modelb", med),
      extract_fitmeasures(fits$fitab, "modelab", med)
    )
    
  })
)

fitmeasures_results_high <- dplyr::bind_rows(
  lapply(names(ModMed_results_high), function(med) {
    
    fits <- ModMed_results_high[[med]]
    
    dplyr::bind_rows(
      extract_fitmeasures(fits$fit0,  "model0", med),
      extract_fitmeasures(fits$fita,  "modela", med),
      extract_fitmeasures(fits$fitb,  "modelb", med),
      extract_fitmeasures(fits$fitab, "modelab", med)
    )
    
  })
)

write_xlsx(fitmeasures_results_low, "//cbsu/data/Group/Camcan/sandbox/cc07/Project3_IAPS_RSA/Results/ModMed_fit_low.xlsx")
write_xlsx(param_results_low, "//cbsu/data/Group/Camcan/sandbox/cc07/Project3_IAPS_RSA/Results/ModMed_params_low.xlsx")
write_xlsx(fitmeasures_results_high, "//cbsu/data/Group/Camcan/sandbox/cc07/Project3_IAPS_RSA/Results/ModMed_fit_high.xlsx")
write_xlsx(param_results_high, "//cbsu/data/Group/Camcan/sandbox/cc07/Project3_IAPS_RSA/Results/ModMed_params_high.xlsx")



# Subsequent Memory Analysis
Trial_RSA <- read_excel("//cbsu/data/Group/Camcan/sandbox/cc07/Project3_IAPS_RSA/Results/trialRSA_results.xlsx")
Trial_RSA$Recall <- ifelse(Trial_RSA$Recall == "Remembered", 1, 0)
tRSA_df <- filter(Trial_RSA, Valence %in% c('Positive','Negative'))
tRSA_df$Valence <- factor(tRSA_df$Valence, levels = c("Negative", "Positive"))
tRSA_df$Age_c <- scale(tRSA_df$Age, center = TRUE, scale = FALSE)
tRSA_df$rsa_t_c <- scale(tRSA_df$rsa_t, center = TRUE, scale = FALSE)
tRSA_df$rsa_u_c <- scale(tRSA_df$rsa_u, center = TRUE, scale = FALSE)
tRSA_df$participant <- factor(tRSA_df$participant)
tRSA_df$stimulus <- factor(tRSA_df$stimulus)
tRSA_df_low <- tRSA_df %>% filter(aro == "low")
tRSA_df_high <- tRSA_df %>% filter(aro == "high")

aro_df <- read_excel("//cbsu/data/Group/Camcan/sandbox/cc07/Project3_IAPS_RSA/IAPS_ratings/Libkuman_2007_Ratings.xls")
aro_df <- aro_df %>%
  select(`1APS#`, Libkuman_arousal_M, Libkuman_val_M) %>%
  rename(stimulus = `1APS#`)
aro_df$stimulus <- factor(aro_df$stimulus) 
new_df <- tRSA_df %>%
  left_join(
    aro_df %>% select(stimulus, Libkuman_arousal_M),
    by = "stimulus"
  )
new_df$Libkuman_arousal_M <- scale(new_df$Libkuman_arousal_M, center = TRUE, scale = FALSE)

sm_mod_t_int <- glmer(Recall ~ rsa_t_c  * Age_c * Valence * Libkuman_arousal_M + (1 | participant) + (1 | stimulus) , data = new_df, family = binomial)
summary(sm_mod_t_int)


sm_mod_t_low <- glmer(Recall ~ rsa_t_c  * Age_c * Valence + (1 | participant) + (1 | stimulus) , data = tRSA_df_low, family = binomial)
summary(sm_mod_t_low)

sm_mod_u_low <- glmer(Recall ~ rsa_u_c  * Age_c * Valence + (1 | participant) + (1 | stimulus) , data = tRSA_df_low, family = binomial)
summary(sm_mod_u_low)

sm_mod_t_high <- glmer(Recall ~ rsa_t_c  * Age_c * Valence + (1 | participant) + (1 | stimulus) , data = tRSA_df_high, family = binomial)
summary(sm_mod_t_high)

sm_mod_u_high <- glmer(Recall ~ rsa_u_c  * Age_c * Valence + (1 | participant) + (1 | stimulus) , data = tRSA_df_high, family = binomial)
summary(sm_mod_u_high)

# assumption checks
vif(lm(rsa_t ~ Age_c * Valence, data = tRSA_df_low))
vif(lm(rsa_u ~ Age_c * Valence, data = tRSA_df_low))
vif(lm(rsa_t ~ Age_c * Valence, data = tRSA_df_high))
vif(lm(rsa_u ~ Age_c * Valence, data = tRSA_df_high))

ggplot(tRSA_df_low, aes(rsa_t, Recall)) +
  geom_smooth(method = "glm", method.args = list(family = "binomial"))
ggplot(tRSA_df_low, aes(rsa_u, Recall)) +
  geom_smooth(method = "glm", method.args = list(family = "binomial"))
ggplot(tRSA_df_high, aes(rsa_t, Recall)) +
  geom_smooth(method = "glm", method.args = list(family = "binomial"))
ggplot(tRSA_df_high, aes(rsa_u, Recall)) +
  geom_smooth(method = "glm", method.args = list(family = "binomial"))

overdisp_fun(sm_mod_t_low)
overdisp_fun(sm_mod_u_low)
overdisp_fun(sm_mod_t_high)
overdisp_fun(sm_mod_u_high)


res_t_low <- extract_glmer_results(sm_mod_t_low, "rsa_t") %>%
  mutate(thres = "t", aro = "low")

res_u_low <- extract_glmer_results(sm_mod_u_low, "rsa_u") %>%
  mutate(thres = "u", aro = "low")

res_t_high <- extract_glmer_results(sm_mod_t_high, "rsa_t") %>%
  mutate(thres = "t", aro = "high")

res_u_high <- extract_glmer_results(sm_mod_u_high, "rsa_u") %>%
  mutate(thres = "u", aro = "high")

glmm_results_aro <- bind_rows(res_t_low, res_u_low, res_t_high, res_u_high)

write_xlsx(glmm_results_aro,
           "//cbsu/data/Group/Camcan/sandbox/cc07/Project3_IAPS_RSA/Results/SubMem_GLMM_results_aroSplit.xlsx"
)


pred_t_low <- ggpredict(sm_mod_t_low, terms = c("rsa_t_c", "Age_c [-1,0,1]", "Valence")) %>%
  as.data.frame() %>%
  mutate(thres = "t", aro = "low")

pred_u_low <- ggpredict(sm_mod_u_low, terms = c("rsa_u_c", "Age_c [-1,0,1]", "Valence")) %>%
  as.data.frame() %>%
  mutate(thres = "u", aro = "low")

pred_t_high <- ggpredict(sm_mod_t_high, terms = c("rsa_t_c", "Age_c [-1,0,1]", "Valence")) %>%
  as.data.frame() %>%
  mutate(thres = "t", aro = "high")

pred_u_high <- ggpredict(sm_mod_u_high, terms = c("rsa_u_c", "Age_c [-1,0,1]", "Valence")) %>%
  as.data.frame() %>%
  mutate(thres = "u", aro = "high")

pred_comb <- bind_rows(pred_t_low, pred_u_low, pred_t_high, pred_u_high)

write_xlsx(as.data.frame(pred_comb),
           "//cbsu/data/Group/Camcan/sandbox/cc07/Project3_IAPS_RSA/Results/SubMem_pred_aroSplit.xlsx"
)


ggplot(pred_t_low, aes(x = x, y = predicted, colour = group)) +
  geom_line(size = 1.2) +
  facet_wrap(~facet) +
  labs(
    x = "Trial-level RSA (centered)",
    y = "Predicted Recall Probability",
    colour = "Age",
    title = "Subsequent Memory Effect (Low Arousal)"
  ) +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5),
    strip.text = element_text(size = 12),
    strip.background = element_blank())


ggplot(pred_t_high, aes(x = x, y = predicted, colour = group)) +
  geom_line(size = 1.2) +
  facet_wrap(~facet) +
  labs(
    x = "Trial-level RSA (centered)",
    y = "Predicted Recall Probability",
    colour = "Age",
    title = "Subsequent Memory Effect (High Arousal)"
  ) +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5),
        strip.text = element_text(size = 12),
        strip.background = element_blank())

ggplot(pred_u_low, aes(x = x, y = predicted, colour = group)) +
  geom_line(size = 1.2) +
  facet_wrap(~facet) +
  labs(
    x = "Trial-level RSA (centered)",
    y = "Predicted Recall Probability",
    colour = "Age",
    title = "Subsequent Memory Effect (Unthresholded, Low Arousal)"
  ) +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5),
        strip.text = element_text(size = 12),
        strip.background = element_blank())

ggplot(pred_u_high, aes(x = x, y = predicted, colour = group)) +
  geom_line(size = 1.2) +
  facet_wrap(~facet) +
  labs(
    x = "Trial-level RSA (centered)",
    y = "Predicted Recall Probability",
    colour = "Age",
    title = "Subsequent Memory Effect (Unthresholded, High Arousal)"
  ) +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5),
        strip.text = element_text(size = 12),
        strip.background = element_blank())

ggplot(tRSA_df_low, aes(x = rsa_t_c, y = Recall, colour = Valence)) +
  geom_smooth(method = "glm", method.args = list(family = "binomial")) +
  theme_classic() +
  labs(
    y = "Probability of Recall",
    x = "RSA (centered)",
    title = "RSA Predicting Memory by Valence, Low Arousal"
  )

ggplot(tRSA_df_high, aes(x = rsa_t_c, y = Recall, colour = Valence)) +
  geom_smooth(method = "glm", method.args = list(family = "binomial")) +
  theme_classic() +
  labs(
    y = "Probability of Recall",
    x = "RSA (centered)",
    title = "RSA Predicting Memory by Valence, High Arousal"
  )



# Exploratory Moderated Mediation with Categorical Representations:
Data_for_ModMed_cat <- read_excel("//cbsu/data/Group/Camcan/sandbox/cc07/Project3_IAPS_RSA/Results/Data_for_ModMed_Cat.xlsx")

mediators_cat <- list(
  beta_amy_u = "beta_u",
  beta_hip_u = "beta_u",
  beta_IPL_t = "beta_t",
  beta_IPL_u = "beta_u",
  beta_dmPFC_t = "beta_t",
  beta_dmPFC_u = "beta_u"
)


ModMed_cat_results <- list()

for (med in names(mediators_cat)) {
  beta <- mediators_cat[[med]]
  m0  <- build_model(model0_template, med, beta)
  ma  <- build_model(modela_template, med, beta)
  mb  <- build_model(modelb_template, med, beta)
  mab <- build_model(modelab_template, med, beta)
  fit0  <- sem(m0,  data = Data_for_ModMed_cat, meanstructure = TRUE, se = "bootstrap", bootstrap = 5000)
  fita  <- sem(ma,  data = Data_for_ModMed_cat, meanstructure = TRUE, se = "bootstrap", bootstrap = 5000)
  fitb  <- sem(mb,  data = Data_for_ModMed_cat, meanstructure = TRUE, se = "bootstrap", bootstrap = 5000)
  fitab <- sem(mab, data = Data_for_ModMed_cat, meanstructure = TRUE, se = "bootstrap", bootstrap = 5000)
  ModMed_cat_results[[med]] <- list(fit0 = fit0, fita = fita, fitb = fitb, fitab = fitab)
}


param_cat_results <- dplyr::bind_rows(
  lapply(names(ModMed_cat_results), function(med) {
    
    fits <- ModMed_cat_results[[med]]
    
    dplyr::bind_rows(
      extract_params(fits$fit0,  "model0", med),
      extract_params(fits$fita,  "modela", med),
      extract_params(fits$fitb,  "modelb", med),
      extract_params(fits$fitab, "modelab", med)
    )
    
  })
)

param_cat_results <- param_cat_results %>%
  dplyr::mutate(
    thres = ifelse(grepl("_t$", mediator), "t",
                   ifelse(grepl("_u$", mediator), "u", NA))
  ) %>%
  dplyr::group_by(thres, Param) %>%
  dplyr::mutate(
    p_holm = p.adjust(pvalue, method = "holm")
  ) %>%
  dplyr::ungroup()


fitmeasures_cat_results <- dplyr::bind_rows(
  lapply(names(ModMed_cat_results), function(med) {
    
    fits <- ModMed_cat_results[[med]]
    
    dplyr::bind_rows(
      extract_fitmeasures(fits$fit0,  "model0", med),
      extract_fitmeasures(fits$fita,  "modela", med),
      extract_fitmeasures(fits$fitb,  "modelb", med),
      extract_fitmeasures(fits$fitab, "modelab", med)
    )
    
  })
)

# For assumption checks, inspect fit indices: 
# we want CFI ≥ .95 (good), ≥ .90 (acceptable); RMSEA ≤ .06 (good), ≤ .08 (acceptable); χ²: non-sig.

write_xlsx(fitmeasures_cat_results, "//cbsu/data/Group/Camcan/sandbox/cc07/Project3_IAPS_RSA/Results/ModMed_fit_cat.xlsx")
write_xlsx(param_cat_results, "//cbsu/data/Group/Camcan/sandbox/cc07/Project3_IAPS_RSA/Results/ModMed_params_cat.xlsx")



