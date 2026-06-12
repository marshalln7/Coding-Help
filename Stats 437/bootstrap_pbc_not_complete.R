# == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == 
# Bootstrap and Survival: Cox model with PBC data
# == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == 
# This script mirrors the slide deck:
#   1) Fit a Cox PH model
#   2) Estimate baseline cumulative hazard
#   3) Turn hazard into survival for a covariate profile
#   4) Compute a treatment-minus-placebo survival difference
#   5) Bootstrap the difference for SE / CI
#   6) Build a null resampling distribution for a p-value
# == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == 

library(tidyverse)
library(survival)
library(survminer)


# ------------------------------------------------------------
# Helper functions - Thank you chatGPT!!!
# ------------------------------------------------------------

# Get the step-function value of the baseline cumulative hazard at time t0
get_H0_at_time = function(basehaz_df, t0) {
  idx = findInterval(t0, basehaz_df$time)
  if (idx ==  0) return(0)
  basehaz_df$hazard[idx]
}

# Helper: survival at a single time point from a survfit object
survival_at_time <- function(sf, t0) {
  s <- summary(sf, times = t0, extend = TRUE)
  as.numeric(s$surv)
}

# Helper: full survival curve from a survfit object
survival_curve <- function(sf, group_label) {
  tibble(
    time = sf$time,
    survival = sf$surv,
    group = group_label
  )
}

# Helper: treatment-minus-placebo survival difference at one time point
survival_difference_at_time <- function(sf_tx, sf_pbo, t0) {
  s_tx <- survival_at_time(sf_tx, t0)
  s_pbo <- survival_at_time(sf_pbo, t0)
  s_tx - s_pbo
}


set.seed(03262026)

# ------------------------------------------------------------
# User-controlled settings -- feel free to change these!
# ------------------------------------------------------------
time0 = 365*7   # seven years in days
B = 1000      # Bootstrap replicates (increase if desired)

# ------------------------------------------------------------
# 1) Load and prepare the PBC data
# ------------------------------------------------------------
# In pbc:
#   status ==  2 means death
#   status ==  1 is transplant
# Here, we treat death and treatment as an event. Everything else is censored.
data(pbc)

pbc_use = pbc %>%
  transmute(
    time =  time,
    event =  ifelse(status == 0, 0, 1), # death (1) or transplant (2) =  1, otherwise censored
    trt =  factor(trt, levels =  c(1, 2), labels =  c("Placebo", "Treatment")),
    age =  age,
    bili =  bili,
    albumin =  albumin
    # edema =  edema,
    # hepato =  hepato
  ) %>%
  drop_na()

# Center continuous covariates at their sample medians.
# After centering, 0 corresponds to the "median patient".
pbc_use = pbc_use %>%
  mutate(
    age_c =  age - median(age),
    log_bili_c =  log(bili) - median(log(bili)),
    albumin_c =  albumin - median(albumin),
    # edema =  factor(edema),
    # hepato =  factor(hepato)
  )

# Quick look at the analysis dataset
print(glimpse(pbc_use))


km_trt = survfit(Surv(time, event) ~ trt, data =  pbc_use)

ggsurvplot(
  km_trt,
  data =  pbc_use,
  conf.int =  TRUE,
  palette =  c("firebrick", "steelblue"),
  xlab =  "Days",
  ylab =  "Recurrence-free survival probability",
  legend.title =  "",
  legend.labs =  c("Placebo", "Treatment"),
  ggtheme =  theme_minimal(base_size =  14)
)

# ------------------------------------------------------------
# 2) Fit the Cox proportional hazards model
# ------------------------------------------------------------

cox_formula = Surv(time, event) ~ trt * (age_c + log_bili_c + albumin_c)
cox_fit = coxph(cox_formula, data =  pbc_use)

# Model output
summary(cox_fit)

# ------------------------------------------------------------
# 3) Define a reference patient profile
# ------------------------------------------------------------
# Because the continuous covariates are centered, 0 means "median".
# This gives a clean profile for comparing Treatment vs Placebo.
profile_pbo = tibble(
  trt =  factor("Placebo", levels =  levels(pbc_use$trt)),
  age_c =  0,
  log_bili_c =  0,
  albumin_c =  0
)

profile_tx = tibble(
  trt =  factor("Treatment", levels =  levels(pbc_use$trt)),
  age_c =  0,
  log_bili_c =  0,
  albumin_c =  0
)

# ------------------------------------------------------------
# 4) Estimate the baseline cumulative hazard
# ------------------------------------------------------------
# Cox does not model h0(t), but we can estimate the cumulative hazard
# after fitting the model.
sf_pbo <- survfit(cox_fit, newdata = profile_pbo)
sf_tx  <- survfit(cox_fit, newdata = profile_tx)

# 7-year survival difference for the chosen profile
delta_obs <- survival_difference_at_time(sf_tx, sf_pbo, time0)

results_point <- tibble(
  time_days = time0,
  survival_placebo = survival_at_time(sf_pbo, time0),
  survival_treatment = survival_at_time(sf_tx, time0),
  delta_survival = delta_obs
)

print(results_point)

# Optional plot: predicted survival curves
curve_pbo <- survival_curve(sf_pbo, "Placebo")
curve_tx  <- survival_curve(sf_tx, "Treatment")

curve_df <- bind_rows(curve_pbo, curve_tx) %>%
  mutate(time = time / 365)

time0_here <- time0 / 365

ggplot(curve_df, aes(x = time, y = survival, color = group)) +
  geom_step(linewidth = 1.1) +
  geom_vline(xintercept = time0_here, linetype = "dashed") +
  labs(
    title = "Predicted survival curves from the Cox model",
    subtitle = "Median patient profile: centered covariates = 0",
    x = "Years",
    y = "Survival probability",
    color = ""
  ) +
  theme_minimal(base_size = 14) +
  theme(legend.position = "bottom")

diff_df <- curve_pbo %>%
  dplyr::select(time, survival_pbo = survival) %>%
  inner_join(curve_tx %>% dplyr::select(time, survival_tx = survival), by = "time") %>%
  mutate(delta = survival_tx - survival_pbo)

ggplot(diff_df, aes(x = time / 365, y = delta)) +
  geom_step(linewidth = 1.1) +
  geom_hline(yintercept = 0, linetype = "dotted") +
  geom_vline(xintercept = time0_here, linetype = "dashed") +
  labs(
    title = "Treatment effect on survival over time",
    subtitle = "Treatment minus placebo",
    x = "Years",
    y = expression(S[Treatment] - S[Placebo])
  ) +
  theme_minimal(base_size = 14)



# ------------------------------------------------------------
# 5) Nonparametric bootstrap for SE and percentile CI
# ------------------------------------------------------------
# Bootstrap idea:
#   Resample rows with replacement from the observed dataset,
#   refit the Cox model, and recompute the 7-year survival difference.
# Write a loop or function that will get the bootstrap null distribution.
# Return delta_boot








boot_se <- sd(delta_boot)
boot_ci <- quantile(delta_boot, probs = c(0.025, 0.975), names = FALSE)

results_boot <- tibble(
  delta_hat = delta_obs,
  boot_se = boot_se,
  ci_lower = boot_ci[1],
  ci_upper = boot_ci[2],
  n_boot = length(delta_boot)
)

print(results_boot)

tibble(delta = delta_boot) %>% 
  ggplot(aes(x = delta)) +
  geom_histogram(bins = 30, color = "white") +
  geom_vline(xintercept = delta_obs, linetype = "dashed", linewidth = 1) +
  geom_vline(xintercept = boot_ci, linetype = "dotted", linewidth = 1) +
  labs(
    title = "Bootstrap distribution of the 7-year treatment effect",
    subtitle = "Treatment minus placebo survival difference",
    x = expression(hat(Delta)(7*365, X)),
    y = "Count"
  ) +
  theme_minimal(base_size = 14)

# ------------------------------------------------------------
# 6) Null resampling distribution for a p-value
# ------------------------------------------------------------
# For a randomized treatment comparison, a simple null resampling scheme is:
#   - pool the data
#   - resample subjects with replacement
#   - assign them to treatment and placebo groups with the original group sizes
# Write a loop or function that will get the bootstrap null distribution.
# Return delta_null



# One-sided p-value for the alternative that treatment improves survival
p_one_sided <- mean(delta_null >= delta_obs)

# Two-sided p-value from the null distribution
p_two_sided <- mean(abs(delta_null) >= abs(delta_obs))

results_p <- tibble(
  delta_hat = delta_obs,
  p_one_sided = p_one_sided,
  p_two_sided = p_two_sided,
  n_null = length(delta_null)
)

print(results_p)

# Optional plot: null resampling distribution
delta_null_df <- tibble(delta = delta_null)

ggplot(delta_null_df, aes(x = delta)) +
  geom_histogram(bins = 30, color = "white") +
  geom_vline(xintercept = delta_obs, linetype = "dashed", linewidth = 1) +
  labs(
    title = "Null resampling distribution of the 7-year treatment effect",
    subtitle = "Used for bootstrap-style p-value calculation",
    x = expression(hat(Delta)[0](7*365, X)),
    y = "Count"
  ) +
  theme_minimal(base_size = 14)

# ------------------------------------------------------------
# 7) Compact summary table for reporting - Thanks chatGPT!!!
# ------------------------------------------------------------

summary_table <- tibble(
  quantity = c(
    "7-year survival: placebo",
    "7-year survival: treatment",
    "Treatment effect at 7 years",
    "Bootstrap SE",
    "95% bootstrap CI lower",
    "95% bootstrap CI upper",
    "One-sided null p-value",
    "Two-sided null p-value"
  ),
  value = c(
    survival_at_time(sf_pbo, time0),
    survival_at_time(sf_tx, time0),
    delta_obs,
    boot_se,
    boot_ci[1],
    boot_ci[2],
    p_one_sided,
    p_two_sided
  )
)

print(summary_table)
