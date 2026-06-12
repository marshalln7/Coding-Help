############################################################
# Propensity Score Matching (ATT) Tutorial: LaLonde / NSW
# Stat 437 — MatchIt workflow
#
# Goal:
#   Estimate the ATT of job training (treat) on earnings in 1978 (re78)
#   using propensity score matching.
#
############################################################

# ---- 0) Packages ----
library(MatchIt)
library(cobalt)    # balance diagnostics (love plot, SMDs)
library(tidyverse)
library(marginaleffects)

set.cobalt.options(binary = "std",continuous = "std")

# ---- 1) Load data ----
# MatchIt includes the classic LaLonde data.
data("lalonde")
str(lalonde)
# Key variables:
#   treat : 1 = participated in job training, 0 = control
#   re78  : earnings in 1978 (outcome)
#   covariates: age, educ, race/ethnicity indicators, married, nodegree, re74, re75

# ---- 2) Define the propensity score model ----
# We model P(treat=1 | X). This is the treatment assignment model.
# Important principle:
#   We include PRE-TREATMENT covariates that predict treatment and/or outcome.

ps_formula <- treat ~ age + educ + race + married + nodegree + re74 + re75

# ---- 3) Run propensity score matching (targets ATT by default) ----
# MatchIt defaults (method = "nearest") do nearest-neighbor matching on the PS.
#
# KNOB 1: ratio = number of controls matched to each treated unit (1 = 1:1)
#
# KNOB 2: caliper
#   caliper restricts matches to be "close enough" in PS distance.
#   Smaller caliper -> better matches but may drop treated units with no good match.
#   Common choice: 0.2 * SD(logit(PS)) (MatchIt can interpret caliper on PS scale;
#   when distance="logit", caliper is on logit scale.)
#
# KNOB 3: replace
#   replace = TRUE allows controls to be reused (can improve match quality when controls are scarce).
#
# KNOB 4: distance
#   distance = "logit" uses logistic regression and matches on logit(PS).
#
# KNOB 5: estimand - ATE normally requires full matching
#   estimand = "ATT" (default for nearest neighbor in MatchIt) targets treated population.
#
m_out <- matchit(
  formula  = ps_formula,
  data     = lalonde,
  method   = "nearest",
  distance = "glm",    # KNOB 4: try "gbm"
  estimand = "ATT",      # KNOB 5: try ATE -- you'll need to change method = "
  ratio    = 2,          # KNOB 1: try 2 or 3
  caliper  = .2,        # KNOB 2: try 0.1 or 0.3; set to NULL for no caliper
  replace  = FALSE       # KNOB 3: try TRUE
)

# MatchIt summary: how many matched, how many dropped, balance summaries
summary(m_out)


# ---- 4) Check covariate balance ----
# We want covariates balanced AFTER matching.
#
# "love.plot" visualizes standardized mean differences (SMDs).
# Rule of thumb: |SMD| < 0.1 or 0.25 is “good” balance (not a law, but a useful target).

love.plot(
  m_out,
  drop.distance = TRUE,
  stats = c("mean.diffs"),
  thresholds = c(m = c(0.25)), ### could change to 0.1
  abs = TRUE
)

# If balance is not good:
# - try adding nonlinear terms/interactions to the PS model,
# - tighten caliper,
# - increase ratio (sometimes helps; sometimes hurts),
# - or allow replacement.

# ---- 5) Extract the matched data ----
# match.data() returns a dataset with weights and matched sample indicators.
m_dat <- match.data(m_out)

# The returned data include:
# - weights: matching weights for each observation
# - distance: estimated propensity score (or logit PS depending on setting)
# - treat, re78, covariates

# ---- 6) Estimate ATT in the matched sample ----
# For 1:1 nearest-neighbor matching without replacement, a simple weighted difference
# in means is a straightforward ATT estimator.
#
# ATT = E[Y | treat=1, matched] - E[Y | treat=0, matched] with matching weights.

att_est <- with(
  m_dat,
  weighted.mean(re78[treat == 1], weights[treat == 1]) -
    weighted.mean(re78[treat == 0], weights[treat == 0])
)
att_est

# ---- 7) A slightly more model-based ATT (optional, still simple) ----
# You can fit an outcome regression in the matched sample using weights.
# This can help with residual imbalance and improves precision sometimes.
#
# This is NOT "required" for matching, but it's a common practical step.
fit_att <- lm(re78 ~ treat, 
              data = m_dat, weights = weights)
summary(fit_att)

avg_comparisons(fit_att,
                variables = "treat",
                vcov = ~subclass,
                newdata = subset(treat == 1))



# The coefficient on treat is an ATT-style estimate in this matched sample.

fit_att2 <- lm(re78 ~ treat * (age + educ + race + married + nodegree + re74 + re75), 
              data = m_dat, weights = weights)
summary(fit_att2)

avg_comparisons(fit_att2,
                variables = "treat",
                vcov = ~subclass,
                newdata = subset(treat == 1))

