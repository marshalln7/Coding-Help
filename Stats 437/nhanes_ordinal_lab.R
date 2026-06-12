# nhanes_ordinal_lab.R
#
# Goals:
# - Fit and interpret cumulative logit models (proportional odds)
# - Explore nonlinearity and interactions
# - Compare nested & non-nested models (LRT, AIC, BIC)
# - Check proportional odds assumption
# - Visualize predicted probabilities and compute simple CV log-loss
#
# Run top-to-bottom; read comments. Sections clearly marked.

# ---------------------------
# 0) Packages and setup
# ---------------------------
library(tidyverse)
library(NHANES)    # dataset
library(ordinal)   # clm()
library(broom)     # tidy()
library(patchwork) # plot layout
library(scales)    # percent_format()
library(splines)
set.seed(03102026)

# ---------------------------
# 1) Load & prepare data
# ---------------------------
data("NHANES")  # loads tibble NHANES

wanted <- c("HealthGen", "BMI", "Age", "Diabetes", "PhysActive",
            "SmokeNow", "Gender", "BPDiaAve", "Education","Poverty",
            "MaritalStatus")

# Select what is available and coerce types
NHANES_work <- NHANES %>%
  dplyr::select(all_of(wanted)) %>% 
  mutate(
    HealthGen = factor(HealthGen,
                       levels = c("Poor","Fair","Good","Vgood","Excellent"),
                       ordered = TRUE)
  ) %>% 
  filter(complete.cases(.))

NHANES_work$MaritalStatus = fct_relevel(NHANES_work$MaritalStatus,"Married")
NHANES_work$Education = fct_relevel(NHANES_work$Education,"8th Grade")

# Quick peek
glimpse(NHANES_work)
table(NHANES_work$HealthGen)

# ---------------------------
# 2) Exploratory plots
# ---------------------------
# Continuous boxplots (HealthGen on x-axis)
p_bmi <- ggplot(NHANES_work, aes(x = HealthGen, y = BMI)) +
  geom_boxplot(fill = "cornflowerblue", outlier.size = 0.6) +
  labs(title = "BMI by Self-Reported Health", x = NULL, y = "BMI") +
  theme(axis.text.x = element_text(angle = 18, hjust = 1))

p_age <- ggplot(NHANES_work, aes(x = HealthGen, y = Age)) +
  geom_boxplot(fill = "lightseagreen", outlier.size = 0.6) +
  labs(title = "Age by Self-Reported Health", x = NULL, y = "Age (years)") +
  theme(axis.text.x = element_text(angle = 18, hjust = 1))

p_pov <- ggplot(NHANES_work, aes(x = HealthGen, y = Poverty)) +
  geom_boxplot(fill = "plum", outlier.size = 0.6) +
  labs(title = "Poverty Index by Self-Reported Health", x = NULL, y = "Poverty Ratio/Score") +
  theme(axis.text.x = element_text(angle = 18, hjust = 1))

(p_bmi + p_age) / p_pov  + plot_layout(heights = c(1, 1))

# Categorical: proportion within each HealthGen (stacked bar normalized)
p_diab <- ggplot(NHANES_work, aes(x = HealthGen, fill = Diabetes)) +
  geom_bar(position = "fill") +
  scale_y_continuous(labels = scales::percent_format()) +
  labs(title = "HealthGen by Diabetes (proportion within HealthGen)", y = "Proportion", x = NULL)

p_phys <- ggplot(NHANES_work, aes(x = HealthGen, fill = PhysActive)) +
  geom_bar(position = "fill") +
  scale_y_continuous(labels = scales::percent_format()) +
  labs(title = "HealthGen by Physical Activity (proportion within HealthGen)", y = "Proportion", x = NULL)

(p_diab + p_phys) + plot_layout(ncol = 2)

# Quick numeric summaries to guide model building
NHANES_work %>% group_by(HealthGen) %>% summarise(n = n(), meanBMI = mean(BMI, na.rm = TRUE), medianAge = median(Age, na.rm=TRUE))

# ---------------------------
# 3) Baseline cumulative logit model
# ---------------------------
fit_clm_base <- clm(HealthGen ~ BMI + Age + Diabetes + PhysActive + Poverty,
                    data = NHANES_work,
                    link = "logit", Hess = TRUE)

summary(fit_clm_base)
tidy(fit_clm_base, conf.int = TRUE)  # coefficients + thresholds

# Interpretation prompt (for students):
# - Which coefficients are negative/positive? What does sign mean for ordered outcome?
# - Compute exp(beta) for a couple of variables:
exp_coef <- exp(coef(fit_clm_base))
exp_coef

# EXERCISE: Interpret at least one categorical variable and one continuous.


# ---------------------------
# 4) Predicted probabilities for profiles
# ---------------------------
# Helper to get predicted probs in tidy format

get_pred_df <- function(fit, newdata) {
  # predict(..., type = "prob") usually returns a list with $fit (matrix) or matrix
  pr <- predict(fit, newdata = newdata, type = "prob")
  if (is.list(pr) && !is.null(pr$fit)) {
    m <- pr$fit
  } else {
    m <- pr
  }
  m <- as.data.frame(m)
  colnames(m) <- colnames(m) # ensure names exist
  out <- cbind(newdata, m) %>%
    pivot_longer(cols = colnames(m), names_to = "category", values_to = "prob")
  out
}

# Example: plot predicted probabilities vs BMI when other covariates fixed
newdat <- expand.grid(
  BMI = seq(18, 40, by = 0.5),
  Age = median(NHANES_work$Age, na.rm = TRUE),
  Diabetes = "Yes",
  PhysActive = "Yes",  # pick first level
  Poverty = 1
)
pred_df <- get_pred_df(fit_clm_base, newdat)

ggplot(pred_df, aes(x = BMI, y = prob, color = category)) +
  geom_line(size = 1) +
  labs(title = "Predicted HealthGen probabilities vs BMI (baseline CLM)",
       y = "Predicted probability", x = "BMI") +
  theme_minimal(base_size = 14)

# EXERCISE: Change variables (e.g., Diabetes = "No") and replot; discuss differences.

# ---------------------------
# 5) Check proportional odds assumption
# ---------------------------
# nominal_test() from ordinal package checks nominal effects (a type of test)
nom_test_res <- nominal_test(fit_clm_base)
print(nom_test_res)

# ---------------------------
# 6) Nonlinearity exploration for BMI
# ---------------------------
# Fit polynomial and spline variants and compare (AIC / LRT)
fit_poly2 <- clm(HealthGen ~ poly(BMI, 2) + Age + Diabetes + PhysActive + Poverty, data = NHANES_work, link="logit", Hess=TRUE)
fit_ns4  <- clm(HealthGen ~ ns(BMI, df = 4) + Age + Diabetes + PhysActive + Poverty, data = NHANES_work, link="logit", Hess=TRUE)

# Compare AIC / BIC
AIC(fit_clm_base, fit_poly2, fit_ns4)
BIC(fit_clm_base, fit_poly2, fit_ns4)

# LRT (nested?)
anova(fit_clm_base, fit_poly2)   
anova(fit_clm_base, fit_ns4)

# Visualize predicted probability curves for each model
newdat2 <- newdat
pred_base <- get_pred_df(fit_clm_base, newdat2) %>% mutate(model = "base")
pred_poly <- get_pred_df(fit_poly2, newdat2) %>% mutate(model = "poly2")
pred_ns   <- get_pred_df(fit_ns4, newdat2) %>% mutate(model = "ns4")

pred_all <- bind_rows(pred_base, pred_poly, pred_ns)

ggplot(pred_all %>% filter(category == "Good"), aes(x = BMI, y = prob, color = model)) +
  geom_line() +
  labs(title = "Probability of 'Good' vs BMI: model comparisons", y = "Pr(Good)")

# EXERCISE: Pick a category (Poor or Excellent) and compare shapes. Which model seems to capture curvature?
# EXERCISE: Pick a different variable to explore nonlinearity and repeat process

# ---------------------------
# 7) Interactions
# ---------------------------
# Example interaction: BMI * PhysActive
fit_int1 <- clm(HealthGen ~ BMI * PhysActive + Age + Diabetes + Poverty, data = NHANES_work, link = "logit", Hess = TRUE)
summary(fit_int1)

# Test whether interaction improves fit
anova(fit_clm_base, fit_int1)  # LRT (nested) — check p-value

# Visualize interaction: predicted probs by BMI for levels of PhysActive
newdat_int <- expand.grid(
  BMI = seq(18, 40, by = 0.5),
  Age = median(NHANES_work$Age, na.rm = TRUE),
  Diabetes = "Yes",
  PhysActive = levels(NHANES_work$PhysActive),
  Poverty = median(NHANES_work$Poverty, na.rm = TRUE)
)
pred_int_df <- get_pred_df(fit_int1, newdat_int)

ggplot(pred_int_df %>% filter(category == "Poor"), aes(x = BMI, y = prob, color = PhysActive)) +
  geom_line() +
  labs(title = "Predicted probability of 'Poor' by BMI and PhysActive (interaction)",
       y = "Pr(Poor)")

# EXERCISE: Try different interaction. Which are meaningful?

# ---------------------------
# 8) More complex model & comparisons
# ---------------------------
fit_complex <- clm(HealthGen ~ BMI + Age + Diabetes + PhysActive + Poverty + SmokeNow + MaritalStatus + Education,
                   data = NHANES_work, link = "logit", Hess = TRUE)

AIC(fit_clm_base, fit_complex)
BIC(fit_clm_base, fit_complex)
anova(fit_clm_base, fit_complex)  # LRT
