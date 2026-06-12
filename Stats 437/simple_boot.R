# -------------------------------------------------------------------------
# GOAL: Estimate the Standard Error of the Ratio of Means using Bootstrap
# Estimator: (Mean Bill Length) / (Mean Bill Depth)
# -------------------------------------------------------------------------

# 1. Setup and Data Cleaning
library(palmerpenguins)
library(tidyverse)

# We'll use Gentoo penguins as our "population" sample
df <- na.omit(penguins[penguins$species == "Gentoo", ])

# 2. Define the Estimator Function
# This takes the data and a set of indices to calculate the ratio
calc_ratio_of_means <- function(data_subset) {
  mean_length <- mean(data_subset$bill_length_mm)
  mean_depth  <- mean(data_subset$bill_depth_mm)
  return(mean_length / mean_depth)
}

# 3. Calculate the "Observed" Point Estimate
# This is our best guess using the original data
observed_ratio <- calc_ratio_of_means(df)

# 4. The Bootstrap Loop (Manual Implementation)
set.seed(03312026)              # For reproducibility
B <- 2000            # Number of bootstrap resamples
n_obs  <- nrow(df)        # Number of observations in our original sample
boot_dist <- numeric(B) # Pre-allocate a vector to store results

for (i in 1:B) {
  
  # STEP A: Resample WITH replacement
  # This is the "magic". We create a new dataset of the same size
  # by drawing from the original. Some rows appear multiple times; some zero.
  idx_boot <- sample(1:n_obs, n_obs, replace = TRUE)
  df_boot <- df[idx_boot,]

  
  # STEP B: Apply the estimator to the resampled data
  
  boot_dist[i] <- calc_ratio_of_means(df_boot)
}

# 5. Quantifying Uncertainty
# The Bootstrap Standard Error is just the standard deviation of our results
boot_se <- sd(boot_dist)

boot_se

# 95% Percentile Confidence Interval
# We look at the 2.5th and 97.5th percentiles of our simulated ratios
conf_int <- quantile(boot_dist, probs = c(0.025, 0.975))
conf_int
