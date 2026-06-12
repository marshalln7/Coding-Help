# Multicollinearity ILLUSTRATION ------------------------------------------
# This is for illustration purposes—you do not need to understand 
# all of this code in detail! The goal is to demonstrate the impact 
# of multicollinearity on standard errors in regression.

library(tidyverse)
library(mvtnorm) # to generate correlated predictors (X variables)

# True Regression Coefficients --------------------------------------------
# We assume the true model follows:
# Y = β0 + β1*X1 + β2*X2 + β3*X3 + ε
b_vec <- c(4,3,-2,1)  # True regression coefficients (including intercept)
sig <- 2  # Standard deviation of the error term

n <- 50     # Sample size
n_samp = 1000  # Number of simulations
set.seed(1)    # For reproducibility

# Simulating Data with Uncorrelated Predictors ----------------------------

# Storage matrices for estimated coefficients and their standard errors
MLR <- matrix(NA, n_samp, length(b_vec) - 1)  
std_err_beta <- matrix(NA, n_samp, length(b_vec) - 1)  

for(j in 1:n_samp){
  # Generate three independent (uncorrelated) predictor variables
  predictors <- cbind(rnorm(n),
                      rnorm(n),
                      rnorm(n))
  design_mat <- cbind(1, predictors)  # Add intercept term
  
  # Generate response variable based on true model
  y <- design_mat %*% b_vec + rnorm(n, sd = sig) 
  
  # Store data in a tibble
  dat <- tibble(y = as.vector(y),
                x1 = predictors[,1],
                x2 = predictors[,2],
                x3 = predictors[,3])
  
  # Fit multiple linear regression (MLR)
  full_mod <- lm(y ~ ., data = dat)
  MLR[j,] <- full_mod$coefficients[2:4]  # Store estimated slopes
  std_err_beta[j,] <- summary(full_mod)$coefficients[2:4,2]  # Store standard errors
}

# Visualize relationships between predictors (should show little to no
# correlation between the X's)
pairs(dat, upper.panel = NULL, pch = 16)

# Plot histograms of estimated coefficients (β̂)
par(mfrow = c(1,3))  # Arrange plots side by side

# Histograms show estimated distributions of β̂ when predictors are uncorrelated
hist(MLR[,1], breaks = 30, main = expression(hat(beta)[1]), col = "gray")
hist(MLR[,2], breaks = 30, main = expression(hat(beta)[2]), col = "gray")
hist(MLR[,3], breaks = 30, main = expression(hat(beta)[3]), col = "gray")


# Simulating Data with Correlated Predictors ------------------------------
# Now we introduce **multicollinearity** by generating predictors that are highly correlated.

MLR_c <- matrix(NA, n_samp, length(b_vec) - 1)
std_err_beta_c <- matrix(NA, n_samp, length(b_vec) - 1)

# Define correlations between predictors
cor_x1_x2 <- .99  # Very strong correlation between X1 and X2
cor_x1_x3 <- .99  # Very strong correlation between X1 and X3
cor_x2_x3 <- .99  # Very strong correlation between X2 and X3

for(j in 1:n_samp){
  # Generate correlated predictors using a covariance matrix
  predictors <- rmvnorm(n = n, sigma = (sig^2) * matrix(c(1, cor_x1_x2, cor_x1_x3,
                                                          cor_x1_x2, 1, cor_x2_x3,
                                                          cor_x1_x3, cor_x2_x3, 1),
                                                        byrow = TRUE,
                                                        nrow = 3,
                                                        ncol = 3))
  design_mat <- cbind(1, predictors)  # Add intercept
  
  # Generate response variable
  y <- design_mat %*% b_vec + rnorm(n, sd = sig)
  
  # Store data
  dat <- tibble(y = as.vector(y),
                x1 = predictors[,1],
                x2 = predictors[,2],
                x3 = predictors[,3])
  
  # Fit MLR with correlated predictors
  full_mod_c <- lm(y ~ ., data = dat)
  MLR_c[j,] <- full_mod_c$coefficients[2:4]  # Store estimated slopes
  std_err_beta_c[j,] <- summary(full_mod_c)$coefficients[2:4,2]  # Store standard errors
}

# Visualize relationships between predictors (should show strong correlation)
pairs(dat, upper.panel = NULL, pch = 16)

# Comparing Coefficients from Uncorrelated vs Correlated Predictors -------

# Histograms of estimated regression coefficients (β̂) before vs after multicollinearity
# Red: Uncorrelated predictors
# Blue: Correlated predictors (overlapping area will be purple)

# Beta 1 (X1 coefficient)
hist(MLR[,1], 
     xlim = range(MLR_c[,1]),
     prob = TRUE,
     main = expression(hat(beta)[1]),
     col = "red",
     border = "white")

hist(MLR_c[,1], 
     breaks = 30, 
     prob = TRUE,
     col = rgb(0,0,.5,.5),  # Transparent blue overlay
     border = "white",
     add = TRUE)

# Beta 2 (X2 coefficient)
hist(MLR[,2], 
     xlim = range(MLR_c[,2]),
     prob = TRUE,
     main = expression(hat(beta)[2]),
     col = "red",
     border = "white")

hist(MLR_c[,2], 
     breaks = 30, 
     prob = TRUE,
     col = rgb(0,0,.5,.5), 
     border = "white",
     add = TRUE)

# Beta 3 (X3 coefficient)
hist(MLR[,3], 
     xlim = range(MLR_c[,3]),
     prob = TRUE,
     main = expression(hat(beta)[3]),
     col = "red",
     border = "white")

hist(MLR_c[,3], 
     breaks = 30, 
     prob = TRUE,
     col = rgb(0,0,.5,.5), 
     border = "white",
     add = TRUE)

# -------------------------------------------------------------------------
# TAKE AWAYS:
#
# 1. Multicollinearity does **not** introduce bias.
#
# 2. With multicollinear predictors (blue histograms), the spread of β̂
# increases significantly. This means the estimates become **more unstable**.
# This increased spread translates into **higher standard errors**, meaning
# greater uncertainty in our coefficient estimates.
#
# 3. Multicollinearity makes it harder to determine which predictor truly
# influences Y, because their effects are confounded with each other.  We need
# more data to ascertain the effect of predictors in the presence of
# multicollinearity


