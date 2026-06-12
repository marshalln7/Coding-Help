# Class - Feb 18, 2025

library("plotly")

## parameters
b0 <- 5.0
b1 <- -1.0
b2 <- 2.0
sig <- 0.5

## function for surface
f <- function(x1, x2) {
  (b0 + b1*x1 + b2*x2)
}

n <- 100

X <- matrix( rnorm(n*2), nrow = n)

######### add correlation among Xs?
# rho <- 0.9
# X <- X %*% chol(matrix(c(1, rho, rho, 1), nrow = 2))
###################################

## generate response according to model
Y <- f(X[,1], X[,2]) + rnorm(n, sd = sig)


## grid of values for plotting the surface
x1 <- seq(-2.5, 2.5, length = 200)
x2 <- seq(-2.5, 2.5, length = 200)
y <- outer(x1, x2, FUN = f) |> t()


fig <- plot_ly(x = x1, y = x2, z = y, width = 900, height = 500) %>% 
  add_surface(opacity = 0.3, showscale = FALSE) %>%
  add_trace(x = X[,1], y = X[,2], z = Y, 
            type = "scatter3d", mode = "markers", 
            marker = list(color = "darkblue", size = 5.0),
            showlegend = FALSE) %>%
  layout(scene = list(aspectmode = "manual", 
                 aspectratio = list(x = 1, y = 1, z = 1),
                 camera = list(eye = list(x = -1.9, y = 1.0, z = 0.5)),
                 xaxis = list(title = "x1"),
                 yaxis = list(title = "x2"),
                 zaxis = list(title = "y")))
fig


## EDA - visualize 
plot(X[,1], Y)
plot(X[,2], Y)
plot(X[, 1], X[, 2])

## pairwise scatter plots
dat <- data.frame(y = Y, x1 = X[,1], x2 = X[,2])

library("GGally")
pairs(dat, upper.panel = NULL)
ggpairs(dat)








