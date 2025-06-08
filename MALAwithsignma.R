source('function.R') ; library(scalreg)
# Simulate data
set.seed(123)
n <- 100
d <- 120
X <- matrix(rnorm(n * d), n, d)
s0 = 20
beta_true <- rep(0, d) ; beta_true[1:s0] <- c( rep(1, s0/2) , rep(-1, s0/2) )
sigma_true <- 1.1
Y <- X %*% theta_true + rnorm(n, sd = sigma_true)
fit <- scalreg(X, Y)
bet_lasso <- fit$coefficients; 
# Run MALA
samples <- mala_fractional_posterior(Y, X, alpha = .99, tau2 = .01, a = .1, b = 20, delta = 0.04, n_iter = 5000, theta_init = bet_lasso)

# Plot posterior means
bet_means <- colMeans(samples$theta)
c( sum((bet_means - theta_true)^2), mean((X %*%bet_means - X %*%theta_true)^2) )
c( sum((bet_lasso - theta_true)^2), mean((X %*%bet_lasso - X %*%theta_true)^2) )
sqrt(mean(samples$sigma2) ) ; samples$accept_rate
(sig_lasso <- fit$hsigma)
