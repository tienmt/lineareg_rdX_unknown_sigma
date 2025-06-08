source('function.R') ; library(scalreg); library(glmnet)
# Simulate data
n <- 100
d <- 120
s0 = d
theta_0 <- rep(0, d) ; 
theta_0[1:s0] <- c( rep(1/10, s0/2) , rep(-1/10, s0/2) )
sigma_true <- 1
lasso_out = bayes_out_1 = bayes_out_09 = bayes_out_07 = bayes_out_05 = bayes_out_03 =  list()
sig_lasso_out = sig_bayes_out = c()

for (ss in 1:30) {
  X <- matrix(rnorm(n * d), n, d)
  Y <- X %*% theta_0 + rnorm(n)
  fit <- scalreg(X, Y) ; bet_lasso <- fit$coefficients; sig_lasso <- fit$hsigma
  sig_lasso_out[ss] = sqrt(sig_lasso)
  fit_glmnet <- cv.glmnet(X, Y, intercept = FALSE,nfolds = 5) ; bet_glmnet <- as.numeric(coef(fit_glmnet)[-1])
  samples <- mala_fractional_posterior(Y, X, alpha = 1,b = 10, tau2 = .1, delta = 0.05, n_iter = 5000,theta_init = bet_glmnet)
  bet_means <- colMeans(samples$theta)
  sig_bayes_out[ss] = sqrt(mean(samples$sigma2))
  bayes_out_1[[ss]] = c( sum((bet_means - theta_0)^2), mean((X %*%bet_means - X %*%theta_0)^2), abs(sqrt(mean(samples$sigma2)) -sigma_true ) )
  lasso_out[[ss]] = c( sum((bet_lasso - theta_0)^2), mean((X %*%bet_lasso - X %*%theta_0)^2), abs(sqrt(sig_lasso) -sigma_true ) )
  
  samples <- mala_fractional_posterior(Y, X, alpha = .9,b = 10, tau2 = .1, delta = 0.055, n_iter = 5000,theta_init = bet_glmnet)
  bet_means <- colMeans(samples$theta)
  sig_bayes_out[ss] = sqrt(mean(samples$sigma2))
  bayes_out_09[[ss]] = c( sum((bet_means - theta_0)^2), mean((X %*%bet_means - X %*%theta_0)^2), abs(sqrt(mean(samples$sigma2)) -sigma_true ) )
  
  samples <- mala_fractional_posterior(Y, X, alpha = .7,b = 10, tau2 = .1, delta = 0.08, n_iter = 5000,theta_init = bet_glmnet)
  bet_means <- colMeans(samples$theta)
  sig_bayes_out[ss] = sqrt(mean(samples$sigma2))
  bayes_out_07[[ss]] = c( sum((bet_means - theta_0)^2), mean((X %*%bet_means - X %*%theta_0)^2), abs(sqrt(mean(samples$sigma2)) -sigma_true ) )
  
  samples <- mala_fractional_posterior(Y, X, alpha = .5,b = 10, tau2 = .1, delta = 0.09, n_iter = 5000,theta_init = bet_glmnet)
  bet_means <- colMeans(samples$theta)
  sig_bayes_out[ss] = sqrt(mean(samples$sigma2))
  bayes_out_05[[ss]] = c( sum((bet_means - theta_0)^2), mean((X %*%bet_means - X %*%theta_0)^2), abs(sqrt(mean(samples$sigma2)) -sigma_true ) )
  
  samples <- mala_fractional_posterior(Y, X, alpha = .3,b = 10, tau2 = .1, delta = 0.09, n_iter = 5000,theta_init = bet_glmnet)
  bet_means <- colMeans(samples$theta)
  sig_bayes_out[ss] = sqrt(mean(samples$sigma2))
  bayes_out_03[[ss]] = c( sum((bet_means - theta_0)^2), mean((X %*%bet_means - X %*%theta_0)^2), abs(sqrt(mean(samples$sigma2)) -sigma_true ) )
  
  print(ss)
}
save.image(file = 'sim_change_alpha.rda')
