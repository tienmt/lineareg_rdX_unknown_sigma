mala_fractional_posterior <- function(Y, X, alpha = 1, tau2 = 1, 
                                      a = 1, b = 1,
                                      delta = 0.01, n_iter = 10000, 
                                      theta_init = NULL, sigma2_init = 1) {
  n <- length(Y)
  d <- ncol(X)
  if (is.null(theta_init)) {
    theta <- rep(0, d)
  } else {  theta <- theta_init
  }
  sigma2 <- sigma2_init
  
  # Store samples
  theta_samples <- matrix(NA, nrow = n_iter, ncol = d)
  sigma2_samples <- numeric(n_iter)
  tX = t(X)
  ac = 0
  for (iter in 1:n_iter) {
    # ===== MALA proposal for theta =====
    resid <- Y - X %*% theta
    grad_logposter <- (alpha / sigma2) * tX %*% resid - 4 * theta / (tau2 + theta^2)
    # MALA proposal
    theta_prop <- theta + (delta^2 / 2) * grad_logposter + delta * rnorm(d)
    
    resid_prop <- Y - X %*% theta_prop
    log_like <- - (alpha / (2 * sigma2)) * sum(resid^2)
    log_like_prop <- - (alpha / (2 * sigma2)) * sum(resid_prop^2)
    
    log_post <- log_like - 2 * sum(log(tau2 + theta^2))
    log_post_prop <- log_like_prop - 2 * sum(log(tau2 + theta_prop^2))
    
    # Proposal densities q(θ' | θ) and q(θ | θ')
    grad_theta_prop <- (alpha / sigma2) * tX %*% resid_prop - 4 * theta_prop / (tau2 + theta_prop^2)
    
    log_q_prop_given_cur <- sum(dnorm(theta_prop, mean = theta + (delta^2 / 2) * grad_logposter, 
                                      sd = delta, log = TRUE))
    log_q_cur_given_prop <- sum(dnorm(theta, mean = theta_prop + (delta^2 / 2) * grad_theta_prop, 
                                      sd = delta, log = TRUE))
    # Metropolis-Hastings acceptance
    log_acc_prob <- (log_post_prop + log_q_cur_given_prop) -
      (log_post + log_q_prop_given_cur)
    if (log(runif(1)) < log_acc_prob) {
      theta <- theta_prop ; ac <- ac +1
    }
    
    # ===== Gibbs update for sigma^2 =====
    shape_post <- a + alpha * n / 2
    rate_post <- b + (alpha / 2) * sum((Y - X %*% theta)^2)
    sigma2 <- 1 / rgamma(1, shape = shape_post, rate = rate_post)
    
    # ===== Store samples =====
    theta_samples[iter, ] <- theta
    sigma2_samples[iter] <- sigma2
    accept_rate <- ac/n_iter
  }
  return(list(theta = theta_samples, sigma2 = sigma2_samples , accept_rate = accept_rate))
}