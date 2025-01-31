library(truncnorm)
# functions

# Covariance function

# base covariance fn

# SE
se <- function(x, alpha_obs, rho, delta) {
  K = matrix(nrow = length(x), ncol = length(x))
  N = length(x)
  for(i in 1:N) {
    K[i, i] = alpha_obs^2 + delta
    if(i == N) {
      break
    }
    for(j in (i + 1):N) {
      K[i, j] = alpha_obs^2 * exp(-(x[i] - x[j])^2 / (2*rho^2))
      K[j, i] = K[i, j]
    }
  }
  return(K)
}

# Matern 3/2
m32 <- function(x, alpha_obs, rho, delta) {
  K = matrix(nrow = length(x), ncol = length(x))
  N = length(x)
  sq_alpha_obs = alpha_obs^2
  r = -1/rho
  
  for(i in 1:N) {
    K[i, i] = sq_alpha_obs + delta
    if(i == N) {
      break
    }
    for(j in (i + 1):N) {
      K[i, j] = (1 - (r * sqrt(3) * abs(x[i] - x[j]))) * 
        exp(r * sqrt(3) * abs(x[i] - x[j])) * sq_alpha_obs
      K[j, i] = K[i, j]
    }
  }
  return(K)
}

# Matern 5/2 
m52 <- function(x, alpha_obs, rho, delta) {
  K = matrix(nrow = length(x), ncol = length(x))
  N = length(x)
  sq_rho = rho^2;
  rho4 = rho^4
  sq_alpha_obs = alpha_obs^2
  r = -1/(2 * sq_rho)
  
  for(i in 1:N) {
    K[i, i] = sq_alpha_obs + delta
    if(i == N) {
      break
    }
    for(j in (i + 1):N) {
      K[i, j] = (1 + (sqrt(5) * abs(x[i] - x[j])/rho) + 
                   ((5 * (x[i] - x[j])^2)/ (3 * sq_rho))) * 
        exp(- sqrt(5) * abs(x[i] - x[j])/rho) * sq_alpha_obs
      K[j, i] = K[i, j]
    }
  }
  return(K)
}

# Derivative cov fns

# Deriv SE
deriv_se <- function(x, derivative, alpha_obs, alpha_grad, rho, delta) {
  N = length(x)
  K = matrix(nrow = length(x), ncol = length(x))
  sq_rho = rho^2
  rho4 = rho^4
  sq_alpha_obs = alpha_obs^2
  sq_alpha_grad = alpha_grad^2
  r = -1/(2 * sq_rho)
  
  for (i in 1:N) {
    if (derivative[i] == 0) {
      K[i, i] = sq_alpha_obs + delta
    } else if (derivative[i] == 1) {
      K[i, i] = (sq_alpha_grad / sq_rho) + delta
    }
    if(i == N) {
      break
    }
    for (j in (i + 1):N) {
      if(derivative[i] == 0 && derivative[j] == 0) {
        K[i, j] = exp(r * (x[i] - x[j])^2) * sq_alpha_obs
      } else if(derivative[i] == 0 && derivative[j] == 1) {
        K[i, j] = exp(r * (x[i] - x[j])^2) * 
          (x[i] - x[j]) * ((alpha_obs * alpha_grad) / sq_rho)
      } else if(derivative[i] == 1 && derivative[j] == 0) {
        K[i, j] = exp(r * (x[i] - x[j])^2) * 
          (x[j] - x[i]) * ((alpha_grad * alpha_obs) / sq_rho)
      } else if(derivative[i] == 1 && derivative[j] == 1) {
        K[i, j] = exp(r * (x[i] - x[j])^2) * 
          (sq_rho - (x[i] - x[j])^2) * (sq_alpha_grad / rho4)
      }
      K[j, i] = K[i, j]
    }
  }
  return(K)
}

# Deriv Matern 3/2
deriv_m32 <- function(x, derivative, alpha_obs, alpha_grad, rho, delta) {
  N = length(x)
  K = matrix(nrow = length(x), ncol = length(x))
  sq_rho = rho^2
  rho4 = rho^4
  sq_alpha_obs = alpha_obs^2
  sq_alpha_grad = alpha_grad^2
  r = -1/rho
  
  for (i in 1:N) {
    if (derivative[i] == 0) {
      K[i, i] = sq_alpha_obs + delta
    } else if (derivative[i] == 1) {
      K[i, i] = (3 * sq_alpha_grad / sq_rho) + delta
    } 
    if(i == N) {
      break
    }
    for (j in (i + 1):N) {
      if(derivative[i] == 0 && derivative[j] == 0) {
        K[i, j] = (1 - (r * sqrt(3) * abs(x[i] - x[j]))) * 
          exp(r * sqrt(3) * abs(x[i] - x[j])) * sq_alpha_obs
      } else if(derivative[i] == 0 && derivative[j] == 1) {
        K[i, j] = (3 * (x[i] - x[j]) / sq_rho) * 
          exp(r * sqrt(3) * abs(x[i] - x[j])) * alpha_obs * alpha_grad
      } else if(derivative[i] == 1 && derivative[j] == 0) {
        K[i, j] = (3 * (x[j] - x[i]) / sq_rho) * 
          exp(r * sqrt(3) * abs(x[i] - x[j])) * alpha_obs * alpha_grad
      } else if(derivative[i] == 1 && derivative[j] == 1) {
        K[i, j] = (1 + (r * sqrt(3) * abs(x[i] - x[j]))) * 
          exp(r * sqrt(3) * abs(x[i] - x[j])) * sq_alpha_grad * (3 / sq_rho)
      }
      K[j, i] = K[i, j]
    }
  }
  return(K)
}

# Deriv Matern 5/2
deriv_m52 <- function(x, d, alpha_obs, alpha_grad, rho, delta) {
  K = matrix(nrow = length(x), ncol = length(x))
  N = length(x)
  sq_rho = rho^2
  sq_alpha_obs = alpha_obs^2
  sq_alpha_grad = alpha_grad^2
  for(i in 1:N) {
    if(d[i] == 0) {
      K[i, i] = (alpha_obs^2) + delta
    } else if(d[i] == 1) {
      K[i, i] = (5 * sq_alpha_grad / (3 * sq_rho)) + delta
    }
    if(i == N) {
      break
    }
    for(j in (i + 1):N) {
      if(d[i] == 0 && d[j] == 0) {
        K[i, j] = (1 + (sqrt(5) * abs(x[i] - x[j])/rho) + 
                     ((5 * (x[i] - x[j])^2)/ (3 * sq_rho))) * 
          exp(- sqrt(5) * abs(x[i] - x[j])/rho) * sq_alpha_obs
      } else if(d[i] == 0 && d[j] == 1) {
        K[i, j] = ((5 * (x[i] - x[j]))/ (3 * sq_rho)) * 
          (1 + (sqrt(5) * abs(x[i] - x[j])/rho)) * 
          exp(- sqrt(5) * abs(x[i] - x[j])/rho) * alpha_obs * alpha_grad
      } else if(d[i] == 1 && d[j] == 0) {
        K[i, j] = ((5 * (x[j] - x[i]))/ (3 * sq_rho)) * 
          (1 + (sqrt(5) * abs(x[i] - x[j])/rho)) * 
          exp(- sqrt(5) * abs(x[i] - x[j])/rho) * alpha_obs * alpha_grad
      } else if(d[i] == 1 && d[j] == 1) {
        K[i, j] = (1 + (sqrt(5) * abs(x[i] - x[j])/rho) - 
                     ((5 * (x[i] - x[j])^2)/ sq_rho)) * 
          exp(- sqrt(5) * abs(x[i] - x[j])/rho) * (5 * sq_alpha_grad / (3 * sq_rho))
      }
      K[j, i] = K[i, j]
    }
  }
  return(K)
}


# Draw GP means from mvnorm
gp_draw <- function(draws, x, Sigma, ...) {
  mu <- rep(0, length(x))
  mvtnorm::rmvnorm(draws, mu, Sigma)
}

# Simulate GP
gp_sim_data <- function(n_obs, dims, corr, true_x, s_x, 
                        rho, alpha, sigma, covfn, delta, ...) {  # dims is no. of multioutput dimensions (>1)
  x_obs <- rnorm(n_obs, true_x, s_x)  # obs_t ~ N(true_t, s)
  K <- list()
  for (j in 1:dims){
    if (covfn == 'm32') {
      K[[j]] <- m32(x_true, alpha_obs = alpha[j], rho = rho[j], delta)
    } else if (covfn == 'm52') {
      K[[j]] <- m52(x_true, alpha_obs = alpha[j], rho = rho[j], delta)
    } else if (covfn == 'se') {
      K[[j]] <- se(x_true, alpha_obs = alpha[j], rho = rho[j], delta)
    }
  }
  # Draw from joint GP (using true_t)
  f0 <- matrix(nrow = dims, ncol = length(x_true))
  for( j in 1:dims){
    f0[j,] <- gp_draw(1, x_true, K[[j]])  # Computing GP mean
  } 
  
  C <- matrix(rep(corr,dims^2), ncol = dims) # Setting within dims output correlation
  diag(C) <- 1
  f <- t(f0) %*% chol(C)     # GP mean with output dims correlations
  y <- matrix(ncol= ncol(f), nrow = nrow(f))
  for (j in 1:ncol(f)) {
    y[,j] <- rnorm((length(f[,j])), mean = f[,j], sd = sigma[j])
  }
  data.frame(y, f, x_true, x_obs)
}

# Simulate deriv GP
deriv_gp_sim_data <- function(n_obs, dims, corr, true_x, s_x, 
                        rho, alpha_obs, alpha_grad, sigma_obs, sigma_grad, covfn, delta...) {  # dims is no. of multioutput dimensions (>1)
  
  x_true <- rep(true_x, 2)   # Rep for obs and grad
  x_obs <- rep(rnorm(length(true_x), true_x, s_x), 2)  # obs_t ~ N(true_t, s)
  deriv <- c(rep(0, n_obs), rep(1, n_obs)) # Indicator for obs and grad
  K <- list()
  for (j in 1:dims){
    if (covfn == 'm32') {
      K[[j]] <- deriv_m32(x = x_true, deriv, alpha_obs = alpha_obs[j], 
                       alpha_grad = alpha_grad[j], rho = rho[j], delta)
    } else if (covfn == 'm52') {
      K[[j]] <- deriv_m52(x_true, deriv, alpha_obs = alpha_obs[j], 
                       alpha_grad = alpha_grad[j], rho = rho[j], delta)
    } else if (covfn == 'se') {
      K[[j]] <- deriv_se(x_true, deriv, alpha_obs = alpha_obs[j], 
                       alpha_grad = alpha_grad[j], rho = rho[j], delta)
    }
  }
  # Draw from joint GP (using true_t)
  f0 <- matrix(nrow = dims, ncol = length(x_true))
  for( j in 1:dims){
    f0[j,] <- gp_draw(1, x_true, K[[j]])  # Computing GP mean
  } 
  C <- matrix(rep(corr,dims^2), ncol = dims) # Setting within dims output correlation
  diag(C) <- 1
  f <- t(f0) %*% chol(C)     # GP mean with output dims correlations
  y <- matrix(ncol= ncol(f), nrow = nrow(f))
  for (j in 1:ncol(f)) {
    y[1:n_obs,j] <- rnorm((length(f[,j])/2), mean = f[1:n_obs,j], sd = sigma_obs[j])        # Original output
    y[(n_obs+1):(2*n_obs),j] <- rnorm((length(f[,j])/2), mean = f[(n_obs+1):(2*n_obs),j], sd = sigma_grad[j])     # Derivative output
  }
  data.frame(y, f, deriv, x_true, x_obs)
}

# Periodic simulation

# For simulating multi-output periodic data
periodic_data <- function(n_obs, dims, true_x, s_x, corr, rho,
                          alpha_obs, alpha_grad, sigma_obs, sigma_grad) {
  x_true <- rep(true_x, 2)   # Rep for obs and grad
  x_obs <- rep(rnorm(length(true_x), true_x, s_x), 2)  # obs_t ~ N(true_t, s)
  deriv <- c(rep(0, n_obs), rep(1, n_obs)) # Indicator for obs and grad
  f_obs <- matrix(nrow = n_obs, ncol = dims)
  for(i in 1:n_obs){
    for(j in 1:dims){
      f_obs[i, j] <- alpha_obs[j] * sin(x_true[i]/rho[j])
    }
  }
  f_grad <- matrix(nrow = n_obs, ncol = dims)
  for(i in 1:n_obs){
    for(j in 1:dims){
      f_grad[i, j] <- (alpha_grad[j] / rho[j]) * cos(x_true[i] / rho[j])
    }
  }
  f0 <- rbind(f_obs, f_grad)
  C <- matrix(rep(corr,dims^2), ncol = dims) # Setting within dims output correlation
  diag(C) <- 1
  f <- f0 %*% chol(C)     # GP mean with output dims correlations
  y <- matrix(nrow = 2*n_obs, ncol = dims)
  for (j in 1:ncol(f)) {
    y[1:n_obs,j] <- rnorm((length(f[,j])/2), mean = f[1:n_obs,j], sd = sigma_obs[j])        # Original output
    y[(n_obs+1):(2*n_obs),j] <- rnorm((length(f[,j])/2), mean = f[(n_obs+1):(2*n_obs),j], sd = sigma_grad[j])     # Derivative output
  }
  data.frame(y, f, deriv, x_true, x_obs)
}

# For simulating multi-output periodic data
periodic_trend_data <- function(n_obs, dims, true_x, s_x, corr, rho,
                          alpha_obs, alpha_grad, sigma_obs, sigma_grad) {
  x_true <- rep(true_x, 2)   # Rep for obs and grad
  x_obs <- rep(rnorm(length(true_x), true_x, s_x), 2)  # obs_t ~ N(true_t, s)
  deriv <- c(rep(0, n_obs), rep(1, n_obs)) # Indicator for obs and grad
  f_obs <- matrix(nrow = n_obs, ncol = dims)
  for(i in 1:n_obs){
    for(j in 1:dims){
      f_obs[i, j] <- (alpha_obs[j] * sin(x_true[i]/rho[j])) + (0.5 * (x_true[i])^2) #exp(p * x_true[i])
    }
  }
  f_grad <- matrix(nrow = n_obs, ncol = dims)
  for(i in 1:n_obs){
    for(j in 1:dims){
      f_grad[i, j] <- (alpha_grad[j] * (1 / rho[j]) * cos(x_true[i] / rho[j])) + (0.25 * x_true[i]) #(p * exp(x_true[i]))
    }
  }
  f0 <- rbind(f_obs, f_grad)
  C <- matrix(rep(corr,dims^2), ncol = dims) # Setting within dims output correlation
  diag(C) <- 1
  f <- f0 %*% chol(C)     # GP mean with output dims correlations
  y <- matrix(nrow = 2*n_obs, ncol = dims)
  for (j in 1:ncol(f)) {
    y[1:n_obs,j] <- rnorm((length(f[,j])/2), mean = f[1:n_obs,j], sd = sigma_obs[j])        # Original output
    y[(n_obs+1):(2*n_obs),j] <- rnorm((length(f[,j])/2), mean = f[(n_obs+1):(2*n_obs),j], sd = sigma_grad[j])     # Derivative output
  }
  data.frame(y, f, deriv, x_true, x_obs)
}

# Exact GP model
gp_model <- function(stanmodel, n_obs, dims, outputs, inputs, latent_sd, latent_inputs, covfn, 
                     rho_prior, ls_param, msd_param, esd_param, msd_param_grad, esd_param_grad, 
                     mean_obs, 
                     sd_obs, 
                     mean_grad, 
                     sd_grad, 
                     nugget, is_deriv, is_scale,is_vary, 
                     is_corr, adapt_delta, iter, warmup, chains, cores, init) {
  if (is_deriv == 1) {
    derivative <- c(rep(0, n_obs), rep(1, n_obs))
    N <- 2* n_obs
  } else {
    derivative <- rep(0, n_obs)
    N <- n_obs
    outputs <- outputs[1:n_obs,]
  }
  gp_data <- list(
    N = N,
    M = n_obs,
    D = dims,
    y = outputs,
    derivative = derivative,
    inputs = inputs,
    s = latent_sd,
    latent = latent_inputs,
    covfn = covfn,
    rho_prior = rho_prior,
    ls_param = ls_param,
    msd_param = msd_param,
    esd_param = esd_param,
    msd_param_grad = msd_param_grad,
    esd_param_grad = esd_param_grad,
    mean_obs = mean_obs,
    sd_obs = sd_obs,
    mean_grad = mean_grad,
    sd_grad = sd_grad,
    nugget = nugget,
    is_deriv = is_deriv,
    is_scale = is_scale,
    is_vary = is_vary,
    is_corr = is_corr
  )
  #model fitting
  fit <- sampling(
    object = stanmodel,
    data = gp_data,
    init = init,
    chains = chains,
    warmup = warmup,
    iter = iter,
    cores = cores,
    control = list(adapt_delta = adapt_delta)
  )
  # fit <- stanmodel$sample(
  #   data = gp_data,
  #   init = init,
  #   chains = chains,
  #   iter_warmup = warmup,
  #   iter_sampling = (iter - warmup),
  #   parallel_chains = cores,
  #   adapt_delta = adapt_delta
  # )
  return(fit)
}

# Specify parameters

# Parameters for GP constant across outputs
true_params <- function(n_obs, dims, msd_pars, esd_pars, ls_pars){
  alpha <- rep(rtruncnorm(1, a = 0.1, b = Inf, mean = msd_pars[1], sd = msd_pars[2]), dims)        # GP marginal SD
  sigma <- rep(rtruncnorm(1, a = 0.1, b = Inf, mean = esd_pars[1], sd = esd_pars[2]), dims)        # Error SD 
  rho <- rep(rtruncnorm(1, a = 0.1, b = Inf, mean = ls_pars[1], sd = ls_pars[2]), dims)        # GP length scale
  data.frame(alpha, sigma, rho)
}

# Parameters for GP varying across outputs
true_params_vary <- function(n_obs, dims, msd_pars, esd_pars, ls_pars){
  alpha <- rtruncnorm(dims, a = 0.1, b = Inf, mean = msd_pars[1], sd = msd_pars[2])      # GP marginal SD
  sigma <- rtruncnorm(dims, a = 0.1, b = Inf, mean = esd_pars[1], sd = esd_pars[2])       # Error SD 
  rho <- rtruncnorm(dims, a = 0.1, b = Inf, mean = ls_pars[1], sd = ls_pars[2])           # GP length scale
  data.frame(alpha, sigma, rho)
}

# Parameters for deriv GP constant across outputs
true_params_deriv <- function(n_obs, dims, deriv_scale, msd_pars, esd_pars, ls_pars){
  alpha_grad <- rep(rtruncnorm(1, a = 0.1, b = Inf, mean = msd_pars[1], sd = msd_pars[2]), dims)   # GP marginal SD for f'
  alpha_obs <- deriv_scale * alpha_grad       # GP marginal SD for f
  sigma_grad <- rep(rtruncnorm(1, a = 0.1, b = Inf, mean = esd_pars[1], sd = esd_pars[2]), dims)   # Error SD for y'
  sigma_obs <- deriv_scale * sigma_grad       # Error SD for y
  rho <- rep(rtruncnorm(1, a = 0.1, b = Inf, mean = ls_pars[1], sd = ls_pars[2]), dims)        # GP length scale
  data.frame(alpha_obs, alpha_grad, sigma_obs, sigma_grad, rho)
}

# Parameters for deriv GP varying across outputs
true_params_deriv_vary <- function(n_obs, dims, deriv_scale, msd_pars, esd_pars, ls_pars){
  alpha_grad <- rtruncnorm(dims, a = 0.1, b = Inf, mean = msd_pars[1], sd = msd_pars[2])           # GP marginal SD for f'
  alpha_obs <- deriv_scale * alpha_grad       # GP marginal SD for f
  sigma_grad <- rtruncnorm(dims, a = 0.1, b = Inf, mean = esd_pars[1], sd = esd_pars[2])          # Error SD for y'
  sigma_obs <- deriv_scale * sigma_grad       # Error SD for y
  rho <- rtruncnorm(dims, a = 0.1, b = Inf, mean = ls_pars[1], sd = ls_pars[2])                # GP length scale
  data.frame(alpha_obs, alpha_grad, sigma_obs, sigma_grad, rho)
}
true_params_deriv_vary_org <- function(n_obs, dims, deriv_scale, msd_pars, esd_pars, ls_pars){
  alpha_grad <- runif(dims, msd_pars[1], msd_pars[2])          # GP marginal SD for f'
  alpha_obs <- deriv_scale * alpha_grad       # GP marginal SD for f
  sigma_grad <- runif(dims, esd_pars[1], esd_pars[2])          # Error SD for y'
  sigma_obs <- deriv_scale * sigma_grad       # Error SD for y
  rho <- runif(dims, ls_pars[1], ls_pars[2])             # GP length scale
  data.frame(alpha_obs, alpha_grad, sigma_obs, sigma_grad, rho)
}

# Other functions for summary
abs_bias_draws <- function(theta_hat, theta) {
  abs(mean(theta_hat) - theta)
}
rmse_draws <- function(theta_hat, theta) {
  sqrt(mean((theta_hat - theta)^2))
}

mae_draws <- function(theta_hat, theta) {
  mean(abs(theta_hat - theta))
}

# Sim summary functions
## Summary table function
compare_summary <- function(model, variable, dims, true_variable, n_obs, 
                            variable_class, sim_id, deriv, scale, vary, corr) {
  out <- as_draws_matrix(model)
  subset <- subset_draws(out, variable = variable)
  summary <- summarise_draws(subset)
  mean <- summary$mean
  rhat <- summary$rhat
  bess <- summary$ess_bulk
  tess <- summary$ess_tail
  sd <- rep(NA, length(variable))
  abs_bias <- rep(NA, length(variable))
  rmse <- rep(NA, length(variable))
  mae <- rep(NA, length(variable))
  for(i in 1:length(variable)) {
    sd[i] <- sd(subset[,i])
    abs_bias[i] <- abs_bias_draws(subset[,i], true_variable[i])
    rmse[i] <- rmse_draws(subset[,i], true_variable[i])
    mae[i] <- mae_draws(subset[,i], true_variable[i])
  }
  true_value <- true_variable
  class <- rep(variable_class, length(variable))
  n <- rep(n_obs, length(variable))
  d <- rep(dims, length(variable))
  pars <- variable
  runtime <- rep(max(rowSums(get_elapsed_time(model))), length(variable))
  sampler_params <- get_sampler_params(model)
  divergent_by_chain <- sapply(sampler_params, function(x) sum(x[, "divergent__"]))
  divergent <- rep(max(divergent_by_chain), length(variable))
  sim_id <- rep(sim_id, length(variable))
  ranks <- colSums(sweep(subset, 2, true_variable, "<"))
  deriv <- rep(deriv, length(variable))
  scale <- rep(scale, length(variable))
  vary <- rep(vary, length(variable))
  corr <- rep(corr, length(variable))
  data.frame(sim_id, 
             n, 
             d, 
             class, 
             pars, 
             true_value, 
             mean, 
             sd, 
             abs_bias, 
             rmse,
             mae, 
             ranks, 
             rhat, 
             bess, 
             tess, 
             deriv, 
             scale, 
             vary, 
             corr, 
             runtime, 
             divergent
             )
}

# plotting functions
gp_plot <- function (dim_id, output, gpfns, input, msd, esd, plot_type, label) {
  plot_data <- data.frame(output[, dim_id], gpfns[, dim_id], input, 
                          rep(msd[dim_id], N), rep(esd[dim_id], N))
  colnames(plot_data) <- c('output', 'gpfns', 'input', 'msd', 'esd')
  if (plot_type == 'gpfns') {
    plot <- ggplot(data = plot_data, aes(x = input, y = output)) +
      theme_bw(base_size = 20, base_family = 'Times') +
      geom_point() +
      geom_line(aes(y = gpfns), colour = c("#0072B2")) +
      geom_line(aes(y = gpfns - msd), colour = 'red', linetype = 'dashed') +
      geom_line(aes(y = gpfns + msd), colour = 'red', linetype = 'dashed') +
      labs(x = 'input', y = 'gpfns', title = label) 
  } else {
    plot <- ggplot(data = plot_data, aes(x = input, y = output)) +
      theme_bw(base_size = 20, base_family = 'Times') +
      geom_point() +
      geom_line(aes(y = output), colour = c("#0072B2")) +
      geom_line(aes(y = output - esd), colour = 'red', linetype = 'dashed') +
      geom_line(aes(y = output + esd), colour = 'red', linetype = 'dashed') +
      labs(x = 'input', y = 'outputs', title = label)
  }
  return(plot)
}
summary_plot <- function(data, x, y, color, facet_cond, ref_point, ref_point_color, scale_color1, scale_color2,
                         title, xlab, ylab, color_lab) {
  ggplot(data, aes(x = x, y = y, color = color)) +
    theme_bw(base_size=65,
             base_family = 'Times') +
    geom_point(size = 5.5 ,
               position = position_dodge(width = 0.7)) +
    geom_errorbar(aes(ymin = lower__, ymax = upper__),
                  width = 0.5,
                  linewidth = 1.5,
                  position = position_dodge(width = 0.7)) +
    geom_point(aes(x=0.5, y=ref_point), color = ref_point_color, size = 5.5) +
    annotate('text', x = 0.7, y = 0.39, label = 'Prior', size = 10, colour = ref_point_color) +
    facet_wrap(~facet_cond) +
    labs(x = xlab, y = ylab, color = color_lab) +
    guides(fill = 'none') + 
    theme(axis.ticks = element_line(linewidth = 3)) +
    scale_colour_manual(values = c(scale_color1, scale_color2)) + ggtitle(title)
}