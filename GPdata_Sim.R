# Here we simulate the GP data scenario for DGP-LVM and other GP model specifications
# We run multiple trials to verify the effects of the modifications in covariance functions
# Ground truth data is simulated from derivative GP with all our proposed modifications

#libraries
library(rstan)
library(parallel)
library(doParallel)
library(foreach)
library(posterior)
# tempdir set (set temp dir so that it doesn't overload disk: primarily for cluster computing)
Sys.setenv(TMPDIR = "/mnt/volume")
unlink(tempdir(), recursive = TRUE)
tempdir(check = TRUE)
#functions
# Specify ground truth covariance function hyperparameters
true_params <- function(n_obs, dims, deriv_scale, ...){
  alpha_grad <- runif(dims,0.4,0.6)          # GP marginal SD for derivative part
  alpha_obs <- alpha_grad * deriv_scale    # GP marginal SD for original part
  sigma_grad <- runif(dims,0.05,0.15)          # Error SD for derivative output
  sigma_obs <- sigma_grad * deriv_scale    # Error SD for original output
  rho <- runif(dims, 0.5, 1)               # GP length scale
  data.frame(alpha_obs, alpha_grad, sigma_obs, sigma_grad, rho)
}
# Generate ground truth for latent x (equidistant true x)
true_x <- function(n_obs, start, end, by){
  x <- seq(start, end, by = by)
  if (length(x) == n_obs) return(x)
  else print('length does not match no. of obs')
}
# derivative covariance function ( to generate simulated data from a derivative GP)
kernel <- function(x, d, alpha_obs, alpha_grad, rho, delta) {
  K = matrix(nrow = length(x), ncol = length(x))
  N = length(x)
  for(i in 1:N) {
    if(d[i] == 0) {
      K[i, i] = (alpha_obs^2) + delta
    } else if(d[i] == 1) {
      K[i, i] = (alpha_grad^2/rho^2) + delta
    }
    if(i == N) {
      break
    }
    for(j in (i + 1):N) {
      if(d[i] == 0 && d[j] == 0) {
        K[i, j] = alpha_obs^2 * exp(-(x[i] - x[j])^2 / (2 * rho^2))
      } else if(d[i] == 0 && d[j] == 1) {
        K[i, j] = (alpha_obs * alpha_grad) / rho^2 * (x[i] - x[j]) * exp(-(x[i] - x[j])^2 / (2*rho^2))
      } else if(d[i] == 1 && d[j] == 0) {
        K[i, j] = (alpha_obs * alpha_grad) / rho^2 * (x[j] - x[j]) * exp(-(x[i] - x[j])^2 / (2*rho^2))
      } else if(d[i] == 1 && d[j] == 1) {
        K[i, j] = alpha_grad^2 / rho^4 * (rho^2 - (x[i] - x[j])^2) * exp(-(x[i] - x[j])^2 / (2*rho^2))
      }
      K[j, i] = K[i, j]
    }
  }
  return(K)
}

# GP f draw
gp_draw <- function(draws, x, Sigma, ...) {
  mu <- rep(0, length(x))
  mvtnorm::rmvnorm(draws, mu, Sigma)
}
# Multi-output gp sim (Simulate the GP data)
gp_sim_data <- function(n_obs, dims, corr, true_x, s_x, 
                        rho, alpha_obs, alpha_grad, sigma_obs, sigma_grad, ...){  # dims is no. of multioutput dimensions (>1)
  x_true <- rep(true_x, 2)   # Rep for obs and grad
  x_obs <- rep(rnorm(length(true_x), true_x, s_x), 2)  # obs_t ~ N(true_t, s)
  deriv <- c(rep(0, n_obs), rep(1, n_obs)) # Indicator for obs and grad
  delta <- 1e-12
  K <- list()
  for (j in 1:dims){
    K[[j]] <- kernel(x_true, deriv, alpha_obs = alpha_obs[j], alpha_grad = alpha_grad[j], rho = rho[j], delta)
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
    y[1:20,j] <- rnorm((length(f[,j])/2), mean = f[1:20,j], sd = sigma_obs[j])        # Original output
    y[21:40,j] <- rnorm((length(f[,j])/2), mean = f[21:40,j], sd = sigma_grad[j])     # Derivative output
  }
  data.frame(y, x_true, x_obs, deriv)
}
# Functions for model summary (to verify recovery of true latent x)
abs_bias_draws <- function(theta_hat, theta) {
  abs(mean(theta_hat) - theta)
}
rmse_draws <- function(theta_hat, theta) {
  sqrt(mean((theta_hat - theta)^2))
}
## Specify model conditions
n_obs <- 20               # Sample size
dims <- c(5, 10, 20)      # Specify output dimensions
s_x <- 0.3                # prior measurement SD for latent x
deriv_scale <- 3          # Scale difference between original and derivative outputs
corr <- 0.5               # between output dimensions true correlation
## model conditions sensible full cross
multi_model_condition_matrix <- matrix(c(1, 1, 1, 1,
                                         1, 1, 1, 0,
                                         1, 1, 0, 1,
                                         1, 0, 1, 1,
                                         1, 1, 0, 0,
                                         1, 0, 1, 0,
                                         1, 0, 0, 1,
                                         1, 0, 0, 0,
                                         0, 0, 1, 1,
                                         0, 0, 1, 0,
                                         0, 0, 0, 1,
                                         0, 0, 0, 0), ncol = 4, byrow = TRUE)
multi_models <- nrow(multi_model_condition_matrix) # Number of sensible multi-output model condition combinations per dims set

# Sort variable names for output
x_names <- sprintf('x[%s]', seq(1:n_obs))
variable_names <- list()
imp_var <- list()
dims_names <- c()
for( k in 1:length(dims)){
  rho_names <- sprintf('rho[%s]', seq(1:dims[k]))
  alpha_obs_names <- sprintf('alpha_obs[%s]', seq(1:dims[k]))
  alpha_grad_names <- sprintf('alpha_grad[%s]', seq(1:dims[k]))
  sigma_obs_names <- sprintf('sigma_obs[%s]', seq(1:dims[k]))
  sigma_grad_names <- sprintf('sigma_grad[%s]', seq(1:dims[k]))
  variable_names[[k]] <- c(rho_names, alpha_obs_names, alpha_grad_names, 
                           sigma_obs_names, sigma_grad_names)
  imp_var[[k]] <- c(variable_names[[k]], x_names)
  dims_names[k] <- c(sprintf('[%s]D', dims[k]))
}
var_names <- list()
model_dim <- list()
is_deriv <- list()
is_scale <- list()
is_vary <- list()
is_corr <- list()
for(k in 1:length(dims)){
  var_names[[k]] <- rep(imp_var[[k]], multi_models)
  model_dim[[k]] <- rep(dims_names[k], each = length(imp_var[[k]])*multi_models)
  # Extended model condition frame
  is_deriv[[k]] <- c(rep(multi_model_condition_matrix[,1], each = length(imp_var[[k]])))
  is_scale[[k]] <- c(rep(multi_model_condition_matrix[,2], each = length(imp_var[[k]])))
  is_vary[[k]] <- c(rep(multi_model_condition_matrix[,3], each = length(imp_var[[k]])))
  is_corr[[k]] <- c(rep(multi_model_condition_matrix[,4], each = length(imp_var[[k]])))
}
# Unlist for final data frame output
var_names <- unlist(var_names)
model_dim <- unlist(model_dim)
is_deriv_list <- unlist(is_deriv)
is_scale_list <- unlist(is_scale)
is_vary_list <- unlist(is_vary)
is_corr_list <- unlist(is_corr)
#start cluster
trials <- 50
temp_summary_list <- list()
# Filename for saved trial results
trial_summary_names <- sprintf('GP_summary_trial[%s]', seq(1:trials))
cores = 25
cl <- parallel::makeCluster(cores, type="FORK")
doParallel::registerDoParallel(cl)
model_comp_list = foreach(i = 1:trials) %dopar% {
  # Check if file exists on disk (only needed when re-running trials)
  if(file.exists(trial_summary_names[i])){
  temp_summary_list[[i]] <- readRDS(trial_summary_names[i],'.rds')
  return(temp_summary_list)
  } else {
  set.seed(i)
  # Declare all the summary variables for final output
  true <- list()
  mdgp <- list()
  gp_out <- list()
  sd <- list()
  bias <- list()
  rmse <- list()
  sd_list <- list()
  bias_list <- list()
  rmse_list <- list()
  # For model diagnostics summary
  draws <- list()
  variable_draws <- list()
  summary <- list()
  runtime <- list()
  mean <- list()
  rhat <- list()
  ess_bulk <- list()
  ess_tail <- list()
  mean_list <- list()
  rhat_list <- list()
  ess_bulk_list <- list()
  ess_tail_list <- list()
  
  # Extended model condition frame
  is_deriv <- c()
  is_scale <- c()
  is_vary <- c()
  is_corr <- c()
  data <- list()
  mdgp_data <- list()
  mdgp_fit <- list()
  sparams <- list()
  sparams_prior_sd <- c()
  params <- list()
  for( k in 1:length(dims)) {
    # Generate simulated data
    params[[k]] <- true_params(n_obs, dims[k], deriv_scale)
    x_true <- true_x(n_obs, 0.5, 10, 0.5)
    data[[k]] <- gp_sim_data(n_obs, dims[[k]], corr, x_true, s_x, 
                             params[[k]]$rho, params[[k]]$alpha_obs, 
                             params[[k]]$alpha_grad, params[[k]]$sigma_obs,
                             params[[k]]$sigma_grad)
    # Sample size; y = [,1:dims[k]]
    N <- nrow(data[[k]][,1:dims[k]])
    # Number of output dimensions
    D <- ncol(data[[k]][,1:dims[k]])
    # Compute prior SD for GP marginal and error SDs
    sparams[[k]] <- apply(data[[k]][,1:dims[k]], 2, sd)
    sparams_prior_sd[k] <- mean(sparams[[k]])
    mdgp_data[[k]] <- list()
    mdgp_fit[[k]] <- list()
    for (j in 1:multi_models){
      mdgp_data[[k]] <- list(N = N,
                             D = D,
                             M = N / 2,
                             y = data[[k]][,1:dims[k]],
                             t = data[[k]]$x_obs,
                             s = s_x,
                             sparams_prior_sd = sparams_prior_sd[k],
                             derivative = data[[k]]$deriv,
                             is_deriv = multi_model_condition_matrix[j,1],
                             is_scale = multi_model_condition_matrix[j,2],
                             is_vary = multi_model_condition_matrix[j,3],
                             is_corr = multi_model_condition_matrix[j,4])
      mdgp_fit[[k]][[j]] <- stan(
        file = 'DGPLVM_se.stan',
        data   = mdgp_data[[k]],
        iter   = 3000,
        warmup = 1000,
        chains = 1,
        cores = 1,
        init = 0,
        control = list(adapt_delta=0.9))
    }
  }
  for (k in 1:length(dims)) {
    # Set true data vector for comparison
    true[[k]] <- c(params[[k]]$rho[1:dims[k]], params[[k]]$alpha_obs[1:dims[k]], params[[k]]$alpha_grad[1:dims[k]], 
                   params[[k]]$sigma_obs[1:dims[k]], params[[k]]$sigma_grad[1:dims[k]], data[[k]]$x_true[1:n_obs])
    # store all summary results
    gp_out[[k]] <- list()
    mdgp[[k]] <- list()
    # latent x summary
    sd[[k]] <- matrix(nrow = multi_models, ncol = length(imp_var[[k]]))
    bias[[k]] <- matrix(nrow = multi_models, ncol = length(imp_var[[k]]))
    rmse[[k]] <- matrix(nrow = multi_models, ncol = length(imp_var[[k]]))
    for(j in 1:multi_models){
      gp_out[[k]][[j]] <- as_draws_matrix(mdgp_fit[[k]][[j]])
      mdgp[[k]][[j]] <- subset_draws(gp_out[[k]][[j]], variable = imp_var[[k]])
      for( l in 1:length(imp_var[[k]])){
        sd[[k]][j, l] <- sd(mdgp[[k]][[j]][,l])
        bias[[k]][j, l] <- abs_bias_draws(mdgp[[k]][[j]][,l], true[[k]][l])
        rmse[[k]][j, l] <- rmse_draws(mdgp[[k]][[j]][,l], true[[k]][l])
      }
    }
    # Model summary
    sd_list[[k]] <- c(t(sd[[k]]))
    bias_list[[k]] <- c(t(bias[[k]]))
    rmse_list[[k]] <- c(t(rmse[[k]]))
    # Model diags
    draws[[k]] <- list()
    variable_draws[[k]] <- list()
    summary[[k]] <- list()
    runtime[[k]] <- list()
    mean[[k]] <- matrix(nrow = length(imp_var[[k]]), ncol = multi_models)
    rhat[[k]] <- matrix(nrow = length(imp_var[[k]]),ncol = multi_models) 
    ess_bulk[[k]] <- matrix(nrow = length(imp_var[[k]]),ncol = multi_models) 
    ess_tail[[k]] <- matrix(nrow = length(imp_var[[k]]),ncol = multi_models)
    for(j in 1:multi_models){
      draws[[k]][[j]] <- as_draws_df(mdgp_fit[[k]][[j]])
      variable_draws[[k]][[j]] <- subset(draws[[k]][[j]], variable = imp_var[[k]])
      summary[[k]][[j]] <- summarise_draws(variable_draws[[k]][[j]])
      runtime[[k]][[j]] <- rowSums(get_elapsed_time(mdgp_fit[[k]][[j]]))
      mean[[k]][,j] <- summary[[k]][[j]]$mean
      rhat[[k]][,j] <- summary[[k]][[j]]$rhat 
      ess_bulk[[k]][,j] <- summary[[k]][[j]]$ess_bulk
      ess_tail[[k]][,j] <- summary[[k]][[j]]$ess_tail 
    }
    runtime[[k]] <- unlist(runtime[[k]])
    # Posterior estimates and model convergence measures
    mean_list[[k]] <- c(mean[[k]])
    rhat_list[[k]] <- c(rhat[[k]])
    ess_bulk_list[[k]] <- c(ess_bulk[[k]])
    ess_tail_list[[k]] <- c(ess_bulk[[k]])
  }
  # Compile for final output
  sd_summary <- c()
  bias_summary <- c()
  rmse_summary <- c()
  mean_summary <- c()
  rhat_summary <- c()
  ess_bulk_summary <- c()
  ess_tail_summary <- c()
  runtime_summary <- list()
  dims_list <- list()
  true_value_list <- list()
  for( k in 1:length(dims)){
    runtime_summary[[k]] <- rep(runtime[[k]], each = length(imp_var[[k]]))
    true_value_list[[k]] <- rep(true[[k]], multi_models)
  }
  sd_summary <- unlist(sd_list)
  bias_summary <- unlist(bias_list)
  rmse_summary <- unlist(rmse_list)
  mean_summary <- unlist(mean_list)
  rhat_summary <- unlist(rhat_list)
  ess_bulk_summary <- unlist(ess_bulk_list)
  ess_tail_summary <- unlist(ess_tail_list)
  runtime_summary <- unlist(runtime_summary)
  true_value_list <- unlist(true_value_list)
  # Final output in a data frame
  temp_summary_list[[i]] <- data.frame(var_names, true_value = true_value_list, mean = mean_summary, sd = sd_summary,
                                       abs_bias = bias_summary, rmse = rmse_summary, rhat = rhat_summary, 
                                       ess_bulk = ess_bulk_summary, ess_tail = ess_tail_summary, runtime = runtime_summary, 
                                       model_dim, is_deriv = is_deriv_list, is_scale = is_scale_list, 
                                       is_vary = is_vary_list, is_corr = is_corr_list)
  saveRDS(temp_summary_list[[i]], file = paste(trial_summary_names[i],'.RDS'))
  return(temp_summary_list)
  }
}
saveRDS(model_comp_list,'Model_Summary')
