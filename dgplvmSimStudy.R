# Here we simulate the GP data scenario for DGP-LVM and other GP model specifications
# We run multiple trials to verify the effects of the modifications in covariance functions
# Ground truth data is simulated from derivative GP with all our proposed modifications

#libraries
library(rstan)
library(parallel)
library(doParallel)
library(foreach)
library(posterior)
library(data.table)
# Source fns file
source('DgplvmSimFns.R')
stanmodel <- stan_model('DerivGPmodels.stan')
# tempdir set (set temp dir so that it doesn't overload disk: primarily for cluster computing)
Sys.setenv(TMPDIR = "/mnt/volume")
unlink(tempdir(), recursive = TRUE)
tempdir(check = TRUE)
#functions

## Specify data conditions
N <- 20  # higher N will take very long. Better to stick to <=30
dims <- c(2, 5, 10)
# Set true input x
x_true <- seq(0.5, 9.5, length.out = N)
s <- 0.3 # latent sd # needed in case of latent GPs
# Generate parameters
marginal_sd_params <- c(3, 0.25)
error_sd_params <- c(1, 0.25)
ls_params <- c(1, 0.05)
corr <- 0.5
delta <- 1e-12
lambda <- 3        # scale between y and y'
covfn <- 'se'

# Generate sim data
n_sim <- 50
params <- list()
simdata <- list()
for (i in 1:n_sim) {
  set.seed(i)
  params[[i]] <- list()
  simdata[[i]] <- list() 
  for (j in 1:length(dims)) {
    params[[i]][[j]] <- true_params_deriv_vary(n_obs = N, 
                                               dims = dims[j], 
                                               deriv_scale = lambda,
                                               msd_pars = marginal_sd_params,
                                               esd_pars = error_sd_params,
                                               ls_pars = ls_params)
    simdata[[i]][[j]] <- deriv_gp_sim_data(n_obs = N, 
                                             dims = dims[j],
                                             corr = corr, 
                                             true_x = x_true, 
                                             s_x = s, rho = params[[i]][[j]]$rho,
                                             alpha_obs = params[[i]][[j]]$alpha_obs, 
                                             alpha_grad = params[[i]][[j]]$alpha_grad,
                                             sigma_obs = params[[i]][[j]]$sigma_obs,
                                             sigma_grad = params[[i]][[j]]$sigma_grad,
                                             covfn = covfn,
                                             delta = delta
    )
  }
}
saveRDS(params, 'TrueParams.rds')
saveRDS(simdata, 'Simdata.rds')
# Model specifications
rho_prior_model <- 0 # 0 = normal; 1 = invgamma;
# Specify priors for the hyperparameters (check if they differ from the simulation setting)
lambda <- 3  # derivative scale factor
marginal_sd_params_grad_model <- c(3, 0.25) # alpha_grad ~ N( , )
error_sd_params_grad_model <- c(1, 0.25)    # sigma_grad ~ N( , )
marginal_sd_params_model <- marginal_sd_params_grad_model * lambda  # alpha_obs ~ N( , )
error_sd_params_model <- error_sd_params_grad_model * lambda  # sigma_obs ~ N( , )
ls_params_model <- c(1, 0.05) # rho ~ N( , ) or InvGamma( , )
covfn_model <- 0  # 0 = SE; 1 = Matern3/2; 2 = Matern 5/2
delta_model <- 1e-6  # set higher for the derivative model if needed
adapt_delta_model <- 0.95
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
                                         0, 0, 0, 0
), ncol = 4, byrow = TRUE)
multi_models <- nrow(multi_model_condition_matrix) # Number of sensible multi-output model condition combinations per dims set

# Sort variable names for output
x_names <- sprintf('x[%s]', seq(1:N))
rho_names <- list()
alpha_obs_names <- list()
alpha_grad_names <- list()
sigma_obs_names <- list()
sigma_grad_names <- list()
for( j in 1:length(dims)){
  rho_names[[j]] <- sprintf('rho[%s]', seq(1:dims[j]))
  alpha_obs_names[[j]] <- sprintf('alpha_obs[%s]', seq(1:dims[j]))
  alpha_grad_names[[j]] <- sprintf('alpha_grad[%s]', seq(1:dims[j]))
  sigma_obs_names[[j]] <- sprintf('sigma_obs[%s]', seq(1:dims[j]))
  sigma_grad_names[[j]] <- sprintf('sigma_grad[%s]', seq(1:dims[j]))
}

#start cluster
n_sim <- 50
# parallel for model fits
out_gp <- list()
out_gp1 <- list()
out_gp2 <- list()
cores <- 50
cl <- parallel::makeCluster(cores, type="PSOCK")
doParallel::registerDoParallel(cl)
# model fit list
model_fit = foreach(i = 1:n_sim) %dopar% {
  library(rstan)
  library(posterior)
  library(data.table)
  # Check if file exists on disk (only needed when re-running trials)
  if(file.exists(paste0('temp_summary_trial','_',i,'.rds'))) {
    out_gp[[i]] <- readRDS(paste0('temp_summary_trial','_',i,'.rds'))
    return(out_gp[[i]])
  } else {
    #set.seed(50+i)
    out_gp1[[i]] <- list()
    out_gp2[[i]] <- list()
    for (j in 1:length(dims)) {
      out_gp2[[i]][[j]] <- list()
      y <- simdata[[i]][[j]][,1:dims[j]]
      x_obs <- simdata[[i]][[j]]$x_obs[1:N]
      x_true <- simdata[[i]][[j]]$x_true[1:N]
      mean_obs <- rep(NA, dims[j])
      sd_obs <- rep(NA, dims[j])
      mean_grad <- rep(NA, dims[j])
      sd_grad <- rep(NA, dims[j])
      for(d in 1:dims[j]) {
        mean_obs[d] <- mean(y[1:N,d])
        sd_obs[d] <- sd(y[1:N,d])
        mean_grad[d] <- mean(y[(N+1):(2*N),d])
        sd_grad[d] <- sd(y[(N+1):(2*N),d])
      }
      for (k in 1:nrow(multi_model_condition_matrix)) {
        gp_fit <- gp_model(stanmodel = stanmodel,
                           n_obs = N, 
                           dims = dims[j], 
                           outputs = y, 
                           inputs = x_obs, 
                           latent_sd = s, 
                           latent_inputs = 1, 
                           covfn = covfn_model, 
                           rho_prior = rho_prior_model,
                           ls_param = ls_params_model, 
                           msd_param = marginal_sd_params_model, 
                           esd_param = error_sd_params_model, 
                           msd_param_grad = marginal_sd_params_grad_model,
                           esd_param_grad = error_sd_params_grad_model,
                           mean_obs = mean_obs,
                           sd_obs = sd_obs,
                           mean_grad = mean_grad,
                           sd_grad = sd_grad,
                           nugget = delta_model, 
                           is_deriv = multi_model_condition_matrix[k,1],
                           is_scale = multi_model_condition_matrix[k,2],
                           is_vary = multi_model_condition_matrix[k,3],
                           is_corr = multi_model_condition_matrix[k,4],
                           iter = 3000, 
                           warmup = 1000, 
                           chains = 1, 
                           cores = 1, 
                           init = 0, 
                           adapt_delta = adapt_delta_model)
        out_gp_x <- compare_summary(model = gp_fit, 
                                    variable = x_names, 
                                    true_variable = x_true, 
                                    n_obs = N, 
                                    dims = dims[j], 
                                    variable_class = 'x',
                                    sim_id = i,
                                    deriv = multi_model_condition_matrix[k,1],
                                    scale = multi_model_condition_matrix[k,2],
                                    vary = multi_model_condition_matrix[k,3],
                                    corr = multi_model_condition_matrix[k,4])
        out_gp_rho <- compare_summary(model = gp_fit, 
                                      variable = rho_names[[j]], 
                                      true_variable = params[[i]][[j]]$rho, 
                                      n_obs = N, 
                                      dims = dims[j],
                                      variable_class = 'rho',  
                                      sim_id = i,
                                      deriv = multi_model_condition_matrix[k,1],
                                      scale = multi_model_condition_matrix[k,2],
                                      vary = multi_model_condition_matrix[k,3],
                                      corr = multi_model_condition_matrix[k,4])
        out_gp_alpha_obs <- compare_summary(model = gp_fit, 
                                            variable = alpha_obs_names[[j]], 
                                            true_variable = params[[i]][[j]]$alpha_obs, 
                                            n_obs = N,
                                            dims = dims[j], 
                                            variable_class = 'alpha_obs',
                                            sim_id = i,
                                            deriv = multi_model_condition_matrix[k,1],
                                            scale = multi_model_condition_matrix[k,2],
                                            vary = multi_model_condition_matrix[k,3],
                                            corr = multi_model_condition_matrix[k,4])
        out_gp_alpha_grad <- compare_summary(model = gp_fit, 
                                             variable = alpha_grad_names[[j]], 
                                             true_variable = params[[i]][[j]]$alpha_grad, 
                                             n_obs = N,
                                             dims = dims[j], 
                                             variable_class = 'alpha_grad',
                                             sim_id = i,
                                             deriv = multi_model_condition_matrix[k,1],
                                             scale = multi_model_condition_matrix[k,2],
                                             vary = multi_model_condition_matrix[k,3],
                                             corr = multi_model_condition_matrix[k,4])
        out_gp_sigma_obs <- compare_summary(model = gp_fit, 
                                            variable = sigma_obs_names[[j]], 
                                            true_variable = params[[i]][[j]]$sigma_obs,
                                            n_obs = N, 
                                            dims = dims[j], 
                                            variable_class = 'sigma_obs', 
                                            sim_id = i,
                                            deriv = multi_model_condition_matrix[k,1],
                                            scale = multi_model_condition_matrix[k,2],
                                            vary = multi_model_condition_matrix[k,3],
                                            corr = multi_model_condition_matrix[k,4])
        out_gp_sigma_grad <- compare_summary(model = gp_fit, 
                                             variable = sigma_grad_names[[j]], 
                                             true_variable = params[[i]][[j]]$sigma_grad,
                                             n_obs = N, 
                                             dims = dims[j], 
                                             variable_class = 'sigma_grad', 
                                             sim_id = i,
                                             deriv = multi_model_condition_matrix[k,1],
                                             scale = multi_model_condition_matrix[k,2],
                                             vary = multi_model_condition_matrix[k,3],
                                             corr = multi_model_condition_matrix[k,4])
        out_gp2[[i]][[j]][[k]] <- rbind(out_gp_x, out_gp_rho, 
                                        out_gp_alpha_obs, out_gp_alpha_grad,
                                        out_gp_sigma_obs, out_gp_sigma_grad)
      }
      out_gp1[[i]][[j]] <- rbindlist(out_gp2[[i]][[j]])
    }
    out_gp[[i]] <- rbindlist(out_gp1[[i]])
    
    saveRDS(out_gp[[i]], paste0('temp_summary_trial','_',i,'.rds'))
    return(out_gp[[i]])
  }
}
compare_table <- rbindlist(model_fit)
## comparison table edit
compare_table$model_id = paste0(compare_table$sim_id, '_', 
                                compare_table$d, '_',
                                compare_table$deriv,
                                compare_table$scale,
                                compare_table$vary,
                                compare_table$corr)
saveRDS(compare_table, 'SimDataOutput.rds')