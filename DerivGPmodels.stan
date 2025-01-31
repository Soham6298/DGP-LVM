functions {
  // for repeating inputs
  vector repeat_vector(vector x, int reps) {
    int N = rows(x);
    array[N * reps] int index;
    for (i in 1:reps) {
      for (j in 1:N) {       
        index[(i-1) * N + j] = j;
      }
    }    
    return x[index];
  }
  // derivative covariance fn
  
  // Deriv SE
  matrix deriv_se(vector x, array[]int derivative, real alpha_obs, real alpha_grad, real rho, real delta) {
    int N = rows(x);
    matrix[N, N] K;
    real sq_rho = square(rho);
    real rho4 = pow(rho, 4);
    real sq_alpha_obs = pow(alpha_obs, 2);
    real sq_alpha_grad = pow(alpha_grad, 2);
    real r = -inv(2 * sq_rho);
    
    for (i in 1:N) {
      if (derivative[i] == 0) {
        K[i, i] = sq_alpha_obs + delta;
      } else if (derivative[i] == 1) {
        K[i, i] = (sq_alpha_grad / sq_rho) + delta;
      }
      for (j in (i + 1):N) {
        if(derivative[i] == 0 && derivative[j] == 0) {
          K[i, j] = exp(r * square(x[i] - x[j])) * sq_alpha_obs;
        } else if(derivative[i] == 0 && derivative[j] == 1) {
          K[i, j] = exp(r * square(x[i] - x[j])) * 
            (x[i] - x[j]) * ((alpha_obs * alpha_grad) / sq_rho);
        } else if(derivative[i] == 1 && derivative[j] == 0) {
          K[i, j] = exp(r * square(x[i] - x[j])) * 
            (x[j] - x[i]) * ((alpha_grad * alpha_obs) / sq_rho);
        } else if(derivative[i] == 1 && derivative[j] == 1) {
          K[i, j] = exp(r * square(x[i] - x[j])) * 
            (sq_rho - square(x[i] - x[j])) * (sq_alpha_grad / rho4);
        }
        K[j, i] = K[i, j];
      }
    }
    return cholesky_decompose(K);
  }

  // Deriv Matern 3/2
  matrix deriv_m32(vector x, array[]int derivative, real alpha_obs, real alpha_grad, real rho, real delta) {
    int N = rows(x);
    matrix[N, N] K;
    real sq_rho = square(rho);
    real sq_alpha_obs = pow(alpha_obs, 2);
    real sq_alpha_grad = pow(alpha_grad, 2);
    real r = -inv(rho);
    
    for (i in 1:N) {
      if (derivative[i] == 0) {
        K[i, i] = sq_alpha_obs + delta;
      } else if (derivative[i] == 1) {
        K[i, i] = (3 * sq_alpha_grad / sq_rho) + delta;
      }
      for (j in (i + 1):N) {
        if(derivative[i] == 0 && derivative[j] == 0) {
          K[i, j] = (1 - (r * sqrt(3) * abs(x[i] - x[j]))) * 
                    exp(r * sqrt(3) * abs(x[i] - x[j])) * sq_alpha_obs;
        } else if(derivative[i] == 0 && derivative[j] == 1) {
          K[i, j] = (3 * (x[i] - x[j]) / sq_rho) * 
                    exp(r * sqrt(3) * abs(x[i] - x[j])) * alpha_obs * alpha_grad;
        } else if(derivative[i] == 1 && derivative[j] == 0) {
          K[i, j] = (3 * (x[j] - x[i]) / sq_rho) * 
                    exp(r * sqrt(3) * abs(x[i] - x[j])) * alpha_obs * alpha_grad;
        } else if(derivative[i] == 1 && derivative[j] == 1) {
          K[i, j] = (1 + (r * sqrt(3) * abs(x[i] - x[j]))) * 
                    exp(r * sqrt(3) * abs(x[i] - x[j])) * sq_alpha_grad * (3 / sq_rho);
        }
        K[j, i] = K[i, j];
      }
    }
    return cholesky_decompose(K);
  }
  
  // Deriv Matern 5/2
  matrix deriv_m52(vector x, array[]int derivative, real alpha_obs, real alpha_grad, real rho, real delta) {
    int N = rows(x);
    matrix[N, N] K;
    real sq_rho = square(rho);
    real rho4 = pow(rho, 4);
    real sq_alpha_obs = pow(alpha_obs, 2);
    real sq_alpha_grad = pow(alpha_grad, 2);
    real r = -inv(2 * sq_rho);
    
    for (i in 1:N) {
      if (derivative[i] == 0) {
        K[i, i] = sq_alpha_obs + delta;
      } else if (derivative[i] == 1) {
        K[i, i] = (5 * sq_alpha_grad / (3 * sq_rho)) + delta;
      }
      for (j in (i + 1):N) {
        if(derivative[i] == 0 && derivative[j] == 0) {
          K[i, j] = (1 + (sqrt(5) * abs(x[i] - x[j])/rho) + 
                    ((5 * square(x[i] - x[j]))/ (3 * sq_rho))) * 
                    exp(- sqrt(5) * abs(x[i] - x[j])/rho) * sq_alpha_obs;
        } else if(derivative[i] == 0 && derivative[j] == 1) {
          K[i, j] = ((5 * (x[i] - x[j]))/ (3 * sq_rho)) * 
                    (1 + (sqrt(5) * abs(x[i] - x[j])/rho)) * 
                    exp(- sqrt(5) * abs(x[i] - x[j])/rho) * alpha_obs * alpha_grad;
        } else if(derivative[i] == 1 && derivative[j] == 0) {
          K[i, j] = ((5 * (x[j] - x[i]))/ (3 * sq_rho)) * 
                    (1 + (sqrt(5) * abs(x[i] - x[j])/rho)) * 
                    exp(- sqrt(5) * abs(x[i] - x[j])/rho) * alpha_obs * alpha_grad;
        } else if(derivative[i] == 1 && derivative[j] == 1) {
          K[i, j] = (1 + (sqrt(5) * abs(x[i] - x[j])/rho) - 
                    ((5 * square(x[i] - x[j]))/ sq_rho)) * 
                    exp(- sqrt(5) * abs(x[i] - x[j])/rho) * (5 * sq_alpha_grad / (3 * sq_rho));
        }
        K[j, i] = K[i, j];
      }
    }
    return cholesky_decompose(K);
  }

  // base covariance fn
  
  // SE
    matrix se(vector x, real alpha_obs, real rho, real delta) {
    int N = rows(x);
    matrix[N, N] K;
    real sq_rho = square(rho);
    real rho4 = pow(rho, 4);
    real sq_alpha_obs = pow(alpha_obs, 2);
    real r = -inv(2 * sq_rho);
    
    for(i in 1:N) {
      K[i, i] = sq_alpha_obs + delta;
      for(j in (i + 1):N) {
        K[i, j] = exp(r * square(x[i] - x[j])) * sq_alpha_obs;
        K[j, i] = K[i, j];
      }
    }
    return cholesky_decompose(K);
  }
  
  // Matern 3/2
  matrix m32(vector x, real alpha_obs, real rho, real delta) {
    int N = rows(x);
    matrix[N, N] K;
    real sq_alpha_obs = pow(alpha_obs, 2);
    real r = -inv(rho);
    
    for(i in 1:N) {
      K[i, i] = sq_alpha_obs + delta;
      for(j in (i + 1):N) {
        K[i, j] = (1 - (r * sqrt(3) * abs(x[i] - x[j]))) * 
                    exp(r * sqrt(3) * abs(x[i] - x[j])) * sq_alpha_obs;
        K[j, i] = K[i, j];
      }
    }
    return cholesky_decompose(K);
  }
  
  //Matern 5/2 
   matrix m52(vector x, real alpha_obs, real rho, real delta) {
    int N = rows(x);
    matrix[N, N] K;
    real sq_rho = square(rho);
    real rho4 = pow(rho, 4);
    real sq_alpha_obs = pow(alpha_obs, 2);
    real r = -inv(2 * sq_rho);
    
    for(i in 1:N) {
      K[i, i] = sq_alpha_obs + delta;
      for(j in (i + 1):N) {
        K[i, j] = (1 + (sqrt(5) * abs(x[i] - x[j])/rho) + 
                    ((5 * square(x[i] - x[j]))/ (3 * sq_rho))) * 
                    exp(- sqrt(5) * abs(x[i] - x[j])/rho) * sq_alpha_obs;
        K[j, i] = K[i, j];
      }
    }
    return cholesky_decompose(K);
  }
}

data {
  // Sample size
  int<lower=1> N;
  // Number of unique x values; usually N / 2
  int<lower=1> M;
  // Number of output dimensions
  int<lower=1> D;
  // Output (Original and derivative concatenated for DGP-LVM)
  matrix[N, D] y; 
  // Indicator variable for original and derivative part
  array[N] int<lower=0, upper=1> derivative;
  // Prior measurement SD (assumed to be known)
  real<lower=0> s;
  // Prior for latent x
  vector[M] inputs;
  // Model conditions (specify which modifications to use ranging from standard GPs to DGP-LVM)
  int<lower=0, upper=1> is_deriv; // 0 = no derivative; 1 = derivative model
  int<lower=0, upper=1> is_scale; // 0 = same scale; 1 = scaled
  int<lower=0, upper=1> is_vary; // 0 = constant rho and alpha for each dims; 1 = varying params
  int<lower=0, upper=1> is_corr; // 0 = no correlation b/w dims; 1 = correlated dims
  int<lower=0, upper=2> covfn; //0 = SE; 1 = M3/2; 2 = M5/2
  int<lower=0, upper=1> latent; //0 = obs inputs; 1 = latent inputs
  int<lower=0, upper=1> rho_prior; //0 = normal; 1 = invgamma;
  // prior hyperparams specification (only for two parameter dist families)
  vector[2] ls_param; 
  vector[2] msd_param;
  vector[2] esd_param;
  // Adjust derivative part priors to account for scale = 1
  vector[2] msd_param_grad; 
  vector[2] esd_param_grad;
  // Input data specific mean and sd for the y and y'
  array[D] real mean_obs;
  array[D] real sd_obs;
  array[D] real mean_grad;
  array[D] real sd_grad;
  // covfn diag constant
  real nugget;
}

transformed data {
// add constant to the diagonal of the covariance matrix for computational stability
  real delta = nugget;
}

parameters {
  // GP Length scale parameter (temporary storage)
  vector<lower=0>[is_vary==1 ? D:1] rho_temp;
  // GP marginal SD parameters (temporary storage)
  vector<lower=0>[is_vary==1 ? D:1] alpha_obs_temp;
  vector<lower=0>[is_scale==1 ? D:0] alpha_grad_temp;
  // Gaussian link parameter
  matrix[is_deriv==1 ? N:M, D] eta;
  // Between dimension correlation parameter (temporary storage)
  cholesky_factor_corr[is_corr==1 ? D:0] L_omega_temp;
  // Error SD parameters (temporary storage)
  vector<lower=0>[is_vary==1 ? D:1] sigma_obs_temp;
  vector<lower=0>[is_scale==1 ? D:0] sigma_grad_temp;
  // only for latent variable models
  vector[latent==1 ? M:0] z;
  //vector<lower=0>[latent==1 ? M:0] x_temp;
  array[D] real intercept_obs;
  array[D] real intercept_grad;
}

transformed parameters {
  // Latent x
  // add bounds on x if known for theoretical reasons (e.g. lower=0)
  vector[M] x;
  if (latent == 1) {
	  x = inputs + z * s;
	  //x = x_temp;
	} else{
	  x = inputs;
	}
  // GP f()
  matrix[is_deriv==1 ? N:M, D] f;
  // Model condition adjustments (temporary to actual parameters)
  vector<lower=0>[D] rho;
  vector<lower=0>[D] alpha_obs;
  vector<lower=0>[D] alpha_grad;
  vector<lower=0>[D] sigma_obs;
  vector<lower=0>[D] sigma_grad;
  cholesky_factor_corr[D] L_omega;
  // if params are constant for dims, it will be repeated
  if (is_vary == 1) {
    rho = rho_temp;
    alpha_obs = alpha_obs_temp;
    sigma_obs = sigma_obs_temp;
  } else {
    for (k in 1:D) {
      rho[k] = rho_temp[1];
      alpha_obs[k] = alpha_obs_temp[1];
      sigma_obs[k] = sigma_obs_temp[1];
    }
  }
  // if deriv scale is 1x, the grad params will be repeated as obs
  if (is_scale == 1) {
    alpha_grad = alpha_grad_temp;
    sigma_grad = sigma_grad_temp;
  } else {
    alpha_grad = alpha_obs;
    sigma_grad = sigma_obs;
  }
  // if no correlation b/w dims, corr mat will be replaced by identity matrix
  if (is_corr==0) {
    L_omega = diag_matrix(rep_vector(1, D));
  } else {
    L_omega = L_omega_temp;
  }
  // Computing covariance matrix for Derivative GP
  if (is_deriv == 1) {
    vector[N] x2 = repeat_vector(x, 2);
    if (covfn == 0) {
      for (k in 1:D) {
        matrix[N, N] K = deriv_se(x2, derivative, alpha_obs[k], alpha_grad[k], rho[k], delta);
        f[, k] = K * eta[, k];
      } 
    } else if (covfn == 1) {
      for (k in 1:D) {
        matrix[N, N] K = deriv_m32(x2, derivative, alpha_obs[k], alpha_grad[k], rho[k], delta);
        f[, k] = K * eta[, k];
      } 
    } else if (covfn == 2) {
      for (k in 1:D) {
        matrix[N, N] K = deriv_m52(x2, derivative, alpha_obs[k], alpha_grad[k], rho[k], delta);
        f[, k] = K * eta[, k];
      } 
    }
    // For correlated outputs
    f = f * L_omega';
  } else {
  // Computing covariance matrix for standard GP
    if (covfn == 0) {
      for (k in 1:D) {
        matrix[M, M] K = se(x, alpha_obs[k], rho[k], delta);
        f[, k] = K * eta[, k];
      } 
    } else if (covfn == 1) {
       for (k in 1:D) {
        matrix[M, M] K = m32(x, alpha_obs[k], rho[k], delta);
        f[, k] = K * eta[, k];
      } 
    } else if (covfn == 2) {
       for (k in 1:D) {
        matrix[M, M] K = m52(x, alpha_obs[k], rho[k], delta);
        f[, k] = K * eta[, k];
      } 
    }
    // For correlated outputs
    f = f * L_omega';
  }
}

model {
  // Priors
  if (rho_prior == 0) {
    rho_temp ~ normal(ls_param[1], ls_param[2]); 
  } else if (rho_prior == 1) {
    rho_temp ~ inv_gamma(ls_param[1], ls_param[2]);
  }
  alpha_obs_temp ~ normal(msd_param[1], msd_param[2]);
  sigma_obs_temp ~ normal(esd_param[1], esd_param[2]);
  if (is_scale == 1) {
    alpha_grad_temp ~ normal(msd_param_grad[1], msd_param_grad[2]);
    sigma_grad_temp ~ normal(esd_param_grad[1], esd_param_grad[2]);
  } else {
    alpha_grad_temp ~ normal(msd_param[1], msd_param[2]);
    sigma_grad_temp ~ normal(esd_param[1], esd_param[2]);
  }
  L_omega_temp ~ lkj_corr_cholesky(1);
  for (k in 1:D) {
    to_vector(eta[, k]) ~ std_normal();
  }
  for (k in 1:D) {
    intercept_obs[k] ~ normal(mean_obs[k], sd_obs[k]);
    intercept_grad[k] ~ normal(mean_grad[k], sd_grad[k]);
  }
  // Likelihood
  if (latent) {
    z ~ std_normal();
    //inputs ~ normal(x_temp, s);
	}
  if (is_deriv == 1) {
    for (k in 1:D) {
      y[1:M, k] ~ normal(intercept_obs[k] + f[1:M, k], sigma_obs[k]);
      y[(M+1):N, k] ~ normal(intercept_grad[k] + f[(M+1):N, k], sigma_grad[k]);
    }
  } else {
    // here the second half of y (derivative info) remains unused
    for (k in 1:D) {
      y[1:M, k] ~ normal(intercept_obs[k] + f[1:M, k], sigma_obs[k]);
    }
  }
}
