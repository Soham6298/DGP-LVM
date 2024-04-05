// DGP-LVM with Matern 3/2 (M32) covariance function 
// GP Paramters:
// GP length scale: rho
// GP marginal SD: alpha_obs (corresponding to the original part of GP); alpha_grad (corresponding to the derivative part of GP)
// Error SD: sigma_obs (corresponding to original output); sigma_grad (corresponding to derivative part)

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
  matrix deriv_m32(vector x, int[] derivative, real alpha_obs, real alpha_grad, real rho, real delta) {
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
  // base covariance fn
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
  int<lower=0, upper=1> derivative[N];
  
  // Prior measurement SD (assumed to be known)
  real<lower=0> s;
  // Prior SD for alphas and sigmas
  real<lower=0> sparams_prior_sd;
  // Prior for latent x
  vector[M] t;
  // Model conditions (specify which modifications to use ranging from standard GPs to DGP-LVM)
  int<lower=0, upper=1> is_deriv; // 0 = no derivative; 1 = derivative model
  int<lower=0, upper=1> is_scale; // 0 = same scale; 1 = scaled
  int<lower=0, upper=1> is_vary; // 0 = constant rho and alpha for each dims; 1 = varying params
  int<lower=0, upper=1> is_corr; // 0 = no correlation b/w dims; 1 = correlated dims
}

transformed data {
// add constant to the diagonal of the covariance matrix for computational stability
  real delta = 1e-6;
}

parameters {
  // Latent x
  // add bounds on x if known for theoretical reasons (e.g. lower=0)
  vector[M] x;
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
}

transformed parameters {
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
    for (k in 1:D) {
      matrix[N, N] K = deriv_se(x2, derivative, alpha_obs[k], alpha_grad[k], rho[k], delta);
      f[, k] = K * eta[, k];
    }
  // For correlated outputs
    f = f * L_omega';
  } else {
  // Computing covariance matrix for standard GP
    for (k in 1:D) {
      matrix[M, M] K = se(x, alpha_obs[k], rho[k], delta);
      f[, k] = K * eta[, k];
    }
  // For correlated outputs
    f = f * L_omega';
  }
}

model {
  // Priors
  rho_temp ~ inv_gamma(5, 5);
  alpha_obs_temp ~ normal(0, sparams_prior_sd);
  alpha_grad_temp ~ normal(0, sparams_prior_sd);
  sigma_obs_temp ~ normal(0, sparams_prior_sd);
  sigma_grad_temp ~ normal(0, sparams_prior_sd);
  L_omega_temp ~ lkj_corr_cholesky(3);
  for (k in 1:D) {
    to_vector(eta[, k]) ~ std_normal();
  }
  // Likelihood
  t ~ normal(x, s); //latent prior
  if (is_deriv == 1) {
    for (k in 1:D) {
      y[1:M, k] ~ normal(f[1:M, k], sigma_obs[k]);
      y[(M+1):N, k] ~ normal(f[(M+1):N, k], sigma_grad[k]);
    }
  } else {
    // here the second half of y (derivative info) remains unused
    for (k in 1:D) {
      y[1:M, k] ~ normal(f[1:M, k], sigma_obs[k]);
    }
  }
}
