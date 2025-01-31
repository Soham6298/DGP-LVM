# We showcase DGP-LVM for a reduced single-cell RNA sequencing data
# The data cellcycle is obtained from Cytopath (https://doi.org/10.1016/j.crmeth.2022.100359)
# We provide recovery of latent input (pseudotime estimation) for DGP-LVM with Squared Exponential and Matern 3/2 covariance functions

#libraries
library(rstan)
library(bayesplot)
library(ggplot2)
library(posterior)
library(dplyr)
library(brms)
library(grid)
library(gridExtra)
# Source functions
source('DgplvmSimFns.R')
stanmodel <- stan_model('DerivGPmodels.stan')
# Read data
cellcycle <- readRDS('DGPLVMcasestudydata.Rdata')
# Data cleanup
cellcycle[cellcycle==0] <- NA
cellcycle <- na.omit(cellcycle)
ngroups <- seq(0.1, 1, by = 0.1)
cellcycle$cellhrsgroup <- cut(cellcycle$cell_cycle_hrs,
                              breaks = ngroups)
# Sub-sample data to manage computation times
set.seed(980)
cellcycle_strat_sample <- cellcycle %>%
  group_by(cellhrsgroup) %>%
  sample_n(size = 2)
# Arrange in concatenated format
cellcycleexp <- cellcycle_strat_sample[,3:14]
gene_names <- colnames(cellcycleexp)
cellcyclevelo <- cellcycle_strat_sample[15:26]
colnames(cellcyclevelo) <- gene_names
# data for model
y <- as.data.frame(rbind(cellcycleexp, cellcyclevelo))
t <- cellcycle_strat_sample$cell_cycle_hrs
#mean_rho <- mean(as.numeric(dist(t)))
#y <- y[,1:2]
#beta <- mean_rho * (5-1)
# Prepare model data
# set for prior measurement SD for latent x
s_x <- 0.03   
# Number of cells (sample size)
N <- nrow(y)/2
# Number of genes (output dimensions)
D <- ncol(y)
# Indicator for original and derivative parts
#deriv <- deriv <- c(rep(0, length(t)), rep(1, length(t)))
# model fitting
mean_obs <- rep(NA, D)
sd_obs <- rep(NA, D)
mean_grad <- rep(NA, D)
sd_grad <- rep(NA, D)
for(j in 1:D) {
    mean_obs[j] <- mean(y[1:N,j])
    sd_obs[j] <- sd(y[1:N,j])
    mean_grad[j] <- mean(y[(N+1):(2*N),j])
    sd_grad[j] <- sd(y[(N+1):(2*N),j])
}
# Model specifications
rho_prior_model <- 0 # 0 = normal; 1 = invgamma;
lambda_gene <- rep(NA, D)
for(i in 1:D){
  lambda_gene[i] <- mean(abs(y[1:N,i])%/%abs(y[(N+1):(2*N),i]))
}
lambda <- mean(lambda_gene)
# Specify priors for the hyperparameters (check if they differ from the simulation setting)
marginal_sd_params_grad_model <- c(1, 0.25)  # alpha ~ N( , )
error_sd_params_grad_model <- c(0.5, 0.25)     # sigma ~ N( , )
marginal_sd_params_model <- marginal_sd_params_grad_model * lambda
error_sd_params_model <- error_sd_params_grad_model * lambda
ls_params_model_se <- c(0.4, 0.1) # rho ~ N( , ) or InvGamma( , )
ls_params_model_m32 <- c(0.8, 0.1) # rho ~ N( , ) or InvGamma( , )
ls_params_model_m52 <- c(0.6, 0.1) # rho ~ N( , ) or InvGamma( , )
latent_model <- 1 # 1 for latent inputs, 0 for manifest
delta_model <- 1e-6  # set higher for the derivative model if needed
adapt_delta_model <- 0.95
is_corr <- 1
# Model fit
# SE
dgp_se_fit <- gp_model(stanmodel,
                   n_obs = N, 
                   dims = D, 
                   outputs = y, 
                   inputs = t, 
                   latent_sd = s_x, 
                   latent_inputs = latent_model, 
                   covfn = 0, #covfn_model <- 0  # 0 - SE; 1 = Matern3/2; 2 = Matern 5/2
                   rho_prior = rho_prior_model,
                   ls_param = ls_params_model_se, 
                   msd_param = marginal_sd_params_model, 
                   msd_param_grad = marginal_sd_params_grad_model,
                   esd_param = error_sd_params_model, 
                   esd_param_grad = error_sd_params_grad_model,
                   mean_obs = mean_obs,
                   sd_obs = sd_obs,
                   mean_grad = mean_grad,
                   sd_grad = sd_grad,
                   nugget = delta_model, 
                   is_deriv = 1, 
                   is_scale = 1, 
                   is_vary = 1, 
                   is_corr = 1, 
                   iter = 2000, 
                   warmup = 1000, 
                   chains = 2, 
                   cores = 2, 
                   init = 0, 
                   adapt_delta = adapt_delta_model)
# Matern 3/2
dgp_m32_fit <- gp_model(stanmodel,
                       n_obs = N, 
                       dims = D, 
                       outputs = y, 
                       inputs = t, 
                       latent_sd = s_x, 
                       latent_inputs = latent_model, 
                       covfn = 1, #covfn_model <- 0  # 0 - SE; 1 = Matern3/2; 2 = Matern 5/2
                       rho_prior = rho_prior_model,
                       ls_param = ls_params_model_m32, 
                       msd_param = marginal_sd_params_model, 
                       msd_param_grad = marginal_sd_params_grad_model,
                       esd_param = error_sd_params_model, 
                       esd_param_grad = error_sd_params_grad_model,
                       mean_obs = mean_obs,
                       sd_obs = sd_obs,
                       mean_grad = mean_grad,
                       sd_grad = sd_grad,
                       nugget = delta_model, 
                       is_deriv = 1, 
                       is_scale = 1, 
                       is_vary = 1, 
                       is_corr = 1, 
                       iter = 2000, 
                       warmup = 1000, 
                       chains = 2, 
                       cores = 2, 
                       init = 0, 
                       adapt_delta = adapt_delta_model)
# Matern 5/2
dgp_m52_fit <- gp_model(stanmodel,
                       n_obs = N, 
                       dims = D, 
                       outputs = y, 
                       inputs = t, 
                       latent_sd = s_x, 
                       latent_inputs = latent_model, 
                       covfn = 2, #covfn_model <- 0  # 0 - SE; 1 = Matern3/2; 2 = Matern 5/2
                       rho_prior = rho_prior_model,
                       ls_param = ls_params_model_m52, 
                       msd_param = marginal_sd_params_model, 
                       msd_param_grad = marginal_sd_params_grad_model,
                       esd_param = error_sd_params_model, 
                       esd_param_grad = error_sd_params_grad_model,
                       mean_obs = mean_obs,
                       sd_obs = sd_obs,
                       mean_grad = mean_grad,
                       sd_grad = sd_grad,
                       nugget = delta_model, 
                       is_deriv = 1, 
                       is_scale = 1, 
                       is_vary = 1, 
                       is_corr = 1, 
                       iter = 2000, 
                       warmup = 1000, 
                       chains = 2, 
                       cores = 2, 
                       init = 0, 
                       adapt_delta = adapt_delta_model)

print(dgp_se_fit, pars = c('rho','alpha_obs', 'alpha_grad', 'sigma_obs', 'sigma_grad', 'x'))
print(dgp_m32_fit, pars = c('rho','alpha_obs', 'alpha_grad', 'sigma_obs', 'sigma_grad', 'x'))
print(dgp_m52_fit, pars = c('rho','alpha_obs', 'alpha_grad', 'sigma_obs', 'sigma_grad', 'x'))

## Model summary
x_names <- sprintf('x[%s]', seq(1:20))
## SE
dgp_se_draws <- as_draws_df(dgp_se_fit)
dgp_se_draws_x <- subset_draws(dgp_se_draws, variable = x_names, chain = 1)
dgp_se_summary <- summarise_draws(dgp_se_draws_x)
dgp_se_summary$obs_t <- t 
cor.test(dgp_se_summary$mean, dgp_se_summary$obs_t, method = 'spearman')
dgp_se_plot <- ggplot(dgp_se_summary, aes(x = obs_t, y = mean - obs_t)) + 
  theme_bw(base_size = 35, base_family = 'Times') +
  geom_errorbarh(aes(xmin = obs_t-(1.96*s_x), xmax = obs_t+(1.96*s_x), color = 'Prior'), linewidth = 1) +
  geom_errorbar(aes(ymin = q5 - obs_t, ymax = q95 - obs_t, color = 'Posterior'), linewidth = 1) + 
  geom_point(size = 3) +
  geom_hline(yintercept = 0) +
  ylim(c(-0.15, 0.15)) +
  labs(x = 'Cell hours', y = 'Pseudotime - Cell hours', color = 'CI') +
  scale_color_manual(values = c("#0072B2", "#E69F00" )) + ggtitle('(a)')
## Matern 3/2
dgp_m32_draws <- as_draws_df(dgp_m32_fit)
dgp_m32_draws_x <- subset_draws(dgp_m32_draws, variable = x_names, chain = 1)
dgp_m32_summary <- summarise_draws(dgp_m32_draws_x)
dgp_m32_summary$obs_t <- t 
dgp_m32_plot <- ggplot(dgp_m32_summary, aes(x = obs_t, y = mean - obs_t)) +
  theme_bw(base_size = 35, base_family = 'Times') +
  geom_errorbarh(aes(xmin = obs_t-(1.96*s_x), xmax = obs_t+(1.96*s_x), color = 'Prior'), linewidth = 1) +
  geom_errorbar(aes(ymin = q5 - obs_t, ymax = q95 - obs_t, color = 'Posterior'), linewidth = 1) + 
  geom_point(size = 3) +
  geom_hline(yintercept = 0) +
  ylim(c(-0.13, 0.13)) +
  labs(x = 'Cell hours', y = 'Pseudotime - Cell hours', color = 'CI') +
  scale_color_manual(values = c("#0072B2", "#E69F00" )) + ggtitle('(c)')
## Matern 5/2
dgp_m52_draws <- as_draws_df(dgp_m52_fit)
dgp_m52_draws_x <- subset_draws(dgp_m52_draws, variable = x_names, chain = 1)
dgp_m52_summary <- summarise_draws(dgp_m52_draws_x)
dgp_m52_summary$obs_t <- t 
dgp_m52_plot <- ggplot(dgp_m52_summary, aes(x = obs_t, y = mean - obs_t)) +
  theme_bw(base_size = 35, base_family = 'Times') +
  geom_errorbarh(aes(xmin = obs_t-(1.96*s_x), xmax = obs_t+(1.96*s_x), color = 'Prior'), linewidth = 1) +
  geom_errorbar(aes(ymin = q5 - obs_t, ymax = q95 - obs_t, color = 'Posterior'), linewidth = 1) + 
  geom_point(size = 3) +
  geom_hline(yintercept = 0) +
  ylim(c(-0.15, 0.15)) +
  labs(x = 'Cell hours', y = 'Pseudotime - Cell hours', color = 'CI') +
  scale_color_manual(values = c("#0072B2", "#E69F00" )) + ggtitle('(b)')
case_study_plot <- (dgp_se_plot + dgp_m52_plot + plot_layout(axes = 'collect', guides = 'collect') &
                      theme(legend.position = 'bottom'))
ggsave('dgp_case_study_plot.pdf',
       case_study_plot,
       dpi = 300,
       width = 30,
       height = 20,
       units = 'cm')
## Parameters
# SE
dgp_se_draws_pars <- subset_draws(dgp_se_draws, variable = c('rho', 'alpha_obs', 'alpha_grad',
                                                             'sigma_obs', 'sigma_grad'), chain = 1)
dgp_se_summary_pars <- summarise_draws(dgp_se_draws_pars)
var_name <- c(rep('rho', D), rep('alpha_obs', D), rep('alpha_grad', D),
              rep('sigma_obs', D), rep('sigma_grad', D))
output_dim <- rep(seq(1:D), 5) # Based on number of hyperparameters
# Design the dataframe for plot
dgp_se_summary_pars$var_name <- var_name
dgp_se_summary_pars$output_dim <- as.factor(output_dim)
dgp_se_summary_rho <- dgp_se_summary_pars[var_name=='rho',]
dgp_se_summary_rho$ftitle <- 'rho'
dgp_se_summary_alpha <- dgp_se_summary_pars[var_name=='alpha_obs'|var_name=='alpha_grad',]
dgp_se_summary_alpha$var_name <- as.factor(dgp_se_summary_alpha$var_name)
dgp_se_summary_alpha$ftitle <- 'alpha'
levels(dgp_se_summary_alpha$var_name) <- c("y'",'y')
dgp_se_summary_sigma <- dgp_se_summary_pars[var_name=='sigma_obs'|var_name=='sigma_grad',]
dgp_se_summary_sigma$var_name <- as.factor(dgp_se_summary_sigma$var_name)
dgp_se_summary_sigma$ftitle <- 'sigma'
levels(dgp_se_summary_sigma$var_name) <- c("y'",'y')
# Matern 3/2
dgp_m32_draws_pars <- subset_draws(dgp_m32_draws, variable = c('rho', 'alpha_obs', 'alpha_grad',
                                                               'sigma_obs', 'sigma_grad'), chain = 1)
dgp_m32_summary_pars <- summarise_draws(dgp_m32_draws_pars)
var_name <- c(rep('rho', D), rep('alpha_obs', D), rep('alpha_grad', D),
              rep('sigma_obs', D), rep('sigma_grad', D))
output_dim <- rep(seq(1:D), 5) # Based on number of hyperparameters
dgp_m32_summary_pars$var_name <- var_name
dgp_m32_summary_pars$output_dim <- as.factor(output_dim)
dgp_m32_summary_rho <- dgp_m32_summary_pars[var_name=='rho',]
dgp_m32_summary_rho$ftitle <- 'rho'
dgp_m32_summary_alpha <- dgp_m32_summary_pars[var_name=='alpha_obs'|var_name=='alpha_grad',]
dgp_m32_summary_alpha$var_name <- as.factor(dgp_m32_summary_alpha$var_name)
dgp_m32_summary_alpha$ftitle <- 'alpha'
levels(dgp_m32_summary_alpha$var_name) <- c("y'",'y')
dgp_m32_summary_sigma <- dgp_m32_summary_pars[var_name=='sigma_obs'|var_name=='sigma_grad',]
dgp_m32_summary_sigma$var_name <- as.factor(dgp_m32_summary_sigma$var_name)
dgp_m32_summary_sigma$ftitle <- 'sigma'
levels(dgp_m32_summary_sigma$var_name) <- c("y'",'y')

# Matern 5/2
dgp_m52_draws_pars <- subset_draws(dgp_m52_draws, variable = c('rho', 'alpha_obs', 'alpha_grad',
                                                               'sigma_obs', 'sigma_grad'), chain = 1)
dgp_m52_summary_pars <- summarise_draws(dgp_m52_draws_pars)
var_name <- c(rep('rho', D), rep('alpha_obs', D), rep('alpha_grad', D),
              rep('sigma_obs', D), rep('sigma_grad', D))
output_dim <- rep(seq(1:D), 5) # Based on number of hyperparameters
dgp_m52_summary_pars$var_name <- var_name
dgp_m52_summary_pars$output_dim <- as.factor(output_dim)
dgp_m52_summary_rho <- dgp_m52_summary_pars[var_name=='rho',]
dgp_m52_summary_rho$ftitle <- 'rho'
dgp_m52_summary_alpha <- dgp_m52_summary_pars[var_name=='alpha_obs'|var_name=='alpha_grad',]
dgp_m52_summary_alpha$var_name <- as.factor(dgp_m52_summary_alpha$var_name)
dgp_m52_summary_alpha$ftitle <- 'alpha'
levels(dgp_m52_summary_alpha$var_name) <- c("y'",'y')
dgp_m52_summary_sigma <- dgp_m52_summary_pars[var_name=='sigma_obs'|var_name=='sigma_grad',]
dgp_m52_summary_sigma$var_name <- as.factor(dgp_m52_summary_sigma$var_name)
dgp_m52_summary_sigma$ftitle <- 'sigma'
levels(dgp_m52_summary_sigma$var_name) <- c("y'",'y')

## SE plot
dgp_se_rho_plot <- ggplot(dgp_se_summary_rho, aes(x = output_dim, y = mean)) +   
  theme_bw(base_size = 65, base_family = 'Times') +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = q5, ymax = q95), linewidth = 1) + 
  labs(x = 'Output dimensions', y = 'Posterior mean') +
  theme(axis.title = element_blank()) +
  facet_wrap(~ftitle) +
  scale_y_log10() +
  ggtitle('(a)')
dgp_se_alpha_plot <- ggplot(dgp_se_summary_alpha, aes(x = output_dim, y = mean, colour = var_name)) +   
  theme_bw(base_size = 65, base_family = 'Times') +
  geom_point(size = 3, position = position_dodge(0.7)) +
  geom_errorbar(aes(ymin = q5, ymax = q95), linewidth = 1, position = position_dodge(0.7)) + 
  labs(x = 'Output dimensions', y = 'Posterior mean') + 
  facet_wrap(~ftitle) +
  scale_y_log10() +
  theme(legend.title = element_blank(), axis.title = element_blank()) +
  scale_colour_manual(values = c("#0072B2", "#E69F00"))
dgp_se_sigma_plot <- ggplot(dgp_se_summary_sigma, aes(x = output_dim, y = mean, colour = var_name)) +   
  theme_bw(base_size = 65, base_family = 'Times') +
  geom_point(size = 3, position = position_dodge(0.7)) +
  geom_errorbar(aes(ymin = q5, ymax = q95), linewidth = 1, position = position_dodge(0.7)) + 
  labs(x = 'Output dimensions', y = 'Posterior mean') + 
  facet_wrap(~ftitle) +
  scale_y_log10() +
  theme(legend.title = element_blank(), axis.title = element_blank()) +
  scale_colour_manual(values = c("#0072B2", "#E69F00"))

## Matern 3/2 plot
dgp_m32_rho_plot <- ggplot(dgp_m32_summary_rho, aes(x = output_dim, y = mean)) +   
  theme_bw(base_size = 65, base_family = 'Times') +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = q5, ymax = q95), linewidth = 1) + 
  labs(x = 'Output dimensions', y = 'Posterior mean') +
  facet_wrap(~ftitle) +
  scale_y_log10() +
  theme(axis.title = element_blank()) +
  ggtitle('(c)')
dgp_m32_alpha_plot <- ggplot(dgp_m32_summary_alpha, aes(x = output_dim, y = mean, colour = var_name)) +   
  theme_bw(base_size = 65, base_family = 'Times') +
  geom_point(size = 3, position = position_dodge(0.7)) +
  geom_errorbar(aes(ymin = q5, ymax = q95), linewidth = 1, position = position_dodge(0.7)) + 
  labs(x = 'Output dimensions', y = 'Posterior mean') + 
  facet_wrap(~ftitle) +
  scale_y_log10() +
  theme(legend.title = element_blank(), axis.title = element_blank()) +
  scale_colour_manual(values = c("#0072B2", "#E69F00"))
dgp_m32_sigma_plot <- ggplot(dgp_m32_summary_sigma, aes(x = output_dim, y = mean, colour = var_name)) +   
  theme_bw(base_size = 65, base_family = 'Times') +
  geom_point(size = 3, position = position_dodge(0.7)) +
  geom_errorbar(aes(ymin = q5, ymax = q95), linewidth = 1, position = position_dodge(0.7)) + 
  labs(x = 'Output dimensions', y = 'Posterior mean') + 
  facet_wrap(~ftitle) +
  scale_y_log10() +
  theme(legend.title = element_blank(), axis.title = element_blank()) +
  scale_colour_manual(values = c("#0072B2", "#E69F00"))
## Matern 5/2 plot
dgp_m52_rho_plot <- ggplot(dgp_m52_summary_rho, aes(x = output_dim, y = mean)) +   
  theme_bw(base_size = 65, base_family = 'Times') +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = q5, ymax = q95), linewidth = 1) + 
  labs(x = 'Output dimensions', y = 'Posterior mean') +
  facet_wrap(~ftitle) +
  scale_y_log10() +
  theme(axis.title = element_blank()) +
  ggtitle('(b)')
dgp_m52_alpha_plot <- ggplot(dgp_m52_summary_alpha, aes(x = output_dim, y = mean, colour = var_name)) +   
  theme_bw(base_size = 65, base_family = 'Times') +
  geom_point(size = 3, position = position_dodge(0.7)) +
  geom_errorbar(aes(ymin = q5, ymax = q95), linewidth = 1, position = position_dodge(0.7)) + 
  labs(x = 'Output dimensions', y = 'Posterior mean') + 
  facet_wrap(~ftitle) +
  scale_y_log10() +
  theme(legend.title = element_blank(), axis.title = element_blank()) +
  scale_colour_manual(values = c("#0072B2", "#E69F00"))
dgp_m52_sigma_plot <- ggplot(dgp_m52_summary_sigma, aes(x = output_dim, y = mean, colour = var_name)) +   
  theme_bw(base_size = 65, base_family = 'Times') +
  geom_point(size = 3, position = position_dodge(0.7)) +
  geom_errorbar(aes(ymin = q5, ymax = q95), linewidth = 1, position = position_dodge(0.7)) + 
  labs(x = 'Output dimensions', y = 'Posterior mean') + 
  facet_wrap(~ftitle) +
  scale_y_log10() +
  theme(legend.title = element_blank(), axis.title = element_blank()) +
  scale_colour_manual(values = c("#0072B2", "#E69F00"))
## Combine
dgp_cellcycle_pars <- (dgp_se_rho_plot + dgp_se_alpha_plot + dgp_se_sigma_plot) /
  (dgp_m52_rho_plot + dgp_m52_alpha_plot + dgp_m52_sigma_plot) +
 # (dgp_m32_rho_plot + dgp_m32_alpha_plot + dgp_m32_sigma_plot) + 
  plot_layout(axis_titles = 'collect', guides = 'collect') &
  theme(axis.text = element_text(size = 55))
dgp_cellcycle_pars
gt <- patchwork::patchworkGrob(dgp_cellcycle_pars)
g <- arrangeGrob(gt, left = textGrob("Estimate", rot = 90, gp = gpar(fontsize=65, fontfamily='Times')), 
                 bottom = textGrob("Gene", vjust = -0.2, gp = gpar(fontsize=65, fontfamily='Times')))
ggsave('dgp_case_study_pars.pdf',
       g,
       dpi = 300,
       width = 100,
       height = 50,
       units = 'cm')

