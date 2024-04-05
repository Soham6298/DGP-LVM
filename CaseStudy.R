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
# Read data
#setwd("") # Set directory
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
y <- rbind(cellcycleexp, cellcyclevelo)
t <- cellcycle_strat_sample$cell_cycle_hrs
mean_rho <- mean(as.numeric(dist(t)))
beta <- mean_rho * (5-1)
# Functions for model summary
abs_bias_draws <- function(theta_hat, theta) {
  abs(mean(theta_hat) - theta)
}

rmse_draws <- function(theta_hat, theta) {
  sqrt(mean((theta_hat - theta)^2))
}
# Prepare model data
# set for prior measurement SD for latent x
s_x <- 0.03   
# Number of cells (sample size)
N <- nrow(y)
# Number of genes (output dimensions)
D <- ncol(y)
# Indicator for original and derivative parts
deriv <- deriv <- c(rep(0, length(t)), rep(1, length(t)))
# Prior SD for GP marginal and error SD
sparams <- apply(y, 2, sd)
sparams_prior_sd <- mean(sparams)
# Data for the model
dgp_data <- list(
  N = N,
  D = D,
  M = N / 2,
  y = y,
  t = t,
  s = s_x,
  sparams_prior_sd = sparams_prior_sd,
  derivative = deriv,
  is_deriv = 1,
  is_scale = 1,
  is_vary = 1,
  is_corr = 1
)
## model fitting
dgp_se_fit <- stan(
  file = 'DGPLVM_se.stan',
  data   = dgp_data,
  iter   = 2000,
  warmup = 1000,
  chains = 2,
  cores = 2,
  init = 0,
  control = list(adapt_delta=0.9)
)
dgp_m32_fit <- stan(
  file = 'DGPLVM_m32.stan',
  data   = dgp_data,
  iter   = 2000,
  warmup = 1000,
  chains = 2,
  cores = 2,
  init = 0,
  control = list(adapt_delta=0.9)
)
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
  ylim(c(-0.13, 0.13)) +
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
  scale_color_manual(values = c("#0072B2", "#E69F00" )) + ggtitle('(b)')
case_study_plot <- (dgp_se_plot + dgp_m32_plot + plot_layout(axes = 'collect', guides = 'collect') &
                      theme(legend.position = 'bottom'))
ggsave('dgp_case_study_plot.pdf',
       case_study_plot,
       dpi = 300,
       width = 50,
       height = 25,
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

## SE plot
dgp_se_rho_plot <- ggplot(dgp_se_summary_rho, aes(x = output_dim, y = mean)) +   
  theme_bw(base_size = 65, base_family = 'Times') +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = q5, ymax = q95), linewidth = 1) + 
  labs(x = 'Output dimensions', y = 'Posterior mean') +
  theme(axis.title = element_blank()) +
  facet_wrap(~ftitle) +
  ggtitle('(a)')
dgp_se_alpha_plot <- ggplot(dgp_se_summary_alpha, aes(x = output_dim, y = mean, colour = var_name)) +   
  theme_bw(base_size = 65, base_family = 'Times') +
  geom_point(size = 3, position = position_dodge(0.7)) +
  geom_errorbar(aes(ymin = q5, ymax = q95), linewidth = 1, position = position_dodge(0.7)) + 
  labs(x = 'Output dimensions', y = 'Posterior mean') + 
  facet_wrap(~ftitle) +
  theme(legend.title = element_blank(), axis.title = element_blank()) +
  scale_colour_manual(values = c("#0072B2", "#E69F00"))
dgp_se_sigma_plot <- ggplot(dgp_se_summary_sigma, aes(x = output_dim, y = mean, colour = var_name)) +   
  theme_bw(base_size = 65, base_family = 'Times') +
  geom_point(size = 3, position = position_dodge(0.7)) +
  geom_errorbar(aes(ymin = q5, ymax = q95), linewidth = 1, position = position_dodge(0.7)) + 
  labs(x = 'Output dimensions', y = 'Posterior mean') + 
  facet_wrap(~ftitle) +
  theme(legend.title = element_blank(), axis.title = element_blank()) +
  scale_colour_manual(values = c("#0072B2", "#E69F00"))

## Matern 3/2 plot
dgp_m32_rho_plot <- ggplot(dgp_m32_summary_rho, aes(x = output_dim, y = mean)) +   
  theme_bw(base_size = 65, base_family = 'Times') +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = q5, ymax = q95), linewidth = 1) + 
  labs(x = 'Output dimensions', y = 'Posterior mean') +
  facet_wrap(~ftitle) +
  theme(axis.title = element_blank()) +
  ggtitle('(b)')
dgp_m32_alpha_plot <- ggplot(dgp_m32_summary_alpha, aes(x = output_dim, y = mean, colour = var_name)) +   
  theme_bw(base_size = 65, base_family = 'Times') +
  geom_point(size = 3, position = position_dodge(0.7)) +
  geom_errorbar(aes(ymin = q5, ymax = q95), linewidth = 1, position = position_dodge(0.7)) + 
  labs(x = 'Output dimensions', y = 'Posterior mean') + 
  facet_wrap(~ftitle) +
  theme(legend.title = element_blank(), axis.title = element_blank()) +
  scale_colour_manual(values = c("#0072B2", "#E69F00"))
dgp_m32_sigma_plot <- ggplot(dgp_m32_summary_sigma, aes(x = output_dim, y = mean, colour = var_name)) +   
  theme_bw(base_size = 65, base_family = 'Times') +
  geom_point(size = 3, position = position_dodge(0.7)) +
  geom_errorbar(aes(ymin = q5, ymax = q95), linewidth = 1, position = position_dodge(0.7)) + 
  labs(x = 'Output dimensions', y = 'Posterior mean') + 
  facet_wrap(~ftitle) +
  theme(legend.title = element_blank(), axis.title = element_blank()) +
  scale_colour_manual(values = c("#0072B2", "#E69F00"))
## Combine
dgp_cellcycle_pars <- (dgp_se_rho_plot + dgp_se_alpha_plot + dgp_se_sigma_plot) /
                      (dgp_m32_rho_plot + dgp_m32_alpha_plot + dgp_m32_sigma_plot) + 
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