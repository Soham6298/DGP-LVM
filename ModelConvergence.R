# We show and check MCMC convergence for DGP-LVM and other fitted models for the simulated data scenarios
# We check Rhats, Bulk-ESS and Tail-ESS for fitted models

# Libraries
library(dplyr)
library(tidyverse)
library(bayesplot)
library(ggplot2)
library(patchwork)
## Set variable names, sample points and number of trials
nobs <- 20
x_names <- sprintf('x[%s]', 1:nobs)
ntrials <- 50
#setwd("") # Set directory for data
## GP sim data
GPtrialnames <- sprintf('GP_summary_trial[%s] .RDS', 1:ntrials)
GPtriallist <- list()
for(i in 1:ntrials){
  GPtriallist[[i]] <- readRDS(GPtrialnames[i])
}
## Periodic data
Pertrialnames <- sprintf('Periodic_summary_trial[%s] .RDS', 1:ntrials)
Pertriallist <- list()
for(i in 1:ntrials){
  Pertriallist[[i]] <- readRDS(Pertrialnames[i])
}
# Prep data frame for convergence plots
GP_comb <- bind_rows(GPtriallist, .id = 'trial')
GP_comb$model_id <- paste(GP_comb$model_dim, GP_comb$is_deriv, GP_comb$is_scale,
                          GP_comb$is_vary, GP_comb$is_corr)
# Convergence for latent x
GP_comb_x <- subset(GP_comb, var_names %in% x_names)
GPsummarydata <- GP_comb_x %>%
  group_by(model_id) %>%
  summarise(mtrue_value = mean(true_value),
            mmean = mean(mean), msd = mean(sd), mabs_bias = mean(abs_bias),
            mrmse = mean(rmse), mrhat = mean(rhat), mess_bulk = mean(ess_bulk),
            mess_tail = mean(ess_tail), mruntime = mean(runtime), model_dim = min(model_dim),
            is_deriv = min(is_deriv), is_scale = min(is_scale), is_vary = min(is_vary),
            is_corr = min(is_corr))
GPsummarydata$model_dim <- as.factor(GPsummarydata$model_dim)
levels(GPsummarydata$model_dim) <- c(10, 20, 5)
GPsummarydata$model_dim <- factor(GPsummarydata$model_dim, levels = c(5, 10, 20))
GPsummarydata$rhat_name <- 'Rhat'
GPsummarydata$bess_name <- 'Bulk-ESS'
GPsummarydata$tess_name <- 'Tail-ESS'
# Rhat plot
x_rhat_summary_plot <- ggplot(GPsummarydata, aes(x = model_dim, y = mrhat)) + 
  theme_bw(base_size = 35, base_family = 'Times') + 
  geom_violin(linewidth = 1) + 
  geom_jitter(width = 0.1) + 
  facet_wrap(~rhat_name) +
  labs(x = 'Output dimensions', y = 'Values') + 
  theme(axis.ticks = element_line(linewidth = 3), legend.position = 'none') +
  ggtitle('(a)')
# Bulk-ESS plot
x_ess_bulk_summary_plot <- ggplot(GPsummarydata, aes(x = model_dim, y = mess_bulk)) + 
  theme_bw(base_size = 35, base_family = 'Times') +
  geom_violin(linewidth = 1) + 
  geom_jitter(width = 0.1) + 
  facet_wrap(~bess_name) +
  labs(x = 'Output dimensions', y = 'Values') + 
  theme(axis.ticks = element_line(linewidth = 3), legend.position = 'none') +
# Tail-ESS plot
x_ess_tail_summary_plot <- ggplot(GPsummarydata, aes(x = model_dim, y = mess_tail)) + 
  theme_bw(base_size = 35, base_family = 'Times') +
  geom_violin(linewidth = 1) + 
  geom_jitter(width = 0.1) + 
  facet_wrap(~tess_name) +
  labs(x = 'Output dimensions', y = 'Values') + 
  theme(axis.ticks = element_line(linewidth = 3),legend.position = 'none') +
# Convergence for hyperparameters
'%!in%'<- function(a,b) ! a %in% b
GP_comb_pars <- subset(GP_comb, var_names %!in% x_names)
GPsummarydata <- GP_comb_pars %>%
  group_by(model_id) %>%
  summarise(mtrue_value = mean(true_value),
            mmean = mean(mean), msd = mean(sd), mabs_bias = mean(abs_bias),
            mrmse = mean(rmse), mrhat = mean(rhat), mess_bulk = mean(ess_bulk),
            mess_tail = mean(ess_tail), mruntime = mean(runtime), model_dim = min(model_dim),
            is_deriv = min(is_deriv), is_scale = min(is_scale), is_vary = min(is_vary),
            is_corr = min(is_corr))
GPsummarydata$model_dim <- as.factor(GPsummarydata$model_dim)
levels(GPsummarydata$model_dim) <- c(10, 20, 5)
GPsummarydata$model_dim <- factor(GPsummarydata$model_dim, levels = c(5, 10, 20))
GPsummarydata$rhat_name <- 'Rhat'
GPsummarydata$bess_name <- 'Bulk-ESS'
GPsummarydata$tess_name <- 'Tail-ESS'
# Rhat plot
pars_rhat_summary_plot <- ggplot(GPsummarydata, aes(x = model_dim, y = mrhat)) + 
  geom_violin(linewidth = 1) + 
  geom_jitter(width = 0.1) + 
  facet_wrap(~rhat_name) +
  theme_bw(base_size = 35, base_family = 'Times') + 
  labs(x = 'Output dimensions', y = 'Values') + 
  theme(axis.ticks = element_line(linewidth = 3), legend.position = 'none') +
  ggtitle('(b)')
# Bulk-ESS plot
pars_ess_bulk_summary_plot <- ggplot(GPsummarydata, aes(x = model_dim, y = mess_bulk)) + 
  geom_violin(linewidth = 1) + 
  geom_jitter(width = 0.1) + 
  facet_wrap(~bess_name) +
  theme_bw(base_size = 35, base_family = 'Times') + 
  labs(x = 'Output dimensions', y = 'Values') + 
  theme(axis.ticks = element_line(linewidth = 3), legend.position = 'none') +
# Tail-ESS plot
pars_ess_tail_summary_plot <- ggplot(GPsummarydata, aes(x = model_dim, y = mess_tail)) + 
  geom_violin(linewidth = 1) + 
  geom_jitter(width = 0.1) + 
  facet_wrap(~tess_name) +
  theme_bw(base_size = 35, base_family = 'Times') + 
  labs(x = 'Output dimensions', y = 'Values') + 
  theme(axis.ticks = element_line(linewidth = 3),legend.position = 'none') +
# Combine the plots
plot_valid_gp <- ((x_rhat_summary_plot + x_ess_bulk_summary_plot + x_ess_tail_summary_plot) /
                 (pars_rhat_summary_plot + pars_ess_bulk_summary_plot + pars_ess_tail_summary_plot) + 
                   plot_layout(axis_titles = 'collect') & theme(axis.title = element_blank()))
gt <- patchwork::patchworkGrob(plot_valid_gp)
g <- arrangeGrob(gt, bottom = textGrob("Number of output dimensions", vjust = -0.5, 
                                       gp = gpar(fontsize=35, fontfamily='Times')))
ggsave('gp_valid.pdf',
       g,
       dpi = 300,
       width = 50,
       height = 30,
       units = 'cm')
## Periodic data
Per_comb <- bind_rows(Pertriallist, .id = 'trial')
Per_comb$model_id <- paste(Per_comb$model_dim, Per_comb$is_deriv, Per_comb$is_scale,
                           Per_comb$is_vary, Per_comb$is_corr)
# Convergence for latent x
Per_comb_x <- subset(Per_comb, var_names %in% x_names)
Persummarydata <- Per_comb_x %>%
  group_by(model_id) %>%
  summarise(mtrue_value = mean(true_value),
            mmean = mean(mean), msd = mean(sd), mabs_bias = mean(abs_bias),
            mrmse = mean(rmse), mrhat = mean(rhat), mess_bulk = mean(ess_bulk),
            mess_tail = mean(ess_tail), mruntime = mean(runtime), model_dim = min(model_dim),
            is_deriv = min(is_deriv), is_scale = min(is_scale), is_vary = min(is_vary),
            is_corr = min(is_corr))
Persummarydata$model_dim <- as.factor(Persummarydata$model_dim)
levels(Persummarydata$model_dim) <- c(10, 20, 5)
Persummarydata$model_dim <- factor(Persummarydata$model_dim, levels = c(5, 10, 20))
Persummarydata$rhat_name <- 'Rhat'
Persummarydata$bess_name <- 'Bulk-ESS'
Persummarydata$tess_name <- 'Tail-ESS'
# Rhat plot
x_rhat_summary_plot <- ggplot(Persummarydata, aes(x = model_dim, y = mrhat)) + 
  geom_violin(linewidth = 1) + 
  geom_jitter(width = 0.1) + 
  facet_wrap(~rhat_name) +
  theme_bw(base_size = 35, base_family = 'Times') + 
  labs(x = 'Output dimensions', y = 'Values') + 
  theme(axis.ticks = element_line(linewidth = 3), legend.position = 'none') +
  ggtitle('(a)')
# Bulk-ESS plot
x_ess_bulk_summary_plot <- ggplot(Persummarydata, aes(x = model_dim, y = mess_bulk)) + 
  geom_violin(linewidth = 1) + 
  geom_jitter(width = 0.1) + 
  facet_wrap(~bess_name) +
  theme_bw(base_size = 35, base_family = 'Times') + 
  labs(x = 'Output dimensions', y = 'Values') + 
  theme(axis.ticks = element_line(linewidth = 3), legend.position = 'none') +
# Tail-ESS plot
x_ess_tail_summary_plot <- ggplot(Persummarydata, aes(x = model_dim, y = mess_tail)) + 
  geom_violin(linewidth = 1) + 
  geom_jitter(width = 0.1) + 
  facet_wrap(~tess_name) +
  theme_bw(base_size = 35, base_family = 'Times') + 
  labs(x = 'Output dimensions', y = 'Values') + 
  theme(axis.ticks = element_line(linewidth = 3), legend.position = 'none') +
# Convergence for hyperparameters
Per_comb_pars <- subset(Per_comb, var_names %!in% x_names)
Persummarydata <- Per_comb_pars %>%
  group_by(model_id) %>%
  summarise(mtrue_value = mean(true_value),
            mmean = mean(mean), msd = mean(sd), mabs_bias = mean(abs_bias),
            mrmse = mean(rmse), mrhat = mean(rhat), mess_bulk = mean(ess_bulk),
            mess_tail = mean(ess_tail), mruntime = mean(runtime), model_dim = min(model_dim),
            is_deriv = min(is_deriv), is_scale = min(is_scale), is_vary = min(is_vary),
            is_corr = min(is_corr))
Persummarydata$model_dim <- as.factor(Persummarydata$model_dim)
levels(Persummarydata$model_dim) <- c(10, 20, 5)
Persummarydata$model_dim <- factor(Persummarydata$model_dim, levels = c(5, 10, 20))
Persummarydata$rhat_name <- 'Rhat'
Persummarydata$bess_name <- 'Bulk-ESS'
Persummarydata$tess_name <- 'Tail-ESS'
# Rhat plot
pars_rhat_summary_plot <- ggplot(Persummarydata, aes(x = model_dim, y = mrhat)) + 
  geom_violin(linewidth = 1) + 
  geom_jitter(width = 0.1) + 
  facet_wrap(~rhat_name) +
  theme_bw(base_size = 35, base_family = 'Times') + 
  labs(x = 'Output dimensions', y = 'Values') + 
  theme(axis.ticks = element_line(linewidth = 3), legend.position = 'none') +
  ggtitle('(b)')
# Bulk-ESS plot
pars_ess_bulk_summary_plot <- ggplot(Persummarydata, aes(x = model_dim, y = mess_bulk)) + 
  geom_violin(linewidth = 1) + 
  geom_jitter(width = 0.1) + 
  facet_wrap(~bess_name) +
  theme_bw(base_size = 35, base_family = 'Times') + 
  labs(x = 'Output dimensions', y = 'Values') + 
  theme(axis.ticks = element_line(linewidth = 3), legend.position = 'none') +
# Tail-ESS plot
pars_ess_tail_summary_plot <- ggplot(Persummarydata, aes(x = model_dim, y = mess_tail)) + 
  geom_violin(linewidth = 1) + 
  geom_jitter(width = 0.1) + 
  facet_wrap(~tess_name) +
  theme_bw(base_size = 35, base_family = 'Times') + 
  labs(x = 'Output dimensions', y = 'Values') + 
  theme(axis.ticks = element_line(linewidth = 3), legend.position = 'none') +
# Combine plots
plot_valid_per <- ((x_rhat_summary_plot + x_ess_bulk_summary_plot + x_ess_tail_summary_plot) /
                      (pars_rhat_summary_plot + pars_ess_bulk_summary_plot + pars_ess_tail_summary_plot) + 
                      plot_layout(axis_titles = 'collect') & theme(axis.title = element_blank()))
gt <- patchwork::patchworkGrob(plot_valid_per)
g <- arrangeGrob(gt, bottom = textGrob("Number of output dimensions", vjust = -0.5, 
                                       gp = gpar(fontsize=35, fontfamily='Times')))
ggsave('per_valid.pdf',
       g,
       dpi = 300,
       width = 50,
       height = 30,
       units = 'cm')

