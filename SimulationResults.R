# Here we generate the summary results from the simulation studies
# A multilevel ANOVA model is used to analyse recovery of latent x through brms (doi.org/10.18637/jss.v080.i01)

# Libraries
library(dplyr)
library(tidyverse)
library(bayesplot)
library(ggplot2)
library(patchwork)
library(brms)
library(lemon)
library(grid)
library(gtable)
## Set variable names, sample points and number of trials
nobs <- 20
x_names <- sprintf('x[%s]', 1:nobs)
ntrials <- 50
#setwd("") # Set directory
## Plots for latent x
## GP sim data
GPtrialnames <- sprintf('GP_summary_trial[%s] .RDS', 1:ntrials)
GPtriallist <- list()
for(i in 1:ntrials){
  GPtriallist[[i]] <- readRDS(GPtrialnames[i])
}
GPx_out <- list()
for(i in 1:ntrials){
  GPx_out[[i]] <- subset(GPtriallist[[i]], var_names %in% x_names)
}
# Process data frame for fitting summary model
GP_x_comb <- bind_rows(GPx_out, .id = 'trial')
GP_x_comb$model_id <- paste(GP_x_comb$model_dim, GP_x_comb$is_deriv,
                            GP_x_comb$is_scale, GP_x_comb$is_vary,
                            GP_x_comb$is_corr)
GP_x_comb$is_vary <- as.factor(GP_x_comb$is_vary)
GP_x_comb$is_corr <- as.factor(GP_x_comb$is_corr)
GP_x_comb$deriv_scale <- as.factor(paste(GP_x_comb$is_deriv, GP_x_comb$is_scale))
GP_x_comb$model_dim <- as.factor(GP_x_comb$model_dim)
levels(GP_x_comb$model_dim) <- c(10, 20, 5)
GP_x_comb$model_dim <- factor(GP_x_comb$model_dim, levels = c(5, 10, 20))
levels(GP_x_comb$deriv_scale) <- c('Observed only','Unscaled derivative','Scaled derivative')
levels(GP_x_comb$is_vary) <- c('No','Yes')
levels(GP_x_comb$is_corr) <- c('No','Yes')
GP_x_comb$model_id <- as.factor(GP_x_comb$model_id)
GP_x_comb$trial <- as.factor(GP_x_comb$trial)
# Fit Multilevel ANOVA for GP data
m_gp <- brm(rmse~(1+deriv_scale + is_vary + is_corr + deriv_scale*is_vary + 
                    deriv_scale*is_corr)*model_dim + (1+deriv_scale+is_vary+is_corr|trial) +
              (1|model_id),
            data = GP_x_comb, chains = 2, cores = 2, file_refit = 'on_change')
## Periodic sim data
Pertrialnames <- sprintf('Periodic_summary_trial[%s] .RDS', 1:ntrials)
Pertriallist <- list()
for(i in 1:ntrials){
  Pertriallist[[i]] <- readRDS(Pertrialnames[i])
}
Perx_out <- list()
for(i in 1:ntrials){
  Perx_out[[i]] <- subset(Pertriallist[[i]], var_names %in% x_names)
}
# Process data frame for fitting summary model
Per_x_comb <- bind_rows(Perx_out, .id = 'trial')
Per_x_comb$model_id <- paste(Per_x_comb$model_dim, Per_x_comb$is_deriv,
                             Per_x_comb$is_scale, Per_x_comb$is_vary,
                             Per_x_comb$is_corr)
Per_x_comb$is_vary <- as.factor(Per_x_comb$is_vary)
Per_x_comb$is_corr <- as.factor(Per_x_comb$is_corr)
Per_x_comb$deriv_scale <- as.factor(paste(Per_x_comb$is_deriv, Per_x_comb$is_scale))
Per_x_comb$model_dim <- ordered(Per_x_comb$model_dim)
levels(Per_x_comb$model_dim) <- c(10, 20, 5)
Per_x_comb$model_dim <- factor(Per_x_comb$model_dim, levels = c(5, 10, 20))
levels(Per_x_comb$deriv_scale) <- c('Observed only','Unscaled derivative','Scaled derivative')
levels(Per_x_comb$is_vary) <- c('No','Yes')
levels(Per_x_comb$is_corr) <- c('No','Yes')
Per_x_comb$model_id <- as.factor(Per_x_comb$model_id)
Per_x_comb$trial <- as.factor(Per_x_comb$trial)
# Fit Multilevel ANOVA for Periodic data
m_per <- brm(rmse~(1 + deriv_scale + is_vary + is_corr + deriv_scale*is_vary + 
                     deriv_scale*is_corr)*model_dim + 
               (1+deriv_scale+is_vary+is_corr|trial) + (1|model_id),
             data = Per_x_comb, chains = 2, cores = 2, 
             file_refit = 'on_change')
# Generate results
# Effect of varying hyperparameters for GP data scenario
m_gp_vary <- conditional_effects(m_gp, effects = 'model_dim:is_vary', 
                                 conditions = make_conditions(m_gp, 'deriv_scale'),
                                 resolution = 300)
# Effect of correlated outputs for GP data scenario
m_gp_corr <- conditional_effects(m_gp, effects = 'model_dim:is_corr', 
                                 conditions = make_conditions(m_gp, 'deriv_scale'),
                                 resolution = 300)
# Effect of varying hyperparameters for Periodic data scenario
m_per_vary <- conditional_effects(m_per, effects = 'model_dim:is_vary', 
                                  conditions = make_conditions(m_per, 'deriv_scale'),
                                  resolution = 300)
# Effect of correlated outputs for Periodic data scenario
m_per_corr <- conditional_effects(m_per, effects = 'model_dim:is_corr', 
                                  conditions = make_conditions(m_per, 'deriv_scale'),
                                  resolution = 300)
# Main effects of derivative information and scaled derivatives for GP data scenario
m_gp_main <- conditional_effects(m_gp, effects = 'model_dim', 
                                 conditions = make_conditions(m_gp, 'deriv_scale'),
                                 resolution = 300)
# Main effects of derivative information and scaled derivatives for Periodic data scenario
m_per_main<- conditional_effects(m_per, effects = 'model_dim', 
                                 conditions = make_conditions(m_per, 'deriv_scale'),
                                 resolution = 300)
## Setup control (prior RMSE)
rmse_draws <- function(theta_hat, theta) {
  sqrt(mean((theta_hat - theta)^2))
}

# compare with naive estimate on the basis of the x_obs prior
n_obs <- 100000
rmse_x_naive <- vector(length = n_obs)
s_x <- 0.3
for (i in 1:n_obs) {
  # as if we only used x_obs to infer x
  x_obs <- rnorm(1, 0, s_x)
  draws_x_prior <- rnorm(2000, x_obs, s_x)
  rmse_x_naive[i] <- rmse_draws(draws_x_prior, 0)
}
hist(rmse_x_naive)
mean(rmse_x_naive)
sd(rmse_x_naive) / sqrt(n_obs)
## Plot for effect of varying hyperparameters for GP data scenario 
cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
df_gp_vary <- as.data.frame(m_gp_vary$`model_dim:is_vary`)
df_gp_vary$cond__ <- ordered(df_gp_vary$cond__)
levels(df_gp_vary$cond__) <- c('Derivative:No | Scaling:No',
                               'Derivative:Yes | Scaling:No',
                               'Derivative:Yes | Scaling:Yes')
df_gp_vary$cond__ <- factor(df_gp_vary$cond__, levels = c('Derivative:No | Scaling:No',
                                                          'Derivative:Yes | Scaling:Yes',
                                                          'Derivative:Yes | Scaling:No'))

p_gp_vary <- ggplot(df_gp_vary, aes(x = effect1__, y = estimate__, color = effect2__)) +
  theme_bw(base_size=65,
           base_family = 'Times') +
  geom_point(size = 5.5 ,
             position = position_dodge(width = 0.7)) +
  geom_errorbar(aes(ymin = lower__, ymax = upper__),
                width = 0.5,
                linewidth = 1.5,
                position = position_dodge(width = 0.7)) +
  geom_point(aes(x=0.5, y=0.407), colour = '#D55E00', size = 5.5) +
  annotate('text', x = 0.7, y = 0.39, label = 'Prior', size = 10, colour = '#D55E00') +
  facet_wrap(~cond__) +
  labs(x = 'Number of output dimensions', y = 'RMSE', colour = 'Varying hyperparameters') +
  guides(fill = 'none') + 
  theme(axis.ticks = element_line(linewidth = 3)) +
  scale_colour_manual(values = c("#E69F00", "#0072B2")) + ggtitle('(a)')

## Plot for effect of varying hyperparameters for Periodic data scenario 
df_per_vary <- as.data.frame(m_per_vary$`model_dim:is_vary`)
df_per_vary$cond__ <- ordered(df_per_vary$cond__)
levels(df_per_vary$cond__) <- c('Derivative:No | Scaling:No',
                                'Derivative:Yes | Scaling:No',
                                'Derivative:Yes | Scaling:Yes')
df_per_vary$cond__ <- factor(df_per_vary$cond__, levels = c('Derivative:No | Scaling:No',
                                                            'Derivative:Yes | Scaling:Yes',
                                                            'Derivative:Yes | Scaling:No'))

p_per_vary <- ggplot(df_per_vary, aes(x = effect1__, y = estimate__, color = effect2__)) +
  theme_bw(base_size=65,
           base_family = 'Times') +
  geom_point(size = 5.5 ,
             position = position_dodge(width = 0.7)) +
  geom_errorbar(aes(ymin = lower__, ymax = upper__),
                width = 0.5,
                linewidth = 1.5,
                position = position_dodge(width = 0.7)) +
  geom_point(aes(x=0.5, y=0.407), colour = '#D55E00', size = 5.5) +
  annotate('text', x = 0.7, y = 0.39, label = 'Prior', size = 10, colour = '#D55E00') +
  facet_wrap(~cond__) +
  labs(x = 'Number of output dimensions', y = 'RMSE', colour = 'Varying hyperparameters') +
  guides(fill = 'none') + 
  theme(axis.ticks = element_line(linewidth = 3)) +
  scale_colour_manual(values = c("#E69F00", "#0072B2")) + ggtitle('(b)')
# Combine varying hyperparameter effects plot
sim_int_eff_full_vary <- p_gp_vary + p_per_vary + 
  plot_layout(axes = 'collect', guides = 'collect') & 
  theme(axis.text = element_text(size = 85), legend.position = 'bottom', 
        axis.title = element_text(size = 85), legend.text = element_text(size = 85),
        legend.title = element_text(size = 85))
ggsave('sim_vary_eff_full.pdf',
       sim_int_eff_full_vary,
       width = 160,
       height = 40,
       units = 'cm',
       dpi = 300,
       limitsize = FALSE)
## Plot for effect of correlated outputs for GP data scenario 
df_gp_corr <- as.data.frame(m_gp_corr$`model_dim:is_corr`)
df_gp_corr$cond__ <- ordered(df_gp_corr$cond__)
levels(df_gp_corr$cond__) <- c('Derivative:No | Scaling:No',
                               'Derivative:Yes | Scaling:No',
                               'Derivative:Yes | Scaling:Yes')
df_gp_corr$cond__ <- factor(df_gp_corr$cond__, levels = c('Derivative:No | Scaling:No',
                                                          'Derivative:Yes | Scaling:Yes',
                                                          'Derivative:Yes | Scaling:No'))

p_gp_corr <- ggplot(df_gp_corr, aes(x = effect1__, y = estimate__, color = effect2__)) +
  theme_bw(base_size=65,
           base_family = 'Times') +
  geom_point(size = 5.5 ,
             position = position_dodge(width = 0.7)) +
  geom_errorbar(aes(ymin = lower__, ymax = upper__),
                width = 0.5,
                linewidth = 1.5,
                position = position_dodge(width = 0.7)) +
  geom_point(aes(x=0.5, y=0.407), colour = '#D55E00', size = 5.5) +
  annotate('text', x = 0.7, y = 0.39, label = 'Prior', size = 10, colour = '#D55E00') +
  facet_wrap(~cond__) +
  labs(x = 'Number of output dimensions', y = 'RMSE', colour = 'Correlated outputs') +
  guides(fill = 'none') + 
  theme(axis.ticks = element_line(linewidth = 3)) +
  scale_colour_manual(values = c("#E69F00", "#0072B2")) + ggtitle('(a)')
## Plot for effect of correlated outputs for Periodic data scenario 
df_per_corr <- as.data.frame(m_per_corr$`model_dim:is_corr`)
df_per_corr$cond__ <- ordered(df_per_corr$cond__)
levels(df_per_corr$cond__) <- c('Derivative:No | Scaling:No',
                                'Derivative:Yes | Scaling:No',
                                'Derivative:Yes | Scaling:Yes')
df_per_corr$cond__ <- factor(df_per_corr$cond__, levels = c('Derivative:No | Scaling:No',
                                                            'Derivative:Yes | Scaling:Yes',
                                                            'Derivative:Yes | Scaling:No'))

p_per_corr <- ggplot(df_per_corr, aes(x = effect1__, y = estimate__, color = effect2__)) +
  theme_bw(base_size=65,
           base_family = 'Times') +
  geom_point(size = 5.5 ,
             position = position_dodge(width = 0.7)) +
  geom_errorbar(aes(ymin = lower__, ymax = upper__),
                width = 0.5,
                linewidth = 1.5,
                position = position_dodge(width = 0.7)) +
  geom_point(aes(x=0.5, y=0.407), colour = '#D55E00', size = 5.5) +
  annotate('text', x = 0.7, y = 0.39, label = 'Prior', size = 10, colour = '#D55E00') +
  facet_wrap(~cond__) +
  labs(x = 'Number of output dimensions', y = 'RMSE', colour = 'Correlated outputs') +
  guides(fill = 'none') + 
  theme(axis.ticks = element_line(linewidth = 3)) +
  scale_colour_manual(values = c("#E69F00", "#0072B2")) + ggtitle('(b)')
# Combine effects of correlated outputs
sim_int_eff_full_corr <- p_gp_corr + p_per_corr + 
  plot_layout(axes = 'collect', guides = 'collect') & 
  theme(legend.title = element_text(size = 85), axis.title = element_text(size=85), 
        axis.text = element_text(size = 85), legend.text = element_text(size = 85),
        legend.position = 'bottom')
ggsave('sim_corr_eff_full.pdf',
       sim_int_eff_full_corr,
       width = 160,
       height = 40,
       units = 'cm',
       dpi = 300,
       limitsize = FALSE)
## Plot for main effects for GP data scenario 
df_gp_main <- as.data.frame(m_gp_main$`model_dim`)
df_gp_main$cond__ <- ordered(df_gp_main$cond__)
levels(df_gp_main$cond__) <- c('Derivative:No | Scaling:No',
                               'Derivative:Yes | Scaling:No',
                               'Derivative:Yes | Scaling:Yes')
df_gp_main$cond__ <- factor(df_gp_main$cond__, levels = c('Derivative:No | Scaling:No',
                                                          'Derivative:Yes | Scaling:Yes',
                                                          'Derivative:Yes | Scaling:No'))

p_gp_main <- ggplot(df_gp_main, aes(x = effect1__, y = estimate__)) +
  theme_bw(base_size=65,
           base_family = 'Times') +
  geom_point(size = 5.5,
             position = position_dodge(width = 0.7)) +
  geom_errorbar(aes(ymin = lower__, ymax = upper__),
                width = 0.5,
                linewidth = 1.5,
                position = position_dodge(width = 0.7)) +
  geom_point(aes(x=0.5, y=0.407), colour = '#D55E00', size = 5.5) +
  annotate('text', x = 0.7, y = 0.39, label = 'Prior', size = 10, colour = '#D55E00') +
  facet_wrap(~cond__) +
  labs(x = 'Number of output dimensions', y = 'RMSE') +
  guides(fill = 'none') + 
  theme(axis.ticks = element_line(linewidth = 3)) +
  scale_colour_manual(values = cbbPalette) + ggtitle('(a)')
## Plot for main effects for Periodic data scenario 
df_per_main <- as.data.frame(m_per_main$`model_dim`)
df_per_main$cond__ <- ordered(df_per_main$cond__)
levels(df_per_main$cond__) <- c('Derivative:No | Scaling:No',
                                'Derivative:Yes | Scaling:No',
                                'Derivative:Yes | Scaling:Yes')
df_per_main$cond__ <- factor(df_per_main$cond__, levels = c('Derivative:No | Scaling:No',
                                                            'Derivative:Yes | Scaling:Yes',
                                                            'Derivative:Yes | Scaling:No'))

p_per_main <- ggplot(df_per_main, aes(x = effect1__, y = estimate__)) +
  theme_bw(base_size = 65,
           base_family = 'Times') +
  geom_point(size = 5.5,
             position = position_dodge(width = 0.7)) +
  geom_errorbar(aes(ymin = lower__, ymax = upper__),
                width = 0.5,
                linewidth = 1.5,
                position = position_dodge(width = 0.7)) +
  geom_point(aes(x=0.5, y=0.407), colour = '#D55E00', size = 5.5) +
  annotate('text', x = 0.7, y = 0.39, label = 'Prior', size = 10, colour = '#D55E00') +
  facet_wrap(~cond__) +
  labs(x = 'Number of output dimensions', y = 'RMSE') +
  guides(fill = 'none') + 
  theme(axis.ticks = element_line(linewidth = 3)) +
  scale_colour_manual(values = cbbPalette) + ggtitle('(b)')
# Combine main effects plot
sim_main <- (p_gp_main + p_per_main + plot_layout(axes = 'collect')) & 
  theme(axis.text = element_text(size = 85), axis.title = element_text(size = 85))
ggsave('sim_main_full.pdf',
       sim_main,
       width = 150,
       height = 40,
       units = 'cm',
       dpi = 300,
       limitsize = FALSE)

## Plot for effect of varying hyperparameters for GP data scenario (removing unscaled derivative models)
df_gp_vary <- as.data.frame(m_gp_vary$`model_dim:is_vary`)
df_gp_vary$cond__ <- ordered(df_gp_vary$cond__)
levels(df_gp_vary$cond__) <- c('Derivative:No | Scaling:No',
                               'Derivative:Yes | Scaling:No',
                               'Derivative:Yes | Scaling:Yes')
df_gp_vary$cond__ <- factor(df_gp_vary$cond__, levels = c('Derivative:No | Scaling:No',
                                                          'Derivative:Yes | Scaling:Yes',
                                                          'Derivative:Yes | Scaling:No'))
df_gp_vary <- df_gp_vary[df_gp_vary$cond__!='Derivative:Yes | Scaling:No',]
p_gp_vary <- ggplot(df_gp_vary, aes(x = effect1__, y = estimate__, colour = effect2__)) +
  theme_bw(base_size=50,
           base_family = 'Times') +
  geom_point(size = 5.5,
             position = position_dodge(width = 0.7)) +
  geom_errorbar(aes(ymin = lower__, ymax = upper__),
                width = 0.5,
                linewidth = 1.5,
                position = position_dodge(width = 0.7)) +
  geom_point(aes(x=0.5, y=0.407), colour = '#D55E00', size = 5.5) +
  annotate('text', x = 0.7, y = 0.39, label = 'Prior', size = 10, colour = '#D55E00') +
  facet_wrap(~cond__) +
  labs(x = 'Number of output dimensions', y = 'RMSE', colour = 'Varying hyperparameters' ) +
  guides(fill = 'none') + 
  theme(axis.ticks = element_line(linewidth = 3)) +
  scale_colour_manual(values = c("#E69F00", "#0072B2")) + ggtitle('(a)')
## Plot for effect of varying hyperparameters for Periodic data scenario (removing unscaled derivative models)
df_per_vary <- as.data.frame(m_per_vary$`model_dim:is_vary`)
df_per_vary$cond__ <- ordered(df_per_vary$cond__)
levels(df_per_vary$cond__) <- c('Derivative:No | Scaling:No',
                                'Derivative:Yes | Scaling:No',
                                'Derivative:Yes | Scaling:Yes')
df_per_vary$cond__ <- factor(df_per_vary$cond__, levels = c('Derivative:No | Scaling:No',
                                                            'Derivative:Yes | Scaling:Yes',
                                                            'Derivative:Yes | Scaling:No'))
df_per_vary <- df_per_vary[df_per_vary$cond__!='Derivative:Yes | Scaling:No',]
p_per_vary <- ggplot(df_per_vary, aes(x = effect1__, y = estimate__, colour = effect2__)) +
  theme_bw(base_size=50,
           base_family = 'Times') +
  geom_point(size = 5.5,
             position = position_dodge(width = 0.7)) +
  geom_errorbar(aes(ymin = lower__, ymax = upper__),
                width = 0.5,
                linewidth = 1.5,
                position = position_dodge(width = 0.7)) +
  geom_point(aes(x=0.5, y=0.407), colour = '#D55E00', size = 5.5) +
  annotate('text', x = 0.7, y = 0.39, label = 'Prior', size = 10, colour = '#D55E00') +
  facet_wrap(~cond__) +
  labs(x = 'Number of output dimensions', y = 'RMSE', colour = 'Varying hyperparameters') +
  guides(fill = 'none') + 
  theme(axis.ticks = element_line(linewidth = 3)) +
  scale_colour_manual(values = c("#E69F00", "#0072B2")) + ggtitle('(b)')
# Combine varying hyperparameter effects plot (removing unscaled derivative models)
sim_int_eff_vary <- p_gp_vary + p_per_vary + 
  plot_layout(axes = 'collect', guides = 'collect') & 
  theme(legend.title = element_text(size=60), axis.title = element_text(size=60), 
        axis.text = element_text(size = 60), legend.text = element_text(size=60),
        legend.position = 'bottom')
ggsave('sim_vary_eff.pdf',
       sim_int_eff_vary,
       width = 90,
       height = 30,
       units = 'cm',
       dpi = 300,
       limitsize = FALSE)
## Plot for effect of correlated outputs for GP data scenario (removing unscaled derivative models)
df_gp_corr <- as.data.frame(m_gp_corr$`model_dim:is_corr`)
df_gp_corr$cond__ <- ordered(df_gp_corr$cond__)
levels(df_gp_corr$cond__) <- c('Derivative:No | Scaling:No',
                               'Derivative:Yes | Scaling:No',
                               'Derivative:Yes | Scaling:Yes')
df_gp_corr$cond__ <- factor(df_gp_corr$cond__, levels = c('Derivative:No | Scaling:No',
                                                          'Derivative:Yes | Scaling:Yes',
                                                          'Derivative:Yes | Scaling:No'))
df_gp_corr <- df_gp_corr[df_gp_corr$cond__!='Derivative:Yes | Scaling:No',]
p_gp_corr <- ggplot(df_gp_corr, aes(x = effect1__, y = estimate__, colour = effect2__)) +
  theme_bw(base_size=50,
           base_family = 'Times') +
  geom_point(size = 5.5,
             position = position_dodge(width = 0.7)) +
  geom_errorbar(aes(ymin = lower__, ymax = upper__),
                width = 0.5,
                linewidth = 1.5,
                position = position_dodge(width = 0.7)) +
  geom_point(aes(x=0.5, y=0.407), colour = '#D55E00', size = 5.5) +
  annotate('text', x = 0.7, y = 0.39, label = 'Prior', size = 10, colour = '#D55E00') +
  facet_wrap(~cond__) +
  labs(x = 'Number of output dimensions', y = 'RMSE', colour = 'Correlated outputs') +
  guides(fill = 'none') + 
  theme(axis.ticks = element_line(linewidth = 3)) +
  scale_colour_manual(values = c("#E69F00", "#0072B2")) + ggtitle('(a)')

## Plot for effect of correlated outputs for Periodic data scenario (removing unscaled derivative models)
df_per_corr <- as.data.frame(m_per_corr$`model_dim:is_corr`)
df_per_corr$cond__ <- ordered(df_per_corr$cond__)
levels(df_per_corr$cond__) <- c('Derivative:No | Scaling:No',
                                'Derivative:Yes | Scaling:No',
                                'Derivative:Yes | Scaling:Yes')
df_per_corr$cond__ <- factor(df_per_corr$cond__, levels = c('Derivative:No | Scaling:No',
                                                            'Derivative:Yes | Scaling:Yes',
                                                            'Derivative:Yes | Scaling:No'))
df_per_corr <- df_per_corr[df_per_corr$cond__!='Derivative:Yes | Scaling:No',]
p_per_corr <- ggplot(df_per_corr, aes(x = effect1__, y = estimate__, colour = effect2__)) +
  theme_bw(base_size=50,
           base_family = 'Times') +
  geom_point(size = 5.5,
             position = position_dodge(width = 0.7)) +
  geom_errorbar(aes(ymin = lower__, ymax = upper__),
                width = 0.5,
                linewidth = 1.5,
                position = position_dodge(width = 0.7)) +
  geom_point(aes(x=0.5, y=0.407), colour = '#D55E00', size = 5.5) +
  annotate('text', x = 0.7, y = 0.39, label = 'Prior', size = 10, colour = '#D55E00') +
  facet_wrap(~cond__) +
  labs(x = 'Number of output dimensions', y = 'RMSE', colour = 'Correlated outputs') +
  guides(fill = 'none') + 
  theme(axis.ticks = element_line(linewidth = 3)) +
  scale_colour_manual(values = c("#E69F00", "#0072B2")) + ggtitle('(b)')
# Combine correlated output effects plot (removing unscaled derivative models)
sim_int_eff_corr <- p_gp_corr + p_per_corr + 
  plot_layout(axes = 'collect', guides = 'collect') & 
  theme(legend.title = element_text(size=60), axis.title = element_text(size=60), 
        axis.text = element_text(size = 60), legend.text = element_text(size=60),
        legend.position = 'bottom')
ggsave('sim_corr_eff.pdf',
       sim_int_eff_corr,
       width = 90,
       height = 30,
       units = 'cm',
       dpi = 300,
       limitsize = FALSE)

## Plot for main effects for GP data scenario (removing unscaled derivative models)
df_gp_main <- as.data.frame(m_gp_main$`model_dim`)
df_gp_main$cond__ <- ordered(df_gp_main$cond__)
levels(df_gp_main$cond__) <- c('Derivative:No | Scaling:No',
                               'Derivative:Yes | Scaling:No',
                               'Derivative:Yes | Scaling:Yes')
df_gp_main$cond__ <- factor(df_gp_main$cond__, levels = c('Derivative:No | Scaling:No',
                                                          'Derivative:Yes | Scaling:Yes',
                                                          'Derivative:Yes | Scaling:No'))
df_gp_main <- df_gp_main[df_gp_main$cond__!='Derivative:Yes | Scaling:No',]
p_gp_main <- ggplot(df_gp_main, aes(x = effect1__, y = estimate__)) +
  theme_bw(base_size=50,
           base_family = 'Times') +
  geom_point(size = 5.5,
             position = position_dodge(width = 0.5)) +
  geom_errorbar(aes(ymin = lower__, ymax = upper__),
                width = 0.5,
                linewidth = 1.5,
                position = position_dodge(width = 0.5)) +
  geom_point(aes(x=0.5, y=0.407), colour = '#D55E00', size = 5.5) +
  annotate('text', x = 0.7, y = 0.39, label = 'Prior', size = 10, colour = '#D55E00') +
  facet_wrap(~cond__) +
  labs(x = 'Number of output dimensions', y = 'RMSE') +
  guides(fill = 'none') + 
  theme(axis.ticks = element_line(linewidth = 3)) +
  ggtitle('(a)')
## Plot for main effects for Periodic data scenario (removing unscaled derivative models)
df_per_main <- as.data.frame(m_per_main$`model_dim`)
df_per_main$cond__ <- ordered(df_per_main$cond__)
levels(df_per_main$cond__) <- c('Derivative:No | Scaling:No',
                                'Derivative:Yes | Scaling:No',
                                'Derivative:Yes | Scaling:Yes')
df_per_main$cond__ <- factor(df_per_main$cond__, levels = c('Derivative:No | Scaling:No',
                                                            'Derivative:Yes | Scaling:Yes',
                                                            'Derivative:Yes | Scaling:No'))
df_per_main <- df_per_main[df_per_main$cond__!='Derivative:Yes | Scaling:No',]
p_per_main <- ggplot(df_per_main, aes(x = effect1__, y = estimate__)) +
  theme_bw(base_size=50,
           base_family = 'Times') +
  geom_point(size = 5.5,
             position = position_dodge(width = 0.7)) +
  geom_errorbar(aes(ymin = lower__, ymax = upper__),
                width = 0.5,
                linewidth = 1.5,
                position = position_dodge(width = 0.7)) +
  geom_point(aes(x=0.5, y=0.407), colour = '#D55E00', size = 5.5) +
  annotate('text', x = 0.7, y = 0.39, label = 'Prior', size = 10, colour = '#D55E00') +
  facet_wrap(~cond__) +
  labs(x = 'Number of output dimensions', y = 'RMSE') +
  guides(fill = 'none') + 
  theme(axis.ticks = element_line(linewidth = 3)) +
  ggtitle('(b)')
# Combine main effects plot (removing unscaled derivative models)
sim_main <- p_gp_main + p_per_main + plot_layout(axes = 'collect') & 
  theme(axis.text = element_text(size = 60), axis.title = element_text(size = 60))
ggsave('sim_main.pdf',
       sim_main,
       width = 90,
       height = 30,
       units = 'cm',
       dpi = 300,
       limitsize = FALSE)

## Plots for hyperparameters
GPtrialnames <- sprintf('GP_summary_trial[%s] .RDS', 1:ntrials)
GPtriallist <- list()
for(i in 1:ntrials){
  GPtriallist[[i]] <- readRDS(GPtrialnames[i])
}
'%!in%'<- function(a,b) ! a %in% b
GPpars_out <- list()
for(i in 1:ntrials){
  GPpars_out[[i]] <- subset(GPtriallist[[i]], var_names %!in% x_names)
}
# Prepare data for plots
GP_pars_comb <- bind_rows(GPpars_out, .id = 'trial')
GP_pars_comb$model_id <- paste(GP_pars_comb$model_dim, GP_pars_comb$is_deriv,
                               GP_pars_comb$is_scale, GP_pars_comb$is_vary,
                               GP_pars_comb$is_corr)
GP_pars_comb$is_vary <- as.factor(GP_pars_comb$is_vary)
GP_pars_comb$is_corr <- as.factor(GP_pars_comb$is_corr)
GP_pars_comb$deriv_scale <- as.factor(paste(GP_pars_comb$is_deriv, GP_pars_comb$is_scale))
GP_pars_comb$model_dim <- as.factor(GP_pars_comb$model_dim)
levels(GP_pars_comb$model_dim) <- c(10, 20, 5)
GP_pars_comb$model_dim <- factor(GP_pars_comb$model_dim, levels = c(5, 10, 20))
levels(GP_pars_comb$deriv_scale) <- c('Observed only','Unscaled derivative','Scaled derivative')
levels(GP_pars_comb$is_vary) <- c('No','Yes')
levels(GP_pars_comb$is_corr) <- c('No','Yes')
dgp_id <- c('[5]D 1 1 1 1', '[10]D 1 1 1 1','[20]D 1 1 1 1')
dgp_sim_pars <- subset(GP_pars_comb, model_id %in% dgp_id)
df <- dgp_sim_pars %>% separate(var_names, c('var_names', 'output_type'))
dgp_sim_rho <- subset(df, var_names=='rho')
dgp_sim_alpha <- subset(df, var_names=='alpha')
dgp_sim_alpha$output_type <- as.factor(dgp_sim_alpha$output_type)
levels(dgp_sim_alpha$output_type) <- c("f '", "f")
dgp_sim_sigma <- subset(df, var_names=='sigma')
dgp_sim_sigma$output_type <- as.factor(dgp_sim_sigma$output_type)
levels(dgp_sim_sigma$output_type) <- c("f '", "f")
# Plot for GP length scale rho
dgp_sim_rho_plot <- ggplot(data = dgp_sim_rho, aes(x=model_dim, y=rmse)) +
  theme_bw(base_size = 35, base_family = 'Times') + 
  geom_violin(linewidth = 1) +
  facet_wrap(~var_names) +
  labs(x='Number of output dimensions', y='RMSE')
# Plot for GP marginal SD alphas
dgp_sim_alpha_plot <- ggplot(data = dgp_sim_alpha, aes(x=model_dim, y=rmse, colour=output_type)) +
  theme_bw(base_size = 35, base_family = 'Times') + 
  geom_violin(linewidth = 1) +
  facet_wrap(~var_names) +
  labs(x='Number of output dimensions', y='RMSE') + theme(axis.ticks = element_line(linewidth=3),
                                                          legend.title = element_blank()) +
  scale_colour_manual(values=c("#0072B2", "#E69F00"))
# Plot for error SD sigmas
dgp_sim_sigma_plot <- ggplot(data = dgp_sim_sigma, aes(x=model_dim, y=rmse, colour=output_type)) +
  theme_bw(base_size = 35, base_family = 'Times') + 
  geom_violin(linewidth = 1) +
  facet_wrap(~var_names) +
  labs(x='Number of output dimensions', y='RMSE') + theme(axis.ticks = element_line(linewidth=3), 
                                                          legend.title = element_blank()) +
  scale_colour_manual(values=c("#0072B2", "#E69F00"))
# Combine plots
dgp_sim_pars_plot <- dgp_sim_rho_plot + dgp_sim_alpha_plot + dgp_sim_sigma_plot +
  plot_layout(axis_titles = 'collect', guides = 'collect')
ggsave('dgp_sim_pars_plot.pdf',
       dgp_sim_pars_plot,
       dpi = 300,
       width = 50,
       height = 15,
       units = 'cm')

## Plots of hyperparameters without varying across output dimensions
dgp_id <- c('[5]D 1 1 0 1', '[10]D 1 1 0 1','[20]D 1 1 0 1')
dgp_sim_pars <- subset(GP_pars_comb, model_id %in% dgp_id)
df <- dgp_sim_pars %>% separate(var_names, c('var_names', 'output_type'))
dgp_sim_rho <- subset(df, var_names=='rho')
dgp_sim_alpha <- subset(df, var_names=='alpha')
dgp_sim_alpha$output_type <- as.factor(dgp_sim_alpha$output_type)
levels(dgp_sim_alpha$output_type) <- c("f '", "f")
dgp_sim_sigma <- subset(df, var_names=='sigma')
dgp_sim_sigma$output_type <- as.factor(dgp_sim_sigma$output_type)
levels(dgp_sim_sigma$output_type) <- c("f '", "f")
# Plot for GP lengthscale rho
dgp_sim_rho_plot <- ggplot(data = dgp_sim_rho, aes(x=model_dim, y=rmse)) +
  theme_bw(base_size = 35, base_family = 'Times') + 
  geom_violin(linewidth = 1) +
  facet_wrap(~var_names) +
  labs(x='Number of output dimensions', y='RMSE') + ggtitle('(a)')
# Plot for GP marginal SD alphas
dgp_sim_alpha_plot <- ggplot(data = dgp_sim_alpha, aes(x=model_dim, y=rmse, colour=output_type)) +
  theme_bw(base_size = 35, base_family = 'Times') + 
  geom_violin(linewidth = 1) +
  facet_wrap(~var_names) +
  labs(x='Number of output dimensions', y='RMSE') + theme(axis.ticks = element_line(linewidth=3),
                                                          legend.title = element_blank()) +
  scale_colour_manual(values=c("#0072B2", "#E69F00"))
# Plot for error SD sigmas
dgp_sim_sigma_plot <- ggplot(data = dgp_sim_sigma, aes(x=model_dim, y=rmse, colour=output_type)) +
  theme_bw(base_size = 35, base_family = 'Times') + 
  geom_violin(linewidth = 1) +
  facet_wrap(~var_names) +
  labs(x='Number of output dimensions', y='RMSE') + theme(axis.ticks = element_line(linewidth=3), 
                                                          legend.title = element_blank()) +
  scale_colour_manual(values=c("#0072B2", "#E69F00"))
# Combine plots
dgp_sim_novary_plot <- dgp_sim_rho_plot + dgp_sim_alpha_plot + dgp_sim_sigma_plot +
  plot_layout(axis_titles = 'collect', guides = 'collect') & theme(axis.title = element_blank(),
                                                                   legend.position = 'none')

## Plots of hyperparameters without correlated output modification
dgp_id <- c('[5]D 1 1 1 0', '[10]D 1 1 1 0','[20]D 1 1 1 0')
dgp_sim_pars <- subset(GP_pars_comb, model_id %in% dgp_id)
df <- dgp_sim_pars %>% separate(var_names, c('var_names', 'output_type'))
dgp_sim_rho <- subset(df, var_names=='rho')
dgp_sim_alpha <- subset(df, var_names=='alpha')
dgp_sim_alpha$output_type <- as.factor(dgp_sim_alpha$output_type)
levels(dgp_sim_alpha$output_type) <- c("f '", "f")
dgp_sim_sigma <- subset(df, var_names=='sigma')
dgp_sim_sigma$output_type <- as.factor(dgp_sim_sigma$output_type)
levels(dgp_sim_sigma$output_type) <- c("f '", "f")
# Plot for GP lengthscale rho
dgp_sim_rho_plot <- ggplot(data = dgp_sim_rho, aes(x=model_dim, y=rmse)) +
  theme_bw(base_size = 35, base_family = 'Times') + 
  geom_violin(linewidth = 1) +
  facet_wrap(~var_names) +
  labs(x='Number of output dimensions', y='RMSE') + ggtitle('(b)')
# Plot for GP marginal SD alphas
dgp_sim_alpha_plot <- ggplot(data = dgp_sim_alpha, aes(x=model_dim, y=rmse, colour=output_type)) +
  theme_bw(base_size = 35, base_family = 'Times') + 
  geom_violin(linewidth = 1) +
  facet_wrap(~var_names) +
  labs(x='Number of output dimensions', y='RMSE') + theme(axis.ticks = element_line(linewidth=3),
                                                          legend.title = element_blank()) +
  scale_colour_manual(values=c("#0072B2", "#E69F00"))
# Plot for error SD sigmas
dgp_sim_sigma_plot <- ggplot(data = dgp_sim_sigma, aes(x=model_dim, y=rmse, colour=output_type)) +
  theme_bw(base_size = 35, base_family = 'Times') + 
  geom_violin(linewidth = 1) +
  facet_wrap(~var_names) +
  labs(x='Number of output dimensions', y='RMSE') + theme(axis.ticks = element_line(linewidth=3), 
                                                          legend.title = element_blank()) +
  scale_colour_manual(values=c("#0072B2", "#E69F00"))
dgp_sim_nocorr_plot <- dgp_sim_rho_plot + dgp_sim_alpha_plot + dgp_sim_sigma_plot +
  plot_layout(axis_titles = 'collect', guides = 'collect') & theme(axis.title.x = element_blank())

## Plots of hyperparameters without varying across output dimensions as well as without correlated output modification
dgp_id <- c('[5]D 1 1 0 0', '[10]D 1 1 0 0','[20]D 1 1 0 0')
dgp_sim_pars <- subset(GP_pars_comb, model_id %in% dgp_id)
df <- dgp_sim_pars %>% separate(var_names, c('var_names', 'output_type'))
dgp_sim_rho <- subset(df, var_names=='rho')
dgp_sim_alpha <- subset(df, var_names=='alpha')
dgp_sim_alpha$output_type <- as.factor(dgp_sim_alpha$output_type)
levels(dgp_sim_alpha$output_type) <- c("f '", "f")
dgp_sim_sigma <- subset(df, var_names=='sigma')
dgp_sim_sigma$output_type <- as.factor(dgp_sim_sigma$output_type)
levels(dgp_sim_sigma$output_type) <- c("f '", "f")
# Plot for GP lengthscale rho
dgp_sim_rho_plot <- ggplot(data = dgp_sim_rho, aes(x=model_dim, y=rmse)) +
  theme_bw(base_size = 35, base_family = 'Times') + 
  geom_violin(linewidth = 1) +
  facet_wrap(~var_names) +
  labs(x='Number of output dimensions', y='RMSE') + ggtitle('(c)')
# Plot for GP marginal SD alphas
dgp_sim_alpha_plot <- ggplot(data = dgp_sim_alpha, aes(x=model_dim, y=rmse, colour=output_type)) +
  theme_bw(base_size = 35, base_family = 'Times') + 
  geom_violin(linewidth = 1) +
  facet_wrap(~var_names) +
  labs(x='Number of output dimensions', y='RMSE') + theme(axis.ticks = element_line(linewidth=3),
                                                          legend.title = element_blank()) +
  scale_colour_manual(values=c("#0072B2", "#E69F00"))
# Plot for error SD alphas
dgp_sim_sigma_plot <- ggplot(data = dgp_sim_sigma, aes(x=model_dim, y=rmse, colour=output_type)) +
  theme_bw(base_size = 35, base_family = 'Times') + 
  geom_violin(linewidth = 1) +
  facet_wrap(~var_names) +
  labs(x='Number of output dimensions', y='RMSE') + theme(axis.ticks = element_line(linewidth=3), 
                                                          legend.title = element_blank()) +
  scale_colour_manual(values=c("#0072B2", "#E69F00"))
# Combine plots
dgp_sim_novarycorr_plot <- dgp_sim_rho_plot + dgp_sim_alpha_plot + dgp_sim_sigma_plot +
  plot_layout(axis_titles = 'collect', guides = 'collect') & theme(axis.title.y = element_blank(),
                                                                   legend.position = 'none')
# Combine final plot
dgp_sim_novarycorr_comb <- dgp_sim_novary_plot / dgp_sim_nocorr_plot / dgp_sim_novarycorr_plot +
  plot_layout(axis_titles = 'collect', guides = 'collect')
ggsave('dgp_sim_novarycorr_plot.pdf',
       dgp_sim_novarycorr_comb,
       dpi = 300,
       width = 50,
       height = 40,
       units = 'cm')

