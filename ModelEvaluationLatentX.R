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
# Source functions
source('DgplvmSimFns.R')
## Set variable names, sample points and number of trials
simdata_out <- readRDS('simulation results/dgplvm_simout_se.rds') #change according to output file name from SimStudy
simdata_out <- subset(simdata_out, class == 'x')
simdata_out$vary <- as.factor(simdata_out$vary)
simdata_out$corr <- as.factor(simdata_out$corr)
simdata_out$deriv_scale <- as.factor(paste(simdata_out$deriv, simdata_out$scale))
simdata_out$d <- as.factor(simdata_out$d)
levels(simdata_out$deriv_scale) <- c('Observed only','Unscaled derivative','Scaled derivative')
levels(simdata_out$vary) <- c('No','Yes')
levels(simdata_out$corr) <- c('No','Yes')
simdata_out$model_id <- as.factor(simdata_out$model_id)
simdata_out$sim_id <- as.factor(simdata_out$sim_id)
# Fit Multilevel ANOVA for GP data
m_rmse <- brm(rmse~(1 + deriv_scale + vary + corr + deriv_scale * vary + 
                     deriv_scale * corr) * d + (1 + deriv_scale + vary + corr | sim_id) + (1 | model_id),
             data = simdata_out, chains = 2, cores = 2, file_refit = 'on_change')
m_mae <- brm(mae~(1 + deriv_scale + vary + corr + deriv_scale * vary + 
                    deriv_scale * corr) * d + (1 + deriv_scale + vary + corr | sim_id) + (1 | model_id),
             data = simdata_out, chains = 2, cores = 2, file_refit = 'on_change')
# Generate results
# Main effects of derivative information and scaled derivatives RMSE
m_rmse_main <- conditional_effects(m_rmse, effects = 'd', 
                                   conditions = make_conditions(m_rmse, 'deriv_scale'),
                                   resolution = 300)
# Effect of varying hyperparameters RMSE
m_rmse_vary <- conditional_effects(m_rmse, effects = 'd:vary', 
                                 conditions = make_conditions(m_rmse, 'deriv_scale'),
                                 resolution = 300)
# Effect of correlated outputs RMSE
m_rmse_corr <- conditional_effects(m_rmse, effects = 'd:corr', 
                                 conditions = make_conditions(m_rmse, 'deriv_scale'),
                                 resolution = 300)
# Main effects of derivative information and scaled derivatives MAE
m_mae_main<- conditional_effects(m_mae, effects = 'd', 
                                 conditions = make_conditions(m_mae, 'deriv_scale'),
                                 resolution = 300)
# Effect of varying hyperparameters MAE
m_mae_vary <- conditional_effects(m_mae, effects = 'd:vary', 
                                  conditions = make_conditions(m_mae, 'deriv_scale'),
                                  resolution = 300)
# Effect of correlated outputs for MAE
m_mae_corr <- conditional_effects(m_mae, effects = 'd:corr', 
                                  conditions = make_conditions(m_mae, 'deriv_scale'),
                                  resolution = 300)

## Setup control (prior RMSE)
rmse_draws <- function(theta_hat, theta) {
  sqrt(mean((theta_hat - theta)^2))
}

# compare with naive estimate on the basis of the x_obs prior
n_obs <- 100000
rmse_x_naive <- vector(length = n_obs)
mae_x_naive <- vector(length = n_obs)
s_x <- 0.3
for (i in 1:n_obs) {
  # as if we only used x_obs to infer x
  x_obs <- rnorm(1, 0, s_x)
  draws_x_prior <- rnorm(2000, x_obs, s_x)
  rmse_x_naive[i] <- rmse_draws(draws_x_prior, 0)
  mae_x_naive[i] <- mae_draws(draws_x_prior, 0)
}
naive_rmse <- mean(rmse_x_naive)
naive_mae <- mean(mae_x_naive)
cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

# Plot for main effects RMSE 
df_rmse_main <- as.data.frame(m_rmse_main$`d`)
df_rmse_main$cond__ <- ordered(df_rmse_main$cond__)
levels(df_rmse_main$cond__) <- c('No derivative',
                                 'Unscaled derivative',
                                 'Scaled derivative')
df_rmse_main$cond__ <- factor(df_rmse_main$cond__, levels = c('No derivative',
                                                              'Scaled derivative',
                                                              'Unscaled derivative'))

p_rmse_main <- ggplot(df_rmse_main, aes(x = effect1__, y = estimate__)) +
  theme_bw(base_size=65,
           base_family = 'Times') +
  geom_point(size = 5.5,
             position = position_dodge(width = 0.7)) +
  geom_errorbar(aes(ymin = lower__, ymax = upper__),
                width = 0.3,
                linewidth = 1,
                position = position_dodge(width = 0.7)) +
  geom_point(aes(x=0.5, y=naive_rmse), colour = '#D55E00', size = 5.5) +
  annotate('text', x = 0.7, y = naive_rmse + 0.01, label = 'Prior', size = 10, colour = '#D55E00') +
  facet_wrap(~cond__) +
  labs(x = 'Number of output dimensions', y = 'RMSE') +
  guides(fill = 'none') + 
  theme(axis.ticks = element_line(linewidth = 3)) +
  scale_colour_manual(values = cbbPalette) + ggtitle('(a)')

## Plot for effect of varying hyperparameters RMSE
df_rmse_vary <- as.data.frame(m_rmse_vary$`d:vary`)
df_rmse_vary$cond__ <- ordered(df_rmse_vary$cond__)
levels(df_rmse_vary$cond__) <- c('No derivative',
                                 'Unscaled derivative',
                                 'Scaled derivative')
df_rmse_vary$cond__ <- factor(df_rmse_vary$cond__, levels = c('No derivative',
                                                              'Scaled derivative',
                                                              'Unscaled derivative'))

p_rmse_vary <- ggplot(df_rmse_vary, aes(x = effect1__, y = estimate__, color = effect2__)) +
  theme_bw(base_size=65,
           base_family = 'Times') +
  geom_point(size = 5.5 ,
             position = position_dodge(width = 0.7)) +
  geom_errorbar(aes(ymin = lower__, ymax = upper__),
                width = 0.5,
                linewidth = 1.5,
                position = position_dodge(width = 0.7)) +
  geom_point(aes(x=0.5, y=naive_rmse), colour = '#D55E00', size = 5.5) +
  annotate('text', x = 0.7, y = naive_rmse + 0.01, label = 'Prior', size = 10, colour = '#D55E00') +
  facet_wrap(~cond__) +
  labs(x = 'Number of output dimensions', y = 'RMSE', colour = 'Varying hyperparameters') +
  guides(fill = 'none') + 
  theme(axis.ticks = element_line(linewidth = 3)) +
  scale_colour_manual(values = c("#E69F00", "#0072B2")) + ggtitle('(b)')

## Plot for effect of correlated outputs RMSE
df_rmse_corr <- as.data.frame(m_rmse_corr$`d:corr`)
df_rmse_corr$cond__ <- ordered(df_rmse_corr$cond__)
levels(df_rmse_corr$cond__) <- c('No derivative',
                                 'Unscaled derivative',
                                 'Scaled derivative')
df_rmse_corr$cond__ <- factor(df_rmse_corr$cond__, levels = c('No derivative',
                                                              'Scaled derivative',
                                                              'Unscaled derivative'))

p_rmse_corr <- ggplot(df_rmse_corr, aes(x = effect1__, y = estimate__, color = effect2__)) +
  theme_bw(base_size=65,
           base_family = 'Times') +
  geom_point(size = 5.5 ,
             position = position_dodge(width = 0.7)) +
  geom_errorbar(aes(ymin = lower__, ymax = upper__),
                width = 0.5,
                linewidth = 1.5,
                position = position_dodge(width = 0.7)) +
  geom_point(aes(x=0.5, y=naive_rmse), colour = '#D55E00', size = 5.5) +
  annotate('text', x = 0.7, y = naive_rmse + 0.01, label = 'Prior', size = 10, colour = '#D55E00') +
  facet_wrap(~cond__) +
  labs(x = 'Number of output dimensions', y = 'RMSE', colour = 'Correlated outputs') +
  guides(fill = 'none') + 
  theme(axis.ticks = element_line(linewidth = 3)) +
  scale_colour_manual(values = c("#CC79A7", "#009E73")) + ggtitle('(c)')

sim_rmse_eff_full <- p_rmse_main + p_rmse_vary + p_rmse_corr +
  plot_layout(axes = 'collect', guides = 'collect') & 
  theme(axis.text = element_text(size = 75),legend.position = 'bottom', 
        axis.title = element_text(size = 75), legend.text = element_text(size = 75),
        legend.title = element_text(size = 75))

ggsave('latentx_rmse_full.pdf',
       sim_rmse_eff_full,
       width = 165,
       height = 40,
       units = 'cm',
       dpi = 300,
       limitsize = FALSE)

## Plot for main effects MAE
df_mae_main <- as.data.frame(m_mae_main$`d`)
df_mae_main$cond__ <- ordered(df_mae_main$cond__)
levels(df_mae_main$cond__) <- c('No derivative',
                                'Unscaled derivative',
                                'Scaled derivative')
df_mae_main$cond__ <- factor(df_mae_main$cond__, levels = c('No derivative',
                                                            'Scaled derivative',
                                                            'Unscaled derivative'))

p_mae_main <- ggplot(df_mae_main, aes(x = effect1__, y = estimate__)) +
  theme_bw(base_size = 65,
           base_family = 'Times') +
  geom_point(size = 5.5,
             position = position_dodge(width = 0.7)) +
  geom_errorbar(aes(ymin = lower__, ymax = upper__),
                width = 0.3,
                linewidth = 1,
                position = position_dodge(width = 0.7)) +
  geom_point(aes(x=0.5, y=naive_mae), colour = '#D55E00', size = 5.5) +
  annotate('text', x = 0.7, y = naive_mae + 0.01, label = 'Prior', size = 10, colour = '#D55E00') +
  facet_wrap(~cond__) +
  labs(x = 'Number of output dimensions', y = 'MAE') +
  guides(fill = 'none') + 
  theme(axis.ticks = element_line(linewidth = 3)) +
  scale_colour_manual(values = cbbPalette) + ggtitle('(a)')

## Plot for effect of varying hyperparameters MAE
df_mae_vary <- as.data.frame(m_mae_vary$`d:vary`)
df_mae_vary$cond__ <- ordered(df_mae_vary$cond__)
levels(df_mae_vary$cond__) <- c('No derivative',
                                'Unscaled derivative',
                                'Scaled derivative')
df_mae_vary$cond__ <- factor(df_mae_vary$cond__, levels = c('No derivative',
                                                            'Scaled derivative',
                                                            'Unscaled derivative'))

p_mae_vary <- ggplot(df_mae_vary, aes(x = effect1__, y = estimate__, color = effect2__)) +
  theme_bw(base_size=65,
           base_family = 'Times') +
  geom_point(size = 5.5 ,
             position = position_dodge(width = 0.7)) +
  geom_errorbar(aes(ymin = lower__, ymax = upper__),
                width = 0.5,
                linewidth = 1.5,
                position = position_dodge(width = 0.7)) +
  geom_point(aes(x=0.5, y= naive_mae), colour = '#D55E00', size = 5.5) +
  annotate('text', x = 0.7, y = naive_mae + 0.01, label = 'Prior', size = 10, colour = '#D55E00') +
  facet_wrap(~cond__) +
  labs(x = 'Number of output dimensions', y = 'MAE', colour = 'Varying hyperparameters') +
  guides(fill = 'none') + 
  theme(axis.ticks = element_line(linewidth = 3)) +
  scale_colour_manual(values = c("#E69F00", "#0072B2")) + ggtitle('(b)')

## Plot for effect of correlated outputs MAE
df_mae_corr <- as.data.frame(m_mae_corr$`d:corr`)
df_mae_corr$cond__ <- ordered(df_mae_corr$cond__)
levels(df_mae_corr$cond__) <- c('No derivative',
                                'Unscaled derivative',
                                'Scaled derivative')
df_mae_corr$cond__ <- factor(df_mae_corr$cond__, levels = c('No derivative',
                                                            'Scaled derivative',
                                                            'Unscaled derivative'))

p_mae_corr <- ggplot(df_mae_corr, aes(x = effect1__, y = estimate__, color = effect2__)) +
  theme_bw(base_size=65,
           base_family = 'Times') +
  geom_point(size = 5.5 ,
             position = position_dodge(width = 0.7)) +
  geom_errorbar(aes(ymin = lower__, ymax = upper__),
                width = 0.5,
                linewidth = 1.5,
                position = position_dodge(width = 0.7)) +
  geom_point(aes(x=0.5, y=naive_mae), colour = '#D55E00', size = 5.5) +
  annotate('text', x = 0.7, y = naive_mae + 0.01, label = 'Prior', size = 10, colour = '#D55E00') +
  facet_wrap(~cond__) +
  labs(x = 'Number of output dimensions', y = 'MAE', colour = 'Correlated outputs') +
  guides(fill = 'none') + 
  theme(axis.ticks = element_line(linewidth = 3)) +
  scale_colour_manual(values = c("#CC79A7", "#009E73")) + ggtitle('(c)')

sim_mae_eff_full <- p_mae_main + p_mae_vary + p_mae_corr + 
  plot_layout(axes = 'collect', guides = 'collect') & 
  theme(axis.text = element_text(size = 75),legend.position = 'bottom', 
        axis.title = element_text(size = 75), legend.text = element_text(size = 75),
        legend.title = element_text(size = 75))

ggsave('latentx_mae_full.pdf',
       sim_mae_eff_full,
       width = 165,
       height = 40,
       units = 'cm',
       dpi = 300,
       limitsize = FALSE)

# Plot for main effects RMSE 
df_rmse_main <- as.data.frame(m_rmse_main$`d`)
df_rmse_main$cond__ <- ordered(df_rmse_main$cond__)
levels(df_rmse_main$cond__) <- c('No derivative',
                                 'Unscaled derivative',
                                 'Scaled derivative')
df_rmse_main$cond__ <- factor(df_rmse_main$cond__, levels = c('No derivative',
                                                              'Scaled derivative',
                                                              'Unscaled derivative'))
df_rmse_main <- df_rmse_main[df_rmse_main$cond__!='Unscaled derivative',]
p_rmse_main <- ggplot(df_rmse_main, aes(x = effect1__, y = estimate__)) +
  theme_bw(base_size=65,
           base_family = 'Times') +
  geom_point(size = 5.5,
             position = position_dodge(width = 0.7)) +
  geom_errorbar(aes(ymin = lower__, ymax = upper__),
                width = 0.3,
                linewidth = 1,
                position = position_dodge(width = 0.7)) +
  geom_point(aes(x=0.5, y=naive_rmse), colour = '#D55E00', size = 5.5) +
  annotate('text', x = 0.7, y = naive_rmse + 0.01, label = 'Prior', size = 10, colour = '#D55E00') +
  facet_wrap(~cond__) +
  labs(x = 'Number of output dimensions', y = 'RMSE') +
  guides(fill = 'none') + 
  theme(axis.ticks = element_line(linewidth = 3)) +
  scale_colour_manual(values = cbbPalette) + ggtitle('(a)')

## Plot for effect of varying hyperparameters RMSE
df_rmse_vary <- as.data.frame(m_rmse_vary$`d:vary`)
df_rmse_vary$cond__ <- ordered(df_rmse_vary$cond__)
levels(df_rmse_vary$cond__) <- c('No derivative',
                                 'Unscaled derivative',
                                 'Scaled derivative')
df_rmse_vary$cond__ <- factor(df_rmse_vary$cond__, levels = c('No derivative',
                                                              'Scaled derivative',
                                                              'Unscaled derivative'))
df_rmse_vary <- df_rmse_vary[df_rmse_vary$cond__!='Unscaled derivative',]
p_rmse_vary <- ggplot(df_rmse_vary, aes(x = effect1__, y = estimate__, color = effect2__)) +
  theme_bw(base_size=65,
           base_family = 'Times') +
  geom_point(size = 5.5 ,
             position = position_dodge(width = 0.7)) +
  geom_errorbar(aes(ymin = lower__, ymax = upper__),
                width = 0.5,
                linewidth = 1.5,
                position = position_dodge(width = 0.7)) +
  geom_point(aes(x=0.5, y=naive_rmse), colour = '#D55E00', size = 5.5) +
  annotate('text', x = 0.7, y = naive_rmse + 0.01, label = 'Prior', size = 10, colour = '#D55E00') +
  facet_wrap(~cond__) +
  labs(x = 'Number of output dimensions', y = 'RMSE', colour = 'Varying hyperparameters') +
  guides(fill = 'none') + 
  theme(axis.ticks = element_line(linewidth = 3)) +
  scale_colour_manual(values = c("#E69F00", "#0072B2")) + ggtitle('(b)')

## Plot for effect of correlated outputs RMSE
df_rmse_corr <- as.data.frame(m_rmse_corr$`d:corr`)
df_rmse_corr$cond__ <- ordered(df_rmse_corr$cond__)
levels(df_rmse_corr$cond__) <- c('No derivative',
                                 'Unscaled derivative',
                                 'Scaled derivative')
df_rmse_corr$cond__ <- factor(df_rmse_corr$cond__, levels = c('No derivative',
                                                              'Scaled derivative',
                                                              'Unscaled derivative'))
df_rmse_corr <- df_rmse_corr[df_rmse_corr$cond__!='Unscaled derivative',]
p_rmse_corr <- ggplot(df_rmse_corr, aes(x = effect1__, y = estimate__, color = effect2__)) +
  theme_bw(base_size=65,
           base_family = 'Times') +
  geom_point(size = 5.5 ,
             position = position_dodge(width = 0.7)) +
  geom_errorbar(aes(ymin = lower__, ymax = upper__),
                width = 0.5,
                linewidth = 1.5,
                position = position_dodge(width = 0.7)) +
  geom_point(aes(x=0.5, y=naive_rmse), colour = '#D55E00', size = 5.5) +
  annotate('text', x = 0.7, y = naive_rmse + 0.01, label = 'Prior', size = 10, colour = '#D55E00') +
  facet_wrap(~cond__) +
  labs(x = 'Number of output dimensions', y = 'RMSE', colour = 'Correlated outputs') +
  guides(fill = 'none') + 
  theme(axis.ticks = element_line(linewidth = 3)) +
  scale_colour_manual(values = c("#CC79A7", "#009E73")) + ggtitle('(c)')

sim_rmse_eff <- p_rmse_main + p_rmse_vary + p_rmse_corr +
  plot_layout(axes = 'collect', guides = 'collect') & 
  theme(axis.text = element_text(size = 70),legend.position = 'bottom', 
        axis.title = element_text(size = 70), legend.text = element_text(size = 70),
        legend.title = element_text(size = 70))

ggsave('latentx_rmse.pdf',
       sim_rmse_eff,
       width = 120,
       height = 40,
       units = 'cm',
       dpi = 300,
       limitsize = FALSE)

## Plot for main effects MAE
df_mae_main <- as.data.frame(m_mae_main$`d`)
df_mae_main$cond__ <- ordered(df_mae_main$cond__)
levels(df_mae_main$cond__) <- c('No derivative',
                                'Unscaled derivative',
                                'Scaled derivative')
df_mae_main$cond__ <- factor(df_mae_main$cond__, levels = c('No derivative',
                                                            'Scaled derivative',
                                                            'Unscaled derivative'))
df_mae_main <- df_mae_main[df_mae_main$cond__!='Unscaled derivative',]
p_mae_main <- ggplot(df_mae_main, aes(x = effect1__, y = estimate__)) +
  theme_bw(base_size = 65,
           base_family = 'Times') +
  geom_point(size = 5.5,
             position = position_dodge(width = 0.7)) +
  geom_errorbar(aes(ymin = lower__, ymax = upper__),
                width = 0.3,
                linewidth = 1,
                position = position_dodge(width = 0.7)) +
  geom_point(aes(x=0.5, y=naive_mae), colour = '#D55E00', size = 5.5) +
  annotate('text', x = 0.7, y = naive_mae + 0.01, label = 'Prior', size = 10, colour = '#D55E00') +
  facet_wrap(~cond__) +
  labs(x = 'Number of output dimensions', y = 'MAE') +
  guides(fill = 'none') + 
  theme(axis.ticks = element_line(linewidth = 3)) +
  scale_colour_manual(values = cbbPalette) + ggtitle('(a)')

## Plot for effect of varying hyperparameters MAE
df_mae_vary <- as.data.frame(m_mae_vary$`d:vary`)
df_mae_vary$cond__ <- ordered(df_mae_vary$cond__)
levels(df_mae_vary$cond__) <- c('No derivative',
                                'Unscaled derivative',
                                'Scaled derivative')
df_mae_vary$cond__ <- factor(df_mae_vary$cond__, levels = c('No derivative',
                                                            'Scaled derivative',
                                                            'Unscaled derivative'))
df_mae_vary <- df_mae_vary[df_mae_vary$cond__!='Unscaled derivative',]
p_mae_vary <- ggplot(df_mae_vary, aes(x = effect1__, y = estimate__, color = effect2__)) +
  theme_bw(base_size=65,
           base_family = 'Times') +
  geom_point(size = 5.5 ,
             position = position_dodge(width = 0.7)) +
  geom_errorbar(aes(ymin = lower__, ymax = upper__),
                width = 0.5,
                linewidth = 1.5,
                position = position_dodge(width = 0.7)) +
  geom_point(aes(x=0.5, y= naive_mae), colour = '#D55E00', size = 5.5) +
  annotate('text', x = 0.7, y = naive_mae + 0.01, label = 'Prior', size = 10, colour = '#D55E00') +
  facet_wrap(~cond__) +
  labs(x = 'Number of output dimensions', y = 'MAE', colour = 'Varying hyperparameters') +
  guides(fill = 'none') + 
  theme(axis.ticks = element_line(linewidth = 3)) +
  scale_colour_manual(values = c("#E69F00", "#0072B2")) + ggtitle('(b)')

## Plot for effect of correlated outputs MAE
df_mae_corr <- as.data.frame(m_mae_corr$`d:corr`)
df_mae_corr$cond__ <- ordered(df_mae_corr$cond__)
levels(df_mae_corr$cond__) <- c('No derivative',
                                'Unscaled derivative',
                                'Scaled derivative')
df_mae_corr$cond__ <- factor(df_mae_corr$cond__, levels = c('No derivative',
                                                            'Scaled derivative',
                                                            'Unscaled derivative'))
df_mae_corr <- df_mae_corr[df_mae_corr$cond__!='Unscaled derivative',]
p_mae_corr <- ggplot(df_mae_corr, aes(x = effect1__, y = estimate__, color = effect2__)) +
  theme_bw(base_size=65,
           base_family = 'Times') +
  geom_point(size = 5.5 ,
             position = position_dodge(width = 0.7)) +
  geom_errorbar(aes(ymin = lower__, ymax = upper__),
                width = 0.5,
                linewidth = 1.5,
                position = position_dodge(width = 0.7)) +
  geom_point(aes(x=0.5, y=naive_mae), colour = '#D55E00', size = 5.5) +
  annotate('text', x = 0.7, y = naive_mae + 0.01, label = 'Prior', size = 10, colour = '#D55E00') +
  facet_wrap(~cond__) +
  labs(x = 'Number of output dimensions', y = 'MAE', colour = 'Correlated outputs') +
  guides(fill = 'none') + 
  theme(axis.ticks = element_line(linewidth = 3)) +
  scale_colour_manual(values = c("#CC79A7", "#009E73")) + ggtitle('(c)')

sim_mae_eff <- p_mae_main + p_mae_vary + p_mae_corr + 
  plot_layout(axes = 'collect', guides = 'collect') & 
  theme(axis.text = element_text(size = 70),legend.position = 'bottom', 
        axis.title = element_text(size = 70), legend.text = element_text(size = 70),
        legend.title = element_text(size = 70))

ggsave('latentx_mae.pdf',
       sim_mae_eff,
       width = 120,
       height = 40,
       units = 'cm',
       dpi = 300,
       limitsize = FALSE)

