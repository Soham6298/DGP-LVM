# We show and check MCMC convergence for DGP-LVM and other fitted models for the simulated data scenarios
# We check Rhats, Bulk-ESS and Tail-ESS for fitted models

# Libraries
library(dplyr)
library(tidyverse)
library(bayesplot)
library(ggplot2)
library(patchwork)
library(lemon)
library(grid)
library(gtable)
library(gridExtra)
## Set variable names, sample points and number of trials
# Source functions
source('DgplvmSimFns.R')
## Set variable names, sample points and number of trials
simdata_out <- readRDS('SimDataOutput.rds')  #change according to output file name from SimStudy
simdata_out$vary <- as.factor(simdata_out$vary)
simdata_out$corr <- as.factor(simdata_out$corr)
simdata_out$deriv_scale <- as.factor(paste(simdata_out$deriv, simdata_out$scale))
simdata_out$d <- as.factor(simdata_out$d)
levels(simdata_out$deriv_scale) <- c('Observed only','Unscaled derivative','Scaled derivative')
levels(simdata_out$vary) <- c('No','Yes')
levels(simdata_out$corr) <- c('No','Yes')
simdata_out$model_id <- as.factor(simdata_out$model_id)
simdata_out$model_id1 <- as.factor(paste0(simdata_out$d, '_',simdata_out$deriv, '_',
                                          simdata_out$scale, '_', simdata_out$vary, '_',
                                          simdata_out$corr))
simdata_out$sim_id <- as.factor(simdata_out$sim_id)

# Convergence for latent x
GP_comb_x <- subset(simdata_out, class == 'x')
GPsummarydata <- GP_comb_x %>%
  group_by(model_id) %>%
  summarise(sim_id = first(sim_id),
            n = first(n),
            d = first(d),
            mtrue_value = mean(true_value),
            mmean = mean(mean), 
            msd = mean(sd), 
            mabs_bias = mean(abs_bias),
            mrmse = mean(rmse), 
            mmae = mean(mae),
            mrhat = mean(rhat),
            mranks = mean(ranks),
            mbess = mean(bess),
            mtess = mean(tess), 
            mruntime = mean(runtime), 
            mdivergent = mean(divergent),
            deriv = first(deriv), 
            scale = first(scale), 
            vary = first(vary),
            corr = first(corr))
GPsummarydata$d<- as.factor(GPsummarydata$d)
GPsummarydata$d <- factor(GPsummarydata$d, levels = c(2, 5, 10))
GPsummarydata$rhat_name <- 'Rhat'
GPsummarydata$bess_name <- 'Bulk-ESS'
GPsummarydata$tess_name <- 'Tail-ESS'
# Rhat plot
x_rhat_summary_plot <- ggplot(GPsummarydata, aes(x = d, y = mrhat)) + 
  theme_bw(base_size = 35, base_family = 'Times') + 
  #geom_boxplot(linewidth = 1) + 
  geom_violin(linewidth = 1) +
  geom_jitter(width = 0.1) + 
  scale_y_log10() +
  facet_wrap(~rhat_name) +
  labs(x = 'Output dimensions', y = 'Values') + 
  theme(axis.ticks = element_line(linewidth = 3), legend.position = 'none') +
  ggtitle('(a)')
# Bulk-ESS plot
x_ess_bulk_summary_plot <- ggplot(GPsummarydata, aes(x = d, y = mbess)) + 
  theme_bw(base_size = 35, base_family = 'Times') +
  #geom_boxplot(linewidth = 1) + 
  geom_violin(linewidth = 1) +
  geom_jitter(width = 0.1) + 
  scale_y_log10() +
  facet_wrap(~bess_name) +
  labs(x = 'Output dimensions', y = 'Values') + 
  theme(axis.ticks = element_line(linewidth = 3), legend.position = 'none') 
# Tail-ESS plot
x_ess_tail_summary_plot <- ggplot(GPsummarydata, aes(x = d, y = mtess)) + 
  theme_bw(base_size = 35, base_family = 'Times') +
  #geom_boxplot(linewidth = 1) + 
  geom_violin(linewidth = 1) +
  geom_jitter(width = 0.1) + 
  scale_y_log10() +
  facet_wrap(~tess_name) +
  labs(x = 'Output dimensions', y = 'Values') + 
  theme(axis.ticks = element_line(linewidth = 3),legend.position = 'none') 

# Convergence for hyperparameters
GP_comb_pars <- subset(simdata_out, class != 'x')
GPsummarydata <- GP_comb_pars %>%
  group_by(model_id) %>%
  summarise(sim_id = first(sim_id),
            n = first(n),
            d = first(d),
            mtrue_value = mean(true_value),
            mmean = mean(mean), 
            msd = mean(sd), 
            mabs_bias = mean(abs_bias),
            mrmse = mean(rmse), 
            mmae = mean(mae),
            mrhat = mean(rhat),
            mranks = mean(ranks),
            mbess = mean(bess),
            mtess = mean(tess), 
            mruntime = mean(runtime), 
            mdivergent = mean(divergent),
            deriv = first(deriv), 
            scale = first(scale), 
            vary = first(vary),
            corr = first(corr))
GPsummarydata$d<- as.factor(GPsummarydata$d)
GPsummarydata$d <- factor(GPsummarydata$d, levels = c(2, 5, 10))
GPsummarydata$rhat_name <- 'Rhat'
GPsummarydata$bess_name <- 'Bulk-ESS'
GPsummarydata$tess_name <- 'Tail-ESS'
# Rhat plot
pars_rhat_summary_plot <- ggplot(GPsummarydata, aes(x = d, y = mrhat)) + 
  geom_violin(linewidth = 1) +
  geom_jitter(width = 0.1) + 
  facet_wrap(~rhat_name) +
  scale_y_log10() +
  theme_bw(base_size = 35, base_family = 'Times') + 
  labs(x = 'Output dimensions', y = 'Values') + 
  theme(axis.ticks = element_line(linewidth = 3), legend.position = 'none') +
  ggtitle('(b)')
# Bulk-ESS plot
pars_ess_bulk_summary_plot <- ggplot(GPsummarydata, aes(x = d, y = mbess)) + 
  geom_violin(linewidth = 1) +
  geom_jitter(width = 0.1) + 
  facet_wrap(~bess_name) +
  scale_y_log10() +
  theme_bw(base_size = 35, base_family = 'Times') + 
  labs(x = 'Output dimensions', y = 'Values') + 
  theme(axis.ticks = element_line(linewidth = 3), legend.position = 'none') 
# Tail-ESS plot
pars_ess_tail_summary_plot <- ggplot(GPsummarydata, aes(x = d, y = mtess)) + 
  geom_violin(linewidth = 1) +
  geom_jitter(width = 0.1) + 
  facet_wrap(~tess_name) +
  scale_y_log10() +
  theme_bw(base_size = 35, base_family = 'Times') + 
  labs(x = 'Output dimensions', y = 'Values') + 
  theme(axis.ticks = element_line(linewidth = 3),legend.position = 'none') 
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