## Plots for hyperparameters
# Prepare data for plots
## Set variable names, sample points and number of trials
simdata_out <- readRDS('simulation results/dgplvm_simout_se.rds')  #change according to output file name from SimStudy
simdata_out$model_cond <- paste0(simdata_out$d,'_',simdata_out$deriv, 
                                 simdata_out$scale, simdata_out$vary, simdata_out$corr)
simdata_out$vary <- as.factor(simdata_out$vary)
simdata_out$corr <- as.factor(simdata_out$corr)
simdata_out$deriv_scale <- as.factor(paste(simdata_out$deriv, simdata_out$scale))
simdata_out$d <- as.factor(simdata_out$d)
levels(simdata_out$deriv_scale) <- c('Observed only','Unscaled derivative','Scaled derivative')
levels(simdata_out$vary) <- c('No','Yes')
levels(simdata_out$corr) <- c('No','Yes')
simdata_out$model_id <- as.factor(simdata_out$model_id)
simdata_out$sim_id <- as.factor(simdata_out$sim_id)
dgp_sim_rho <- subset(simdata_out, class =='rho')
dgp_sim_alpha_obs <- subset(simdata_out, class =='alpha_obs')
dgp_sim_alpha_obs$output_type <- 'f'
dgp_sim_alpha_grad <- subset(simdata_out, class =='alpha_grad')
dgp_sim_alpha_grad$output_type <- "f'"
dgp_sim_alpha <- rbind(dgp_sim_alpha_obs, dgp_sim_alpha_grad)
dgp_sim_alpha$output_type <- as.factor(dgp_sim_alpha$output_type)
dgp_sim_alpha$var_names <- 'alpha'
dgp_sim_sigma_obs <- subset(simdata_out, class =='sigma_obs')
dgp_sim_sigma_obs$output_type <- 'f'
dgp_sim_sigma_grad <- subset(simdata_out, class =='sigma_grad')
dgp_sim_sigma_grad$output_type <- "f'"
dgp_sim_sigma <- rbind(dgp_sim_sigma_obs, dgp_sim_sigma_grad)
dgp_sim_sigma$output_type <- as.factor(dgp_sim_sigma$output_type)
dgp_sim_sigma$var_names <- 'sigma'

## Plots of hyperparameters with scaling across output dimensions
dgp_id <- c('2_1111', '5_1111','10_1111')
dgp_sim_rho1111 <- subset(dgp_sim_rho, model_cond %in% dgp_id)
dgp_sim_alpha1111 <- subset(dgp_sim_alpha, model_cond %in% dgp_id)
dgp_sim_sigma1111 <- subset(dgp_sim_sigma, model_cond %in% dgp_id)
# Plot for GP lengthscale rho
dgp_sim_rho_plot <- ggplot(data = dgp_sim_rho1111, aes(x=d, y=rmse)) +
  theme_bw(base_size = 35, base_family = 'Times') + 
  geom_violin(linewidth = 1) +
  facet_wrap(~class) +
  labs(x='Number of output dimensions', y='RMSE') + ggtitle('(a)')
# Plot for GP marginal SD alphas
dgp_sim_alpha_plot <- ggplot(data = dgp_sim_alpha1111, aes(x=d, y=rmse, colour=output_type)) +
  theme_bw(base_size = 35, base_family = 'Times') + 
  geom_violin(linewidth = 1) +
  facet_wrap(~var_names) +
  labs(x='Number of output dimensions', y='RMSE') + theme(axis.ticks = element_line(linewidth=3),
                                                          legend.title = element_blank()) +
  scale_colour_manual(values=c("#0072B2", "#E69F00"))
# Plot for error SD sigmas
dgp_sim_sigma_plot <- ggplot(data = dgp_sim_sigma1111, aes(x=d, y=rmse, colour=output_type)) +
  theme_bw(base_size = 35, base_family = 'Times') + 
  geom_violin(linewidth = 1) +
  facet_wrap(~var_names) +
  labs(x='Number of output dimensions', y='RMSE') + theme(axis.ticks = element_line(linewidth=3), 
                                                          legend.title = element_blank()) +
  scale_colour_manual(values=c("#0072B2", "#E69F00"))
# Combine plots
dgp_sim_scaled_plot <- dgp_sim_rho_plot + dgp_sim_alpha_plot + dgp_sim_sigma_plot +
  plot_layout(axis_titles = 'collect', guides = 'collect') & theme(axis.title.x = element_blank())

## Plots of hyperparameters without scaling output modification
dgp_id <- c('2_1011', '5_1011','10_1011')
dgp_sim_rho1011 <- subset(dgp_sim_rho, model_cond %in% dgp_id)
dgp_sim_alpha1011 <- subset(dgp_sim_alpha, model_cond %in% dgp_id)
dgp_sim_sigma1011 <- subset(dgp_sim_sigma, model_cond %in% dgp_id)
# Plot for GP lengthscale rho
dgp_sim_rho_plot <- ggplot(data = dgp_sim_rho1011, aes(x=d, y=rmse)) +
  theme_bw(base_size = 35, base_family = 'Times') + 
  geom_violin(linewidth = 1) +
  facet_wrap(~class) +
  labs(x='Number of output dimensions', y='RMSE') + ggtitle('(b)')
# Plot for GP marginal SD alphas
dgp_sim_alpha_plot <- ggplot(data = dgp_sim_alpha1011, aes(x=d, y=rmse, colour=output_type)) +
  theme_bw(base_size = 35, base_family = 'Times') + 
  geom_violin(linewidth = 1) +
  facet_wrap(~var_names) +
  labs(x='Number of output dimensions', y='RMSE') + theme(axis.ticks = element_line(linewidth=3),
                                                          legend.title = element_blank()) +
  scale_colour_manual(values=c("#0072B2", "#E69F00"))
# Plot for error SD sigmas
dgp_sim_sigma_plot <- ggplot(data = dgp_sim_sigma1011, aes(x=d, y=rmse, colour=output_type)) +
  theme_bw(base_size = 35, base_family = 'Times') + 
  geom_violin(linewidth = 1) +
  facet_wrap(~var_names) +
  labs(x='Number of output dimensions', y='RMSE') + theme(axis.ticks = element_line(linewidth=3), 
                                                          legend.title = element_blank()) +
  scale_colour_manual(values=c("#0072B2", "#E69F00"))
dgp_sim_noscale_plot <- dgp_sim_rho_plot + dgp_sim_alpha_plot + dgp_sim_sigma_plot +
  plot_layout(axis_titles = 'collect', guides = 'collect') #& theme(axis.title.x = element_blank())
dgp_sim_pars_scale_plot <- (dgp_sim_scaled_plot / dgp_sim_noscale_plot) + 
  plot_layout(axis_titles = 'collect', guides = 'collect')
ggsave('sim_pars_scale_plot.pdf',
       dgp_sim_pars_scale_plot,
       dpi = 300,
       width = 60,
       height = 30,
       units = 'cm')

## Plots of hyperparameters without varying across output dimensions
dgp_id <- c('2_1101', '5_1101','10_1101')
dgp_sim_rho1101 <- subset(dgp_sim_rho, model_cond %in% dgp_id)
dgp_sim_alpha1101 <- subset(dgp_sim_alpha, model_cond %in% dgp_id)
dgp_sim_sigma1101 <- subset(dgp_sim_sigma, model_cond %in% dgp_id)
# Plot for GP lengthscale rho
dgp_sim_rho_plot <- ggplot(data = dgp_sim_rho1101, aes(x=d, y=rmse)) +
  theme_bw(base_size = 35, base_family = 'Times') + 
  geom_violin(linewidth = 1) +
  facet_wrap(~class) +
  labs(x='Number of output dimensions', y='RMSE') + ggtitle('(a)')
# Plot for GP marginal SD alphas
dgp_sim_alpha_plot <- ggplot(data = dgp_sim_alpha1101, aes(x=d, y=rmse, colour=output_type)) +
  theme_bw(base_size = 35, base_family = 'Times') + 
  geom_violin(linewidth = 1) +
  facet_wrap(~var_names) +
  labs(x='Number of output dimensions', y='RMSE') + theme(axis.ticks = element_line(linewidth=3),
                                                          legend.title = element_blank()) +
  scale_colour_manual(values=c("#0072B2", "#E69F00"))
# Plot for error SD alphas
dgp_sim_sigma_plot <- ggplot(data = dgp_sim_sigma1101, aes(x=d, y=rmse, colour=output_type)) +
  theme_bw(base_size = 35, base_family = 'Times') + 
  geom_violin(linewidth = 1) +
  facet_wrap(~var_names) +
  labs(x='Number of output dimensions', y='RMSE') + theme(axis.ticks = element_line(linewidth=3), 
                                                          legend.title = element_blank()) +
  scale_colour_manual(values=c("#0072B2", "#E69F00"))
# Combine plots
dgp_sim_novary_plot <- dgp_sim_rho_plot + dgp_sim_alpha_plot + dgp_sim_sigma_plot +
  plot_layout(axis_titles = 'collect', guides = 'collect') & theme(axis.title.y = element_blank(),
                                                                   legend.position = 'none')
## Plots of hyperparameters without correlated output modification
dgp_id <- c('2_1110', '5_1110','10_1110')
dgp_sim_rho1110 <- subset(dgp_sim_rho, model_cond %in% dgp_id)
dgp_sim_alpha1110 <- subset(dgp_sim_alpha, model_cond %in% dgp_id)
dgp_sim_sigma1110 <- subset(dgp_sim_sigma, model_cond %in% dgp_id)
# Plot for GP lengthscale rho
dgp_sim_rho_plot <- ggplot(data = dgp_sim_rho1110, aes(x=d, y=rmse)) +
  theme_bw(base_size = 35, base_family = 'Times') + 
  geom_violin(linewidth = 1) +
  facet_wrap(~class) +
  labs(x='Number of output dimensions', y='RMSE') + ggtitle('(b)')
# Plot for GP marginal SD alphas
dgp_sim_alpha_plot <- ggplot(data = dgp_sim_alpha1110, aes(x=d, y=rmse, colour=output_type)) +
  theme_bw(base_size = 35, base_family = 'Times') + 
  geom_violin(linewidth = 1) +
  facet_wrap(~var_names) +
  labs(x='Number of output dimensions', y='RMSE') + theme(axis.ticks = element_line(linewidth=3),
                                                          legend.title = element_blank()) +
  scale_colour_manual(values=c("#0072B2", "#E69F00"))
# Plot for error SD alphas
dgp_sim_sigma_plot <- ggplot(data = dgp_sim_sigma1110, aes(x=d, y=rmse, colour=output_type)) +
  theme_bw(base_size = 35, base_family = 'Times') + 
  geom_violin(linewidth = 1) +
  facet_wrap(~var_names) +
  labs(x='Number of output dimensions', y='RMSE') + theme(axis.ticks = element_line(linewidth=3), 
                                                          legend.title = element_blank()) +
  scale_colour_manual(values=c("#0072B2", "#E69F00"))
# combine plots
dgp_sim_nocorr_plot <- dgp_sim_rho_plot + dgp_sim_alpha_plot + dgp_sim_sigma_plot +
  plot_layout(axis_titles = 'collect', guides = 'collect') & theme(axis.title.y = element_blank())
## Plots of hyperparameters without varying across output dimensions as well as without correlated output modification
dgp_id <- c('2_1100', '5_1100','10_1100')
dgp_sim_rho1100 <- subset(dgp_sim_rho, model_cond %in% dgp_id)
dgp_sim_alpha1100 <- subset(dgp_sim_alpha, model_cond %in% dgp_id)
dgp_sim_sigma1100 <- subset(dgp_sim_sigma, model_cond %in% dgp_id)
# Plot for GP lengthscale rho
dgp_sim_rho_plot <- ggplot(data = dgp_sim_rho1100, aes(x=d, y=rmse)) +
  theme_bw(base_size = 35, base_family = 'Times') + 
  geom_violin(linewidth = 1) +
  facet_wrap(~class) +
  labs(x='Number of output dimensions', y='RMSE') + ggtitle('(c)')
# Plot for GP marginal SD alphas
dgp_sim_alpha_plot <- ggplot(data = dgp_sim_alpha1100, aes(x=d, y=rmse, colour=output_type)) +
  theme_bw(base_size = 35, base_family = 'Times') + 
  geom_violin(linewidth = 1) +
  facet_wrap(~var_names) +
  labs(x='Number of output dimensions', y='RMSE') + theme(axis.ticks = element_line(linewidth=3),
                                                          legend.title = element_blank()) +
  scale_colour_manual(values=c("#0072B2", "#E69F00"))
# Plot for error SD alphas
dgp_sim_sigma_plot <- ggplot(data = dgp_sim_sigma1100, aes(x=d, y=rmse, colour=output_type)) +
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
       height = 45,
       units = 'cm')
