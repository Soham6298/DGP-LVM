setwd("F:/Work/DGP/plots")
m_gp_vary <- conditional_effects(m_gp, effects = 'model_dim:is_vary', 
                                 conditions = make_conditions(m_gp, 'deriv_scale'),
                                 resolution = 300)
m_gp_corr <- conditional_effects(m_gp, effects = 'model_dim:is_corr', 
                                 conditions = make_conditions(m_gp, 'deriv_scale'),
                                 resolution = 300)
m_per_vary <- conditional_effects(m_per, effects = 'model_dim:is_vary', 
                                  conditions = make_conditions(m_per, 'deriv_scale'),
                                  resolution = 300)
m_gp_main <- conditional_effects(m_gp, effects = 'model_dim', 
                                 conditions = make_conditions(m_gp, 'deriv_scale'),
                                 resolution = 300)
m_per_main<- conditional_effects(m_per, effects = 'model_dim', 
                                 conditions = make_conditions(m_per, 'deriv_scale'),
                                 resolution = 300)
## gp_vary_full
df_gp_vary <- as.data.frame(m_gp_vary$`model_dim:is_vary`)
df_gp_vary$cond__ <- ordered(df_gp_vary$cond__)
levels(df_gp_vary$cond__) <- c('Derivative:No|Scaling:No',
                       'Derivative:Yes|Scaling:No',
                       'Derivative:Yes|Scaling:Yes')
df_gp_vary$cond__ <- factor(df_gp_vary$cond__, levels = c('Derivative:No|Scaling:No',
                                          'Derivative:Yes|Scaling:Yes',
                                          'Derivative:Yes|Scaling:No'))

p_gp_vary <- ggplot(df_gp_vary, aes(x = effect1__, y = estimate__, color = effect2__)) +
     theme_bw(base_size=30) +
     geom_point(size = 4 ,
                position = position_dodge(width = 0.5)) +
     geom_errorbar(aes(ymin = lower__, ymax = upper__),
                   width = 0.3,
                   linewidth = 1,
                   position = position_dodge(width = 0.5)) +
     facet_wrap(~cond__) +
     labs(x = 'Output dimensions', y = 'RMSE', colour = 'Varying hyperparameters') +
     guides(fill = 'none') + 
     theme(legend.position = 'bottom') +
     scale_colour_manual(values = cbbPalette)
ggsave('gp_vary_full.png',
       p_gp_vary,
       width = 40,
       height = 25,
       units = 'cm',
       dpi = 300)
       
## gp_corr_full
df_gp_corr <- as.data.frame(m_gp_corr$`model_dim:is_corr`)
df_gp_corr$cond__ <- ordered(df_gp_corr$cond__)
levels(df_gp_corr$cond__) <- c('Derivative:No|Scaling:No',
                               'Derivative:Yes|Scaling:No',
                               'Derivative:Yes|Scaling:Yes')
df_gp_corr$cond__ <- factor(df_gp_corr$cond__, levels = c('Derivative:No|Scaling:No',
                                                  'Derivative:Yes|Scaling:Yes',
                                                  'Derivative:Yes|Scaling:No'))

p_gp_corr <- ggplot(df_gp_corr, aes(x = effect1__, y = estimate__, color = effect2__)) +
  theme_bw(base_size=30) +
  geom_point(size = 4 ,
             position = position_dodge(width = 0.5)) +
  geom_errorbar(aes(ymin = lower__, ymax = upper__),
                width = 0.3,
                linewidth = 1,
                position = position_dodge(width = 0.5)) +
  facet_wrap(~cond__) +
  labs(x = 'Output dimensions', y = 'RMSE', colour = 'Correlated outputs') +
  guides(fill = 'none') + 
  theme(legend.position = 'bottom') +
  scale_colour_manual(values = cbbPalette)
ggsave('gp_corr_full.png',
       p_gp_corr,
       width = 40,
       height = 25,
       units = 'cm',
       dpi = 300)
## per_vary_full
df_per_vary <- as.data.frame(m_per_vary$`model_dim:is_vary`)
df_per_vary$cond__ <- ordered(df_per_vary$cond__)
levels(df_per_vary$cond__) <- c('Derivative:No|Scaling:No',
                               'Derivative:Yes|Scaling:No',
                               'Derivative:Yes|Scaling:Yes')
df_per_vary$cond__ <- factor(df_per_vary$cond__, levels = c('Derivative:No|Scaling:No',
                                                          'Derivative:Yes|Scaling:Yes',
                                                          'Derivative:Yes|Scaling:No'))

p_per_vary <- ggplot(df_per_vary, aes(x = effect1__, y = estimate__, color = effect2__)) +
  theme_bw(base_size=30) +
  geom_point(size = 4,
             position = position_dodge(width = 0.5)) +
  geom_errorbar(aes(ymin = lower__, ymax = upper__),
                width = 0.3,
                linewidth = 1,
                position = position_dodge(width = 0.5)) +
  facet_wrap(~cond__) +
  labs(x = 'Output dimensions', y = 'RMSE', colour = 'Varying hyperparameters') +
  guides(fill = 'none') + 
  theme(legend.position = 'bottom') +
  scale_colour_manual(values = cbbPalette)
ggsave('per_vary_full.png',
       p_per_vary,
       width = 40,
       height = 25,
       units = 'cm',
       dpi = 300)
## gp_main_full
df_gp_main <- as.data.frame(m_gp_main$`model_dim`)
df_gp_main$cond__ <- ordered(df_gp_main$cond__)
levels(df_gp_main$cond__) <- c('Derivative:No|Scaling:No',
                                'Derivative:Yes|Scaling:No',
                                'Derivative:Yes|Scaling:Yes')
df_gp_main$cond__ <- factor(df_gp_main$cond__, levels = c('Derivative:No|Scaling:No',
                                                            'Derivative:Yes|Scaling:Yes',
                                                            'Derivative:Yes|Scaling:No'))

p_gp_main <- ggplot(df_gp_main, aes(x = effect1__, y = estimate__)) +
  theme_bw(base_size=30) +
  geom_point(size = 4,
             position = position_dodge(width = 0.5)) +
  geom_errorbar(aes(ymin = lower__, ymax = upper__),
                width = 0.3,
                linewidth = 1,
                position = position_dodge(width = 0.5)) +
  facet_wrap(~cond__) +
  labs(x = 'Output dimensions', y = 'RMSE') +
  guides(fill = 'none') + 
  theme(legend.position = 'bottom') +
  scale_colour_manual(values = cbbPalette)
ggsave('gp_main_full.png',
       p_gp_main,
       width = 40,
       height = 25,
       units = 'cm',
       dpi = 300)
## per_main_full
df_per_main <- as.data.frame(m_per_main$`model_dim`)
df_per_main$cond__ <- ordered(df_per_main$cond__)
levels(df_per_main$cond__) <- c('Derivative:No|Scaling:No',
                               'Derivative:Yes|Scaling:No',
                               'Derivative:Yes|Scaling:Yes')
df_per_main$cond__ <- factor(df_per_main$cond__, levels = c('Derivative:No|Scaling:No',
                                                          'Derivative:Yes|Scaling:Yes',
                                                          'Derivative:Yes|Scaling:No'))

p_per_main <- ggplot(df_per_main, aes(x = effect1__, y = estimate__)) +
  theme_bw(base_size=30) +
  geom_point(size = 4,
             position = position_dodge(width = 0.5)) +
  geom_errorbar(aes(ymin = lower__, ymax = upper__),
                width = 0.3,
                linewidth = 1,
                position = position_dodge(width = 0.5)) +
  facet_wrap(~cond__) +
  labs(x = 'Output dimensions', y = 'RMSE') +
  guides(fill = 'none') + 
  theme(legend.position = 'bottom') +
  scale_colour_manual(values = cbbPalette)
ggsave('per_main_full.png',
       p_per_main,
       width = 40,
       height = 25,
       units = 'cm',
       dpi = 300)


## gp_vary
df_gp_vary <- as.data.frame(m_gp_vary$`model_dim:is_vary`)
df_gp_vary$cond__ <- ordered(df_gp_vary$cond__)
levels(df_gp_vary$cond__) <- c('Derivative:No|Scaling:No',
                               'Derivative:Yes|Scaling:No',
                               'Derivative:Yes|Scaling:Yes')
df_gp_vary$cond__ <- factor(df_gp_vary$cond__, levels = c('Derivative:No|Scaling:No',
                                                          'Derivative:Yes|Scaling:Yes',
                                                          'Derivative:Yes|Scaling:No'))
df_gp_vary <- df_gp_vary[df_gp_vary$cond__!='Derivative:Yes|Scaling:No',]
p_gp_vary <- ggplot(df_gp_vary, aes(x = effect1__, y = estimate__, color = effect2__)) +
  theme_bw(base_size=30) +
  geom_point(size = 4 ,
             position = position_dodge(width = 0.5)) +
  geom_errorbar(aes(ymin = lower__, ymax = upper__),
                width = 0.3,
                linewidth = 1,
                position = position_dodge(width = 0.5)) +
  facet_wrap(~cond__) +
  labs(x = 'Output dimensions', y = 'RMSE', colour = 'Varying hyperparameters') +
  guides(fill = 'none') + 
  theme(legend.position = 'bottom') +
  scale_colour_manual(values = cbbPalette)
ggsave('gp_vary.png',
       p_gp_vary,
       width = 30,
       height = 20,
       units = 'cm',
       dpi = 300)

## gp_corr
df_gp_corr <- as.data.frame(m_gp_corr$`model_dim:is_corr`)
df_gp_corr$cond__ <- ordered(df_gp_corr$cond__)
levels(df_gp_corr$cond__) <- c('Derivative:No|Scaling:No',
                               'Derivative:Yes|Scaling:No',
                               'Derivative:Yes|Scaling:Yes')
df_gp_corr$cond__ <- factor(df_gp_corr$cond__, levels = c('Derivative:No|Scaling:No',
                                                          'Derivative:Yes|Scaling:Yes',
                                                          'Derivative:Yes|Scaling:No'))
df_gp_corr <- df_gp_corr[df_gp_corr$cond__!='Derivative:Yes|Scaling:No',]
p_gp_corr <- ggplot(df_gp_corr, aes(x = effect1__, y = estimate__, color = effect2__)) +
  theme_bw(base_size=30) +
  geom_point(size = 4 ,
             position = position_dodge(width = 0.5)) +
  geom_errorbar(aes(ymin = lower__, ymax = upper__),
                width = 0.3,
                linewidth = 1,
                position = position_dodge(width = 0.5)) +
  facet_wrap(~cond__) +
  labs(x = 'Output dimensions', y = 'RMSE', colour = 'Correlated outputs') +
  guides(fill = 'none') + 
  theme(legend.position = 'bottom') +
  scale_colour_manual(values = cbbPalette)
ggsave('gp_corr.png',
       p_gp_corr,
       width = 30,
       height = 20,
       units = 'cm',
       dpi = 300)
## per_vary
df_per_vary <- as.data.frame(m_per_vary$`model_dim:is_vary`)
df_per_vary$cond__ <- ordered(df_per_vary$cond__)
levels(df_per_vary$cond__) <- c('Derivative:No|Scaling:No',
                                'Derivative:Yes|Scaling:No',
                                'Derivative:Yes|Scaling:Yes')
df_per_vary$cond__ <- factor(df_per_vary$cond__, levels = c('Derivative:No|Scaling:No',
                                                            'Derivative:Yes|Scaling:Yes',
                                                            'Derivative:Yes|Scaling:No'))
df_per_vary <- df_per_vary[df_per_vary$cond__!='Derivative:Yes|Scaling:No',]
p_per_vary <- ggplot(df_per_vary, aes(x = effect1__, y = estimate__, color = effect2__)) +
  theme_bw(base_size=30) +
  geom_point(size = 4,
             position = position_dodge(width = 0.5)) +
  geom_errorbar(aes(ymin = lower__, ymax = upper__),
                width = 0.3,
                linewidth = 1,
                position = position_dodge(width = 0.5)) +
  facet_wrap(~cond__) +
  labs(x = 'Output dimensions', y = 'RMSE', colour = 'Varying hyperparameters') +
  guides(fill = 'none') + 
  theme(legend.position = 'bottom') +
  scale_colour_manual(values = cbbPalette)
ggsave('per_vary.png',
       p_per_vary,
       width = 30,
       height = 20,
       units = 'cm',
       dpi = 300)
## gp_main
df_gp_main <- as.data.frame(m_gp_main$`model_dim`)
df_gp_main$cond__ <- ordered(df_gp_main$cond__)
levels(df_gp_main$cond__) <- c('Derivative:No|Scaling:No',
                               'Derivative:Yes|Scaling:No',
                               'Derivative:Yes|Scaling:Yes')
df_gp_main$cond__ <- factor(df_gp_main$cond__, levels = c('Derivative:No|Scaling:No',
                                                          'Derivative:Yes|Scaling:Yes',
                                                          'Derivative:Yes|Scaling:No'))
df_gp_main <- df_gp_main[df_gp_main$cond__!='Derivative:Yes|Scaling:No',]
p_gp_main <- ggplot(df_gp_main, aes(x = effect1__, y = estimate__)) +
  theme_bw(base_size=30) +
  geom_point(size = 4,
             position = position_dodge(width = 0.5)) +
  geom_errorbar(aes(ymin = lower__, ymax = upper__),
                width = 0.3,
                linewidth = 1,
                position = position_dodge(width = 0.5)) +
  facet_wrap(~cond__) +
  labs(x = 'Output dimensions', y = 'RMSE') +
  guides(fill = 'none') + 
  theme(legend.position = 'bottom') +
  scale_colour_manual(values = cbbPalette)
ggsave('gp_main.png',
       p_gp_main,
       width = 30,
       height = 20,
       units = 'cm',
       dpi = 300)
## per_main
df_per_main <- as.data.frame(m_per_main$`model_dim`)
df_per_main$cond__ <- ordered(df_per_main$cond__)
levels(df_per_main$cond__) <- c('Derivative:No|Scaling:No',
                                'Derivative:Yes|Scaling:No',
                                'Derivative:Yes|Scaling:Yes')
df_per_main$cond__ <- factor(df_per_main$cond__, levels = c('Derivative:No|Scaling:No',
                                                            'Derivative:Yes|Scaling:Yes',
                                                            'Derivative:Yes|Scaling:No'))
df_per_main <- df_per_main[df_per_main$cond__!='Derivative:Yes|Scaling:No',]
p_per_main <- ggplot(df_per_main, aes(x = effect1__, y = estimate__)) +
  theme_bw(base_size=30) +
  geom_point(size = 4,
             position = position_dodge(width = 0.5)) +
  geom_errorbar(aes(ymin = lower__, ymax = upper__),
                width = 0.3,
                linewidth = 1,
                position = position_dodge(width = 0.5)) +
  facet_wrap(~cond__) +
  labs(x = 'Output dimensions', y = 'RMSE') +
  guides(fill = 'none') + 
  theme(legend.position = 'bottom') +
  scale_colour_manual(values = cbbPalette)
ggsave('per_main.png',
       p_per_main,
       width = 30,
       height = 20,
       units = 'cm',
       dpi = 300)

