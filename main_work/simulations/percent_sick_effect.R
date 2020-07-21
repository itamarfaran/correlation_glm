source('main_work/code/01_general_functions.R')
source('main_work/code/02_simulation_functions.R')
source('main_work/code/03_estimation_functions.R')
source('main_work/code/04_inference_functions.R')

n_sim <- 200
n = 80
p = 32
percent_alpha = 0

examples <- expand.grid(p_s = c(0.2, 0.35, 0.5, 0.65, 0.8), ARMA = c(0, 0.5))

toplot <- do.call(rbind, pbmclapply(1:nrow(examples), function(i) create_variance_estimates(
  n_sim = 7, n = 50, p = 25, percent_alpha = 0, range_alpha = 0, p_s = examples[i, 1], ARMA = examples[i, 2]), mc.cores = ncores))

toplot[,`:=`(ratio = est/emp)]
p <- ggplot(toplot, aes(x = p_s, y = ratio, group = p_s)) + 
  geom_hline(yintercept = 1, size = 1, linetype = 3) + 
  geom_boxplot(alpha = 0.2, fill = 'grey') + 
  scale_y_continuous(labels = function(x) paste0('x', x)) + 
  facet_grid(autocorrelated ~ .) + 
  theme_user() + 
  labs(
    title = 'Effect of % Diagnosed Subjects on Variance Estimation',
    x = TeX('$n_{d}/\\left(n_{d}+n_{h}\\right)$'), y = 'Ratio Between Estimated & Empirical\nStandard Deviations'
  )

custom_ggsave('percent_sick_effect.png', p, width=2)
