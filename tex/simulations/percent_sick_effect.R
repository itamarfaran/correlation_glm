source('tex/simulations/aux_.R')

n_sim <- 200
n <- 80
p <- 32
percent_alpha <- 0

examples <- expand.grid(p_s = c(0.2, 0.35, 0.5, 0.65, 0.8), ARMA = c(0, 0.5))

toplot <- do.call(rbind, pbmclapply(seq_len(nrow(examples)), function(i) create_variance_estimates(
  n_sim = n_sim, n = n, p = p, percent_alpha = percent_alpha,
  range_alpha = c(1, 1), p_s = examples[i, 1], ARMA = examples[i, 2]), mc.cores = ncores))

toplot[,`:=`(ratio = est/emp)]
p <- ggplot(toplot, aes(x = p_s, y = ratio, group = p_s)) + 
  geom_hline(yintercept = 1, size = 1, linetype = 3) + 
  geom_boxplot(alpha = 0.2, fill = 'grey', outlier.shape = NA) + 
  scale_y_continuous(labels = function(x) paste0('x', x)) + 
  facet_grid(autocorrelated ~ .) + 
  theme_user() + 
  labs(
    title = 'Effect of Relative Group Sizes on the Bias of the SD Estimate',
    x = TeX('$n_{d}/\\left(n_{d}+n_{h}\\right)$'), y = 'Ratio Between Estimated & Empirical\nStandard Deviations'
  )

save(toplot, file='tex/simulations/percent_sick_effect.RData')
custom_ggsave('percent_sick_effect.png', p, width=2)
