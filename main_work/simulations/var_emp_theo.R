source('main_work/simulations/auxilary_functions.R')

n_sim <- 100
range_alpha <- c(0.7, 1.1)

examples <- rbind(
  expand.grid(n = c(40, 80, 120), p = c(25, 40, 55), percent_alpha = c(0), ARMA = c(0)),
  expand.grid(n = c(40, 80, 120), p = c(25, 40, 55), percent_alpha = c(0.4), ARMA = c(0.5))
)
toplot <- do.call(rbind, lapply(1:nrow(examples), function(i) create_variance_estimates(
  n_sim = n_sim, n = examples[i, 1], p = examples[i, 2], p_s = 0.5, percent_alpha = examples[i, 3],
  range_alpha = range_alpha, ARMA = examples[i, 4], ncores = ncores
  )))

toplot[,`:=`(is_null = ifelse(alpha == 1, 'Null', 'Non Null'))]
p1 <- ggplot(toplot[case == 'No Effect' & autocorrelated == 'Not Autocorrelated'],
       aes(x = emp, y = est)) + 
  geom_point(shape = 21, col = 'black', fill = 'darkgrey', alpha = 0.7) + 
  facet_grid(p ~ n, labeller = label_both) + 
  geom_abline(slope = 1, intercept = 0) + 
  theme_user() + 
  labs(
    title = 'GEE Estimate of Variance',
    subtitle = 'Null cases, No auto-correlation',
    x = 'Empirical Variance', y = 'GEE Estimated Variance'
    )


p2 <- ggplot(toplot[case == 'Effect' & autocorrelated == 'Autocorrelated'],
       aes(x = emp, y = est, shape = is_null)) + 
  geom_point(shape = 21, col = 'black', fill = 'darkgrey', alpha = 0.7) + 
  facet_grid(p ~ n, labeller = label_both) + 
  geom_abline(slope = 1, intercept = 0) + 
  theme_user() + 
  theme(legend.position = 'bottom') + 
  labs(
    title = 'GEE Estimate of Variance',
    subtitle = 'Non null cases, With auto-correlation',
    x = 'Empirical Variance', y = 'GEE Estimated Variance',
    shape = '\u03B1 Value:'
  )
save(toplot, file = 'main_work/simulations/var_emp_theo.RData')

custom_ggsave('variance_emp_theo_null.png', p1, width=2, height=1)
custom_ggsave('variance_emp_theo_nonnull.png', p2, width=2, height=1)
