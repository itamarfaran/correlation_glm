source('tex/simulations/aux_.R')

n_sim <- 100
range_alpha <- c(0.7, 1.1)

examples <- rbind(
  expand.grid(n = c(40, 80, 120), p = c(25, 40, 55), percent_alpha = c(0), ARMA = c(0)),
  expand.grid(n = c(40, 80, 120), p = c(25, 40, 55), percent_alpha = c(0.4), ARMA = c(0.5))
)
toplot <- do.call(rbind, lapply(seq_len(nrow(examples)), function(i) create_variance_estimates(
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
    title = 'GEE Framework for Variance',
    x = 'Empirical Standard Deviations', y = 'GEE Estimated Standard Deviations'
    )


p2 <- ggplot(toplot[case == 'Effect' & autocorrelated == 'Autocorrelated'],
       aes(x = emp, y = est, shape = is_null)) + 
  geom_point(shape = 21, col = 'black', fill = 'darkgrey', alpha = 0.7) + 
  facet_grid(p ~ n, labeller = label_both) + 
  geom_abline(slope = 1, intercept = 0) + 
  theme_user() + 
  theme(legend.position = 'bottom') + 
  labs(
    title = 'GEE Framework for Variance',
    x = 'Empirical Standard Deviations', y = 'GEE Estimated Standard Deviations'
  )

toplot[,bias:=est - emp]
toplot[,per_bias:=est/emp - 1]
p3 <- ggplot(toplot, aes(x='.',y=per_bias)) +
  geom_hline(yintercept = 0) +
  geom_boxplot(outlier.shape = NA) +
  scale_y_continuous(labels = scales::percent) + 
  theme_user() +
  labs(title=TeX('StD Bias'), x='', y='')

out2 <- arrangeGrob(p2, p3, layout_matrix = matrix(c(1, 1, 1, 2), nr=1))
save(toplot, file = 'tex/simulations/var_emp_theo.RData')

custom_ggsave('variance_emp_theo_null.png', p1, width=2, height=1.2)
custom_ggsave('variance_emp_theo_nonnull.png', out2, width=2, height=1.2)
