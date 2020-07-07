source('main_work/simulations/auxilary_functions.R')

n_sim <- 1
n <- 100
p <- 50
percent_alpha <- 0.4
range_alpha <- c(0.7, 1.1)


out <- do.call(rbind, list(
  create_estimates(n_sim = n_sim, n = n, p = p, percent_alpha = percent_alpha, range_alpha = range_alpha, ARMA = 0),
  create_estimates(n_sim = n_sim, n = n, p = p, percent_alpha = percent_alpha, range_alpha = range_alpha, ARMA = 0.5)
  ))

p1 <- ggplot(out[type == 'Theta'], aes(x = Parameter, y = Estimate)) + 
  geom_point(alpha = 0.3) + 
  facet_grid(. ~ autocorrelated) + 
  geom_abline(intercept = 0, slope = 1) + 
  ggtitle('Bias of Theta Estimate') + 
  theme_user()

p2 <- ggplot(mapping = aes(x = Parameter, y = Estimate)) + 
  geom_point(data = out[type == 'Alpha' & Value != 'Null'], shape = 17) + 
  geom_boxplot(data = out[type == 'Alpha' & Value == 'Null'], width = 0.1, outlier.shape = 1) + 
  facet_grid(. ~ autocorrelated) + 
  geom_abline(intercept = 0, slope = 1) + 
  ggtitle('Bias of Alpha Estimate', 'Showing Box Plot for Null Values (\u03B1 = 1)') + 
  theme_user()

out2 <- arrangeGrob(p1, p2, nrow=2)
custom_ggsave('bias_small_sample.png', out2, 1.5, 1.2)

