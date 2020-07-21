source('main_work/simulations/auxilary_functions.R')

n_sim <- 1
n <- 100
p <- 50
percent_alpha <- 0.8
range_alpha <- c(0.7, 1.1)


out <- do.call(rbind, list(
  create_estimates(n_sim = n_sim, n = n, p = p, percent_alpha = percent_alpha, range_alpha = range_alpha, ARMA = 0),
  create_estimates(n_sim = n_sim, n = n, p = p, percent_alpha = percent_alpha, range_alpha = range_alpha, ARMA = 0.5)
  ))

p1 <- ggplot(out[type == 'Theta'], aes(x = Parameter, y = Estimate)) + 
  geom_point(alpha = 0.2, shape = 21, color = 'black', fill = 'grey') + 
  facet_grid(. ~ autocorrelated) + 
  geom_abline(intercept = 0, slope = 1) + 
  ggtitle('Real Theta vs. Estimate') + 
  theme_user()

p2 <- ggplot(out[type == 'Alpha'], mapping = aes(x = Parameter, y = Estimate, shape = Value)) + 
  geom_point(alpha = 0.6, color = 'black', fill = 'grey') + 
  scale_shape_manual(values = c('Null'=23, 'Non-Null'=16)) + 
  facet_grid(. ~ autocorrelated) + 
  geom_abline(intercept = 0, slope = 1) + 
  ggtitle('Real Alpha vs. Estimate') + 
  theme_user() + 
  theme(legend.position = 'bottom')

out2 <- arrangeGrob(p1, p2, nrow=2)
custom_ggsave('bias_small_sample.png', out2, 1.5, 1.2)

