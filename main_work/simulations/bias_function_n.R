source('main_work/simulations/auxilary_functions.R')

n_sim <- 1
percent_alpha <- 0.7
range_alpha <- c(0.7, 1.1)

examples <- expand.grid(n = c(20, 40, 60, 80, 100), p = c(10, 20, 30, 40))
# examples <- expand.grid(n = c(20, 40), p = c(10, 20))

toplot <- do.call(rbind, lapply(1:nrow(examples), function(i) create_estimates(
  n_sim = n_sim, n = examples[i, 1], p = examples[i, 2], percent_alpha = 0, range_alpha = range_alpha, ARMA = 0,
  verbose = FALSE
  )))

toplot[,Bias := Estimate - Parameter]
toplot_samp <- toplot[,.SD[sample(.N, min(.N, 20))],by = .(n, p, type)]

pow = 1.414
S_sqrt <- function(x) sign(x)*abs(x)^(1/pow)
IS_sqrt <- function(x) (abs(x)^pow)*sign(x)
S_sqrt_trans <- function() scales::trans_new("S_sqrt",S_sqrt,IS_sqrt)


p <- ggplot(mapping = aes(x = n, y = Bias, group = n)) + 
  geom_hline(yintercept = 0) +
  geom_point(data=toplot_samp, position = 'jitter', alpha = 0.4) +
  geom_boxplot(data=toplot, alpha = 0.33) +
  scale_y_continuous(trans="S_sqrt") +
  facet_grid(p ~ type) +
  labs(title = 'Bias as a Function of Sample Size and Dimension of Theta', x = 'Sample Size') + 
  theme_user()

custom_ggsave('bias_function_n.png', p, width = 2, height = 1.2)
