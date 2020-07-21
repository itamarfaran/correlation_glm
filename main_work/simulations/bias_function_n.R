source('main_work/simulations/auxilary_functions.R')

n_sim <- 1
percent_alpha <- 0.7
range_alpha <- c(0.7, 1.1)

examples <- expand.grid(n = c(20, 40, 60, 80, 100), p = c(15, 30, 45, 60))
# examples <- expand.grid(n = c(20, 40), p = c(10, 20))
examples <- rbind(examples, examples, examples)

toplot <- do.call(rbind, lapply(1:nrow(examples), function(i) create_estimates(
  n_sim = n_sim, n = examples[i, 1], p = examples[i, 2], percent_alpha = 0, range_alpha = range_alpha, ARMA = 0,
  verbose = FALSE
  )))

toplot_ <- copy(toplot)

toplot[,`:=`(Bias = Estimate - Parameter, RMSE = (Estimate - Parameter)^2)]
# toplot <- melt(
#   toplot,
#   id.vars = c('type', 'Estimate', 'Parameter', 'case', 'autocorrelated', 'n', 'p', 'Value'),
#   variable.name = 'Variable',
#   value.name = 'Value_'
#   )

pow = 1.414
S_sqrt <- function(x) sign(x)*abs(x)^(1/pow)
IS_sqrt <- function(x) (abs(x)^pow)*sign(x)
S_sqrt_trans <- function() scales::trans_new("S_sqrt",S_sqrt,IS_sqrt)


p <- ggplot(toplot[type == 'Alpha'], mapping = aes(x = n, fill = p, y = Bias, group = interaction(n, p))) + 
  geom_hline(yintercept = 0) +
  geom_boxplot() +
  scale_y_continuous(trans="S_sqrt") +
  scale_fill_gradient(low='#606060', high='#E0E0E0') + 
  labs(title = 'Bias as a Function of Sample Size and Dimension of Theta', x = 'Sample Size') + 
  theme_user() + 
  theme(legend.position = 'bottom')

custom_ggsave('bias_function_n.png', p, width = 2, height = 1.2)
