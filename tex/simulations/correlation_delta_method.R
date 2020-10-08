source('tex/simulations/aux_.R')

p <- 20
df <- 40
n <- 500

cov_mat <- matrix(runif(2*p^2, -1, 1), 2*p, p)
cov_mat <- t(cov_mat) %*% cov_mat / (2*p)
cor_mat <- force_symmetry(cov2cor(cov_mat))

sample <- rWishart(n, df, cov_mat)/df
sample <- sapply(1:n, function(i) force_symmetry(cov2cor(sample[,,i])), simplify = 'array')
sample <- convert_corr_array_to_data_matrix(sample)

cov_emp <- cov(sample)
cov_theo <- corrmat_covariance(cor_mat)/df

toplot <- data.table(theo = triangle2vector(cov_theo), emp = triangle2vector(cov_emp))
toplot[, diff := theo - emp]

toplot_samp <- toplot[sample(.N, 1250, prob = sqrt(abs(emp)))]

p1 <- ggplot(toplot_samp, aes(x = theo, y = emp)) + 
  geom_point(shape=21, alpha=0.5, fill='darkgrey') + 
  geom_abline(intercept = 0, slope = 1) + 
  labs(
    # title='Empirical and Theoretical Covariance of the Correlation Matrix',
    x='Delta Method Derived Estimates',
    y='Parametric Bootstrap Estimates'
    ) + theme_user()

# p2 <- ggplot(toplot_samp, aes(x = theo, y = diff)) + 
#   geom_point(shape=21, alpha=0.5, fill='darkgrey') + 
#   geom_hline(yintercept = 0) + 
#   labs(
#     # title='Empirical and Theoretical Covariance of the Correlation Matrix',
#     x='Delta Method Derived Estimates',
#     y='Difference Between Bootstrap & Delta Estimates'
#   ) + theme_user()
# 
# out <- arrangeGrob(p1, p2, ncol=1)
out <- p1
custom_ggsave('correlation_delta_method.png', out, width = .75, height = .75)
