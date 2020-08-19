source('tex/simulations/aux_.R')

p <- 20
df <- 4000
n <- 50000

cov_mat <- matrix(runif(2*p^2, -1, 1), 2*p, p)
cov_mat <- t(cov_mat) %*% cov_mat / (2*p)
cor_mat <- force_symmetry(cov2cor(cov_mat))

sample <- rWishart(n, df, cov_mat)/df
sample <- sapply(1:n, function(i) force_symmetry(cov2cor(sample[,,i])), simplify = 'array')
sample <- convert_corr_array_to_data_matrix(sample)

cov_emp <- cov(sample)
cov_theo <- corrmat_covariance(cor_mat)/df

plot(
  diag(cov_theo),
  (1 - triangle2vector(cor_mat)^2)^2/df
  
)


fisher_z <- function(r) 0.5*log((1 + r)/(1 - r))
cov_emp_z <- cov(fisher_z(sample))
diags_ <- cor_mat %>% triangle2vector() %>% (function(r) 1/(1 - r^2)) %>% diag()
cov_theo_z <- diags_ %*% cov_theo %*% diags_

toplot <- data.table(theo = triangle2vector(cov_theo_z), emp = triangle2vector(cov_emp_z))
toplot <- data.table(theo = triangle2vector(cov_theo), emp = triangle2vector(cov_emp))
toplot[, diff := theo - emp]

toplot_samp <- toplot[sample(.N, 1250, prob = sqrt(abs(emp)))]

p1 <- ggplot(toplot, aes(x = theo, y = emp)) + 
  geom_point(shape=21, alpha=0.5, fill='darkgrey') + 
  geom_abline(intercept = 0, slope = 1) + 
  labs(
    # title='Empirical and Theoretical Covariance of the Correlation Matrix',
    x='Delta Method Derived Estimates',
    y='Parametric Bootstrap Estimates'
    ) + theme_user()
plot(p1)
