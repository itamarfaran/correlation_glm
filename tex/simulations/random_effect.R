source('tex/simulations/aux.R')

corr_mat <- force_symmetry(build_parameters(12, 0, 0:1)$corr_mat)

arma_vec <- c(-.8, -.5, -.2, 0, .2, .5, .8)
rand_vec <- c(-1, 10, 50, 100, 200, 500, 1000, Inf)
examples <- data.table(expand.grid(
  df = 12*(1:10),
  sample_size = 20*(1:6),
  AR = arma_vec,
  MA = arma_vec,
  random_effect = rand_vec
))
simulate_rms <- Vectorize(function(df, sample_size, AR, MA, random_effect, ncores = 1){
  create_correlation_matrices(
    real_corr = corr_mat,
    df = 100,
    sample_size = 40,
    AR = 0.5,
    MA = -0.2,
    random_effect = -1,
    ncores = ncores) %>% 
    calculate_mean_matrix() %>% 
    (function(x) {diag(x) <- 1; x}) %>% 
    efrons_rms() %>%
    return()
})

if(ncores > 1){
  tt <- function(i) with(examples[i], simulate_rms(
    sample_size = sample_size,
    AR = AR, MA = MA,
    random_effect = random_effect))
  res <- simplify(pbmclapply(1:nrow(examples), tt, mc.cores = ncores))
  examples[,rms := res]
} else {
  examples[,rms := simulate_rms(
    sample_size = sample_size,
    AR = AR, MA = MA,
    random_effect = random_effect
  )]
}
examples[random_effect < 0, random_effect := 0]
save(examples, file='main_work/simulations/random_effect.RData')

plot_fun <- function(col){
  p <- ggplot(examples[sample(.N, 2500)], aes_string(x = col, y = 'rms')) + 
    geom_point(fill = 'grey', position = 'jitter', alpha = .4) +
    ylim(examples[,range(rms)]*c(0.97, 1.03)) +
    labs(x='', y='', title = col) + 
    theme_user()
  if(col == 'random_effect') p <- p + scale_x_sqrt()
  return(p)
}

p <- arrangeGrob(grobs = lapply(colnames(examples), plot_fun), nrow=2, ncol=3)
custom_ggsave('random_effect.png', p)
plot(p)
