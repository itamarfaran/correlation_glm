source('tex/simulations/aux_.R')

p <- 24
infect_regions <- 10
n <- 50
ARMA <- 0.5
n_sim <- 20


data_conf <- list(
  link = 'lib/data/Amnesia_all_AAL.mat',
  corr_matrix_name = 'corrmats',
  healthy_index_name = 'CONTROLS',
  sick_index_name = 'TGA'
)

subset <- sort(sample(86, size = p))
sample_data <- prepare_corrmat_data(
  link = data_conf$link,
  corr_matrix_name = data_conf$corr_matrix_name,
  healthy_index_name = data_conf$healthy_index_name,
  sick_index_name = data_conf$sick_index_name,
  subset = subset
)
corrmat_h <- calculate_mean_matrix(sample_data$samples$healthy)
corrmat_s_old <- calculate_mean_matrix(sample_data$samples$sick)

for(try_ in seq_len(1000)){
  infected_regions <- sort(sample(p, infect_regions))
  corrmat_s_new <- corrmat_h
  corrmat_s_new[infected_regions,infected_regions] <- corrmat_s_old[infected_regions,infected_regions]
  if(is.positive.definite(corrmat_s_new))
    # if(max(abs(corrmat_s_new - corrmat_h)) > 0.22)
    if(min(corrmat_s_new - corrmat_h) < -0.17)
      break()
}

corrplot(corrmat_s_new - corrmat_h, is.corr = F)

create_variance_estimates_custom_matrix <- function(n_sim, linkFun = linkFunctions$multiplicative_identity){
  if (ARMA == 0) ARMA <- NULL
  n_s <- ceiling(0.5*n)
  n_h <- n - n_s
  
  samples <- create_samples(n_sim = n_sim, n_h = n_h, n_s = n_s, p = p,
                            real_theta = corrmat_h, real_sick = corrmat_s_new, linkFun = linkFun,
                            ARsick = ARMA, ARhealth = ARMA, MAsick = ARMA, MAhealth = ARMA)
  results <- lapply(
    1:n_sim, function(i) estimate_alpha(
      healthy_dt = samples$samples[[i]]$healthy,
      sick_dt = samples$samples[[i]]$sick,
      linkFun = linkFun, verbose = FALSE)
  )
  
  gee_vars <- sapply(1:n_sim, function(i) compute_gee_variance(
    cov_obj = results[[i]],
    healthy_dt = samples$samples[[i]]$healthy,
    sick_dt = samples$samples[[i]]$sick
  ), simplify = 'array')
  
  sds <- lapply(1:n_sim, function(i) sqrt_diag(gee_vars[,,i]))
  
  
  out <- data.table(
    sims_num = rep(1:n_sim, each = p),
    real = as.vector(samples$alpha),
    estimate = as.vector(sapply(transpose(results)$alpha, as.vector)),
    sd = do.call(c, sds)
  )
  
  is_null <- rep(TRUE, nrow(out))
  is_null[infected_regions] <- FALSE
  out[,`:=`(
    z_value = (estimate - linkFun$NULL_VAL)/sd,
    is_null = is_null
    )]
  
  return(out)
}

res <- create_variance_estimates_custom_matrix(2)

res[,`:=`(
  p_value = 2*pnorm(abs(z_value), lower.tail = F),
  is_null_char = ifelse(is_null, 'Null Value', 'Non Null Value')
)]

inference_plot <- ggplot(res, aes(x = z_value, y = ..density.., fill = is_null_char)) + 
  geom_histogram(bins = 8*sqrt(length(res)), alpha = .5, position = 'identity', col = 'white') +
  # facet_grid(is_null_char~.) + 
  geom_hline(yintercept = 0) + 
  labs(
    x = 'Z scores \n of misspecified link function',
    y = 'Density',
    fill = ''
  ) +
  theme_user() + 
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    legend.position = c(.8,.75),
    legend.background = element_blank(),
    legend.title = element_blank(),
    legend.box.background = element_rect(colour = "black")
  ) + 
  scale_fill_manual(
    values = c('Null Value' = '#A6ACAF', 'Non Null Value' = '#212F3D')
  )


custom_ggsave('robustness_tga.png', inference_plot, width = 2, height = 1.1)
