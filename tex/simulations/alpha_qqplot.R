source('tex/simulations/aux_.R')

n_sim <- 100
p <- 24
percent_alpha <- 0
range_alpha <- c(1,1)


data <- create_samples(n_sim = n_sim, n_h = 25, n_s = 25, p = p, Tlength = 115,
               percent_alpha = 0, range_alpha = c(1,1), ncores = ncores)

estimates <- pbmclapply(data$samples, function(lst) estimate_alpha(lst$healthy, lst$sick, verbose=FALSE)$alpha, mc.cores = ncores)

toplot <- data.table(sim = rep(seq_len(n_sim), each = p),
                     index = rep(seq_len(p), times = n_sim),
                     estimate = do.call(c, lapply(estimates, as.vector)))

toplot %>% 
  ggplot(aes(sample=estimate)) +
  stat_qq(size = 1) +
  stat_qq_line() + 
  facet_wrap(vars(index))
