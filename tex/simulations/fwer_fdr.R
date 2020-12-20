source('tex/simulations/aux_.R')

# todo: check on strung null with FWER = P(ANY(alpha_null  <= 1))

n_sim <- 60
n <- 100
p <- 32
range_alpha <- c(0.7, 0.9)
ARMA <- 0.6
reps <- 20
sig_level <- .05

samples <- rbind(
  rbindlist(lapply(
    1:reps, create_sample_estimates,
    n_sim = n_sim, n = n, p = p, p_s = 0.5, percent_alpha = 0, range_alpha = range_alpha, ARMA = 0,
    ncores = ncores
  )),
  rbindlist(lapply(
    1:reps, create_sample_estimates,
    n_sim = n_sim, n = n, p = p, p_s = 0.5, percent_alpha = 0, range_alpha = range_alpha, ARMA = 0.6,
    ncores = ncores
  )),
  rbindlist(lapply(
    1:reps, create_sample_estimates,
    n_sim = n_sim, n = n, p = p, p_s = 0.5, percent_alpha = 0.5, range_alpha = range_alpha, ARMA = 0,
    ncores = ncores
  )),
  rbindlist(lapply(
    1:reps, create_sample_estimates,
    n_sim = n_sim, n = n, p = p, p_s = 0.5, percent_alpha = 0.5, range_alpha = range_alpha, ARMA = 0.6,
    ncores = ncores
  ))
)

samples[,z_value := (estimated_alpha - 1)/sd]
samples[,p_value := 2*pnorm(abs(z_value), lower.tail = F)]
samples[,`:=`(
  holm_p = p.adjust(p_value, 'bonferroni'),
  bh_p = p.adjust(p_value, 'BH')
), by = .(autocorrelated, p_s, case, sim_num, sim)]

out <- rbind(
  samples %>%
    filter(real_alpha == 1) %>%
    group_by(autocorrelated, p_s, case, sim, sim_num) %>%
    summarise(value = mean(p_value < sig_level)) %>% 
    group_by(autocorrelated, p_s, case, sim) %>% 
    summarise(type = 'No Correction\n(False Positive Rate\nper Sample)', value = mean(value)),

  samples %>%
    filter(real_alpha == 1) %>%
    group_by(autocorrelated, p_s, case, sim, sim_num) %>%
    summarise(value = any(holm_p < sig_level)) %>% 
    group_by(autocorrelated, p_s, case, sim) %>%
    summarise(type = 'FWER\n(Bonferroni)', value = mean(value)),
  
  samples %>% 
    group_by(autocorrelated, p_s, case, sim, sim_num) %>%
    summarise(V = sum(bh_p < sig_level & real_alpha == 1),
              R = sum(bh_p < sig_level)) %>% 
    mutate(value = V/remove_zeros(R)) %>%
    group_by(autocorrelated, p_s, case, sim) %>% 
    summarise(type = 'FDR\n(BH)', value = mean(value))
  )


save(samples, out, file = 'tex/simulations/fwer_fdr.RData')

pl <- ggplot(out, aes(x = type, y = value)) + 
  geom_boxplot(fill = 'lightgrey') +
  facet_grid(case~autocorrelated ) + 
  geom_hline(yintercept = 0) + 
  geom_hline(yintercept = sig_level, linetype = 3) + 
  ylim(0, 0.2) + theme_user() + 
  labs(title = 'FDR & FWER', subtitle=paste0('Error Rate set to ', sig_level), y = 'Rate') + 
  theme(axis.title.x=element_blank())

custom_ggsave('fwer_fdr.png', pl)
