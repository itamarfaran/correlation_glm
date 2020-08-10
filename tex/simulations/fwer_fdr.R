source('tex/simulations/aux.R')

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
  holm_p = p.adjust(p_value, 'holm'),
  bh_p = p.adjust(p_value, 'BH')
), by = .(autocorrelated, p_s, case, sim_num, sim)]

out <- rbind(
  samples[real_alpha == 1,.(value = sum(p_value < sig_level)/.N),
          by = .(autocorrelated, p_s, case, sim_num, sim)][
            ,.(type = 'No Correction\n(False Positive Rate\nper Sample)', value = mean(value))
            , by = .(autocorrelated, p_s, case, sim)],
  
  samples[real_alpha == 1,.(value = sum(holm_p < sig_level)/.N),
          by = .(autocorrelated, p_s, case, sim_num, sim)][
            ,.(type = 'FWER\n(Holm)', value = mean(value))
            , by = .(autocorrelated, p_s, case, sim)],
  
  samples[,.(value = sum(bh_p < sig_level & real_alpha == 1)/sum(bh_p < sig_level)),
          by = .(autocorrelated, p_s, case, sim_num, sim)][
            ,.(type = 'FDR\n(BH)', value = mean(value))
            , by = .(autocorrelated, p_s, case, sim)]
)

save(samples, out, file = 'main_work/simulations/fwer_fdr.RData')

pl <- ggplot(out, aes(x = type, y = value)) + 
  geom_boxplot(fill = 'lightgrey') +
  facet_grid(case~autocorrelated ) + 
  geom_hline(yintercept = 0) + 
  geom_hline(yintercept = sig_level, linetype = 3) + 
  ylim(0, 0.2) + theme_user() + 
  labs(title = 'FDR & FWER', subtitle=paste0('Error Rate set to ', sig_level), y = 'Rate') + 
  theme(axis.title.x=element_blank())

custom_ggsave('fwer_fdr.png', pl)
