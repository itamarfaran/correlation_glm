source('main_work/simulations/auxilary_functions.R')

timeit <- function(p, r){
  raw <- function(p){
    samples <- create_samples(n_sim = 1, n_h = 50, n_s = 50, p = p, Tlength = 115,
                              percent_alpha = 0, range_alpha = c(1,1),
                              ARsick = NULL, ARhealth = NULL, MAsick = NULL, MAhealth = NULL)
    
    t1 <- Sys.time()
    results <- estimate_alpha(healthy_dt = samples$samples$healthy, sick_dt = samples$samples$sick, verbose = FALSE)
    t1 <- as.double(Sys.time() - t1, units = 'secs')

    t2 <- Sys.time()
    compute_gee_variance(results, samples$samples$healthy, samples$samples$sick)
    t2 <- as.double(Sys.time() - t2, units = 'secs')
    
    return(c(t1, t2))
  }
  
  cat(paste0(p, ', '))
  out <- matrix(0, nr = r, nc = 2)
  for(i in 1:r) out[i,] <- raw(p)
  out <- data.table(out)
  colnames(out) <- c('est', 'gee')
  out[, p := p]
  return(out)
}

p <- 5*(2:12)
r <- 10
out <- data.table(p = numeric(), est = numeric(), gee = numeric())

for(p_ in p) out <- rbind(out, timeit(p_, r))
out[,`:=`(Estimation = est, GEE = gee, Total = est + gee, est = NULL, gee = NULL)]

cols = colnames(out)[-1]

results <- out[,.(
  Type = cols,
  mean = sapply(.SD, mean),
  sd = sapply(.SD, sd),
  median = sapply(.SD, median),
  max = sapply(.SD, max)
  ), .SDcols = cols, by = p]

save(out, results, file = 'main_work/simulations/performance_data.RData')

# results[,.(p, Type, Median = median, Max = max)] %>% 
#   melt(id.vars = c('p', 'Type'), variable.name = 'Measure', value.name = 'Value') %>%
p <- results %>%
  ggplot(aes(x = p, y = mean, linetype = Type, size = Type)) +
  geom_line() +
  scale_color_manual(values = c('Median' = 'black', 'Max' = 'grey')) + 
  scale_linetype_manual(values = c('Estimation' = 5, 'GEE' = 6, 'Total' = 1)) + 
  scale_size_manual(values = c('Estimation' = 0.8, 'GEE' = 0.8, 'Total' = 1.2)) +
  labs(
    title = 'Computation Time',
    x = 'P',
    y = 'Average Time (in Seconds)'
    ) + 
  theme_user() + 
  theme(legend.position = 'bottom')

custom_ggsave('performance_sim.png', p)

