load('tex/tga_analysis/analysis_tga_multiplicative_identitynobiascorr.RData')
orig <- list(sample_data=sample_data,
             results=results,
             gee_var=gee_var)

load('tex/tga_analysis/analysis_tga_log_multiplicative_identity_log.RData')
log_transform <- list(sample_data=sample_data,
                      results=results,
                      gee_var=gee_var)

load('tex/tga_analysis/analysis_tga_multiplicative_identity_rev.RData')
reverse_samples <- list(sample_data=sample_data,
                        results=results,
                        gee_var=gee_var)

index <- seq_along(orig$results$alpha)
toplt_log <- rbind(
  tibble(link='identity', index=index, values=as.vector(orig$results$alpha), null=1),
  tibble(link='log', index=index, values=as.vector(log_transform$results$alpha), null=0)
)

toplt_long %>% 
  filter(link %in% c('identity', 'log')) %>%
  ggplot(aes(x = index, y = values)) + 
  geom_line(aes(y=null)) + 
  geom_point() + 
  facet_grid(link ~ ., scales = 'free_y')

toplt_rev <- rbind(
  tibble(type='estimate', index=index, orig=as.vector(orig$results$alpha), rev=as.vector(reverse_samples$results$alpha)),
  tibble(type='sd', index=index, orig=sqrt(diag(orig$gee_var)), rev=sqrt(diag(reverse_samples$gee_var)))
)




toplt_long %>% 
  filter(link %in% c('identity', 'reverse')) %>%
  ggplot(aes(x = index, y = values)) + 
  geom_line(aes(y=null)) + 
  geom_point() + 
  facet_grid(link ~ ., scales = 'free_y')

toplt_rev %>%
  filter(type=='estimate') %>%
  mutate(orig = orig - 1,
         rev = rev - 1) %>%
  ggplot(aes(x=orig, y=rev)) + 
  geom_point() + 
  geom_abline(slope=-1, intercept = 0)

toplt_rev %>%
  filter(type=='estimate') %>%
  mutate(orig = orig - 1,
         rev = rev - 1) %>%
  lm(rev ~ orig, .) %>%
  summary()

toplt_rev %>% 
  gather('link', 'value', -index, -type) %>%
  spread(type, value) %>% 
  mutate(z_val = (estimate - 1)/sd,
         p_val = 2*pnorm(abs(z_val), lower.tail = F)) %>%
  # group_by(link) %>% 
  # mutate(p_adj = p.adjust(p_val, 'BH')) %>%
  # ungroup() %>%
  select(index, link, z_val, p_val) %>%
  mutate(z_val = abs(z_val)) %>%
  gather('type', 'value', -index, -link) %>% 
  spread(link, value) %>% 
  ggplot(aes(x = orig, y=rev)) + 
  geom_point() + 
  facet_wrap(vars(type), scales = 'free')
