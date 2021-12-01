load('data/adaptMCMC_final.RData')

# prepare mcmc results --------------------
n = nrow(samp[[1]]$samples) # number of iterations for each chain

results = rbind.data.frame(
  samp[[1]]$samples,
  samp[[2]]$samples,
  samp[[3]]$samples,
  samp[[4]]$samples) %>%
  cbind(chain = c(rep(c(1:4), each=n)))
names(results) = c(theta_names, 'chain')

res_long = results %>%
  group_by(chain) %>%
  mutate(i = 1:n()) %>%
  pivot_longer(-c(chain, i)) %>%
  mutate(name = factor(name, levels=names(results)))

burn_in = 5E4
thin = 10

# Posterior mean by chain ------------------
res_long %>%
  filter(i %in% seq(burn_in, n, thin)) %>%
  group_by(chain, name) %>%
  summarise(m = mean(value)) %>%
  pivot_wider(names_from = name, values_from = m) %>%
  mutate(across(alpha:y_init, ~exp(.x)))

# Summary -------------------
res_long %>%
  filter(i %in% seq(burn_in, n, thin)) %>%
  group_by(name) %>%
  summarise(m = mean(value),
            l = hdi(value)[1],
            u = hdi(value)[2]) %>%
  mutate(across(m:u, ~exp(.x)))

# Traceplots -----------------------
p1 = res_long %>%
  filter(i %in% seq(burn_in, n, thin)) %>%
  ggplot(aes(x=i, y=value, col=factor(chain))) + 
  geom_line(alpha = 0.8) +
  facet_wrap(~name, scales='free_y') +
  scale_x_continuous(labels = scales::label_number_si()) +
  theme_bw() + 
  labs(x = 'iteration', col='chain') +
  theme(legend.position='bottom',
        legend.margin=margin(0,0,0,0),
        legend.box.margin=margin(-10,-10,0,-10)) + 
  guides(color = guide_legend(override.aes = list(alpha = 1)))

ggsave('output/traceplot.pdf', plot = p1, device='pdf')

# Posterior density -------------------------
p2 = res_long %>%
  filter(i %in% seq(burn_in, n, thin)) %>%
  mutate(value = exp(value)) %>%
  ggplot() + 
  geom_histogram(aes(x=value,fill=factor(chain)), 
                 position='identity', 
                 alpha = 0.5, bins=50) + 
  facet_wrap(~name, scales='free') +
  labs(y = 'frequency', fill = 'chain') +
  theme_bw() + 
  theme(legend.position='bottom',
        legend.margin=margin(0,0,0,0),
        legend.box.margin=margin(-10,-10,0,-10))

p2
ggsave('output/postdens.pdf', plot = p2, device='pdf')

# Posterior density with prior (TODO) -------------------------
ag_prior = function(x) dnorm(x, 1, 0.5)
bd_prior = function(x) dnorm(x, 0.05, 0.05)
x0_prior = function(x) dlnorm(x, log(19.58), 0.05)
y0_prior = function(x) dlnorm(x, log(30.09), 0.05)

xs = c(seq(0, 2, length.out = 1E3), seq(10, 60, length.out = 1E3))

data.frame(x = xs,
           a_t = sapply(xs, ag_prior),
           b_t = sapply(xs, bd_prior),
           g_t = sapply(xs, ag_prior),
           d_t = sapply(xs, bd_prior),
           x_t = sapply(xs, x0_prior),
           y_t = sapply(xs, y0_prior)) %>%
  pivot_longer(-x) %>%
  ggplot(aes(x=x, y=value)) + 
  geom_line() + 
  facet_wrap(~name)

plot(xs, sapply(xs, bd_prior))

# Trajectory Plot ----------------------
theta_fit = log(c(0.396,0.0151,1.14,0.0334,55.0,29.9))
z_fit = sim_z(theta_fit, times_full)
fit_plt(z_fit, lynx_hare_df)

times_full = lynx_hare_df$year
times_pred = 1913:1935

z_fit$test = ifelse(z_fit$time %in% times_pred, 1, 0)

p3 = left_join(lynx_hare_df, z_fit, by=c('year'='time'),
          suffix = c('_dat','_mod')) %>%
  pivot_longer(-c(year,test),
               names_to = c('species','type'),
               names_sep = '_',
               values_to = 'population') %>%
  pivot_wider(names_from = type,
              values_from = population) %>%
  ggplot(aes(x = year, col = species)) + 
  geom_point(aes(y = dat)) +
  geom_line(aes(y = mod, lty=factor(test))) + 
  theme_bw() + 
  theme(legend.position='bottom',
        plot.margin=grid::unit(c(0.5,0.5,0,0), "mm"),
        legend.margin=margin(0,0,0,0),
        legend.box.margin=margin(-10,-10,0,-10)) + 
  labs(y = 'pelts collected (in thousands)') +
  scale_linetype_discrete('type', labels=c('train','prediction')) +
  scale_color_discrete(labels=c('lynx','hare')) 
p3
ggsave('output/traj.pdf', plot = p3, device='pdf')
