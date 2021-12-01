# Set up MCMC
theta_names = c('alpha', 'beta', 'gamma', 'delta', 'x_init', 'y_init')
theta = log(c(0.5, 0.02, 0.9, 0.02, 19.58, 30.09))

# differential equation for deSolve
dz_dt <- function (t, z, theta) {
  theta = lapply(theta, exp)
  with(as.list(c(z, theta)), {
    dx = prey*(alpha - beta*pred)
    dy = -pred*(gamma - delta*prey)
    return(list(c(dx, dy)))
  })
}

# simulate one run given log parameters
sim_z <- function(logtheta, times){
  params = as.list(logtheta[1:4])
  names(params) = theta_names[1:4]
  z_init = as.numeric(logtheta[5:6])
  names(z_init) = c('prey','pred')
  mod <- as.data.frame(
    ode(func = dz_dt, 
        y = exp(z_init), 
        parms = params, 
        times = times,
        method="ode45"))
  return(mod)
}

# plot results
fit_plt <- function(mod, data){
  left_join(data, mod, by=c('year'='time'),
            suffix = c('_dat','_mod')) %>%
    pivot_longer(-year,
                 names_to = c('species','type'),
                 names_sep = '_',
                 values_to = 'population') %>%
    pivot_wider(names_from = type,
                values_from = population) %>%
    ggplot(aes(x = year, col = species)) + 
    geom_line(aes(y = mod)) + 
    geom_point(aes(y = dat)) +
    theme_bw() + 
    theme(legend.position='bottom',
          plot.margin=grid::unit(c(0,0,0,0), "mm"),
          legend.margin=margin(0,0,0,0),
          legend.box.margin=margin(-10,-10,0,-10)) + 
    labs(y = 'population (in thousands)') +
    scale_color_discrete(labels=c('lynx','hare'))
}

# log prior 
log_prior = function(theta){
  theta = exp(theta)
  sum(dnorm(c(theta[1], theta[3]), 1, 0.25, log=T), 
      dnorm(c(theta[2], theta[4]), 0.05, 0.05, log=T),
      dlnorm(c(theta[5:6]), log(c(19.58, 30.09)), 0.05, log=T))
}

# log likelihood
log_lik = function(theta){
  # Run model
  mod = sim_z(theta)
  
  # calculate loglik
  sum(dlnorm(mod$prey, log(data$prey), 0.01, log=T),
      dlnorm(mod$pred, log(data$pred), 0.01, log=T))
}

# combined log density
log_posterior = function(theta){log_lik(theta) + log_prior(theta)}

# metropolis step for regular MCMC (not used)
metropolis_step <- function(theta, log_posterior, proposal) {
  propTheta <- proposal(theta)
  a <- log_posterior(propTheta) - log_posterior(theta)
  if (a > 0){
    propTheta
  } else if (log(runif(1)) < a){
    propTheta
  } else {
    theta
  }
}

metropolis <- function(theta, log_posterior, proposal, m) {
  out = matrix(NA_real_, nrow = m, ncol = length(theta))
  out[1, ] = theta
  for (i in 2:m) {
    out[i, ] <- metropolis_step(out[i-1, ], log_posterior, proposal)
  }
  out
}

proposal <- function(x) {
  z = runif(length(x), -1E-3, 1E-3)
  x + z
}