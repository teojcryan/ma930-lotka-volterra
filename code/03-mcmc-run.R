# Run MCMC
n_iter = 4E5

samp <- MCMC.parallel(log_posterior, n=n_iter, init=theta,
                      adapt=TRUE, acc.rate=0.234,
                      n.chain = 4, n.cpu = 8,
                      packages=c('deSolve'))