# Load packages
library(tidyverse)
library(deSolve)
library(adaptMCMC)
library(HDInterval)

source('code/01-data.R') # Load data
source('code/02-mcmc-setup.R') # set up mcmc
source('code/03-mcmc-run.R') # run mcmc
source('code/04-analysis.R') # analyse mcmc results