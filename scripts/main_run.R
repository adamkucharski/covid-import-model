# Analysis of importations from India

library(tidyverse)
library(covidregionaldata)
library(forecast)
library(mgcv)
library(doMC)
library(mvtnorm)
library(data.table)


# Set directory
setwd("~/Documents/GitHub/covid-import-model")
data_path <- "~/Documents/COVID_data/B_617_2/"

registerDoMC(4)  #change to your number of CPU cores

# Load data ----------------------------------------------

source("R/load_data.R") # Check slow import function enable

# Load model functions ----------------------------------------------

source("R/model_functions.R")

source("R/model_mcmc.R")

# Run inference  -------------------------------------------------------------

# Define parameters
kk_pick <- 0.2 # overdispersion
cap_outbreak_size <- 1e4 # Cap simulations to prevent runaway epidemics
priorScale <- function(x){ifelse(abs(x)<=1,1,0)}

# Fit model
run_transmission_mcmc(MCMC.runs = 1e5)

# Simulate & plot outputs  -------------------------------------------------------------

# Extract posteriors
for( iiM in 2){
  source("R/load_posterior.R",local=TRUE)
  
  # Run fitted model
  theta_mle <- data.frame(thetatab[pick.max,])
  #theta_mle <- c(rr=1.6,r_scale=0.9,r_scale_2=0.9,decline=0.2,dt_decline=0.5,imp=1.5)
  output1 <- fit_R_deterministic(theta_mle,add_days = 0)
  
  source("R/plot_outputs.R")

}

# plot_post()

# compare_R_fits()


