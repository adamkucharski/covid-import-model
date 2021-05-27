# Analysis of importations from India

library(tidyverse)
library(covidregionaldata)
library(forecast)
library(mgcv)
library(doMC)
library(mvtnorm)
library(data.table)


# Set directories
# setwd("~/Documents/GitHub/covid-import-model")
# data_path <- "~/Documents/COVID_data/B_617_2/"

registerDoMC(4)  #change to your number of CPU cores

# Load data ----------------------------------------------

source("R/load_data.R") # Load and format remaining data

# Load model functions ----------------------------------------------

source("R/model_functions.R")

source("R/model_mcmc.R")

# Run inference  -------------------------------------------------------------

# Define parameters
cap_outbreak_size <- 1e4 # Cap simulations to prevent runaway epidemics
priorScale <- function(x){ifelse(abs(x)<=1,1,0)} # Prior on relative R values for non-travellers
priorTime <- function(x){ifelse(x< total_days,1,0)} # Prior on timing of surge effect

# Fit model
run_transmission_mcmc(MCMC.runs = 1e4) # Specify number of MCMC iterations: >1e5 recommended

# Simulate & plot outputs  -------------------------------------------------------------

# Extract posteriors
for( iiM in 2){ # Fit two levels (i.e. traveller/non-traveller)
  source("R/load_posterior.R",local=TRUE)
  
  # Run fitted model
  theta_mle <- data.frame(thetatab[pick.max,]) # Extract max likelihood
  output1 <- fit_R_deterministic(theta_mle,add_days = 0) # Simulate epidemic
  
  source("R/plot_outputs.R") # Plot main figure

}

# Plot posteriors
# plot_post()


