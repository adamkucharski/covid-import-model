# Analysis of importations from India

library(tidyverse)
library(covidregionaldata)
library(forecast)
library(mgcv)
library(doMC)
library(mvtnorm)
library(data.table)


# Set directories
setwd("~/Documents/GitHub/covid-import-model")
data_path <- "~/Documents/COVID_data/B_617_2/"

registerDoMC(4)  #change to your number of CPU cores

# Load data ----------------------------------------------
local_run <- F
location_ID <- "UK"
local_nn <- 1

source("R/load_data.R") # Load and format remaining data

# Load model functions ----------------------------------------------

source("R/model_functions.R")

source("R/model_mcmc.R")

# Run inference  -------------------------------------------------------------

# Define parameters
cap_outbreak_size <- 1e4 # Cap simulations to prevent runaway epidemics
priorScale <- function(x){ifelse(abs(x)<=1,1,0)} # Prior on relative R values for non-travellers
priorScale2 <- function(x){1} #ifelse(abs(x)<=1,1,0)} # Prior on relative R values for non-travellers
priorRep <- function(x){ifelse(x>=1,1,0)} #ifelse(abs(x)<=1,1,0)} # Prior on reporting value
priorTime <- function(x){ifelse(x>= (voc_n) & x < (voc_n+7),1,0)} # Prior on timing of surge effect
priorTime2 <- function(x){ifelse(x>= (step_3_n) & x < (step_3_n+7),1,0)} # Prior on timing of end change in R

# Fit model
run_transmission_mcmc(MCMC.runs = 1e5) # Specify number of MCMC iterations: >1e5 recommended

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


