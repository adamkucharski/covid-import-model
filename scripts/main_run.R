# Analysis of importations from India

library(tidyverse)
library(covidregionaldata)
library(forecast)
library(mgcv)
library(doMC)


# Set directory
setwd("~/Documents/GitHub/covid-import-model")
data_path <- "~/Dropbox/LSHTM/2020_nCoV_main_db/B117_variant/relative_advantage/"

registerDoMC(4)  #change to your number of CPU cores

# Set up data ----------------------------------------------

source("R/load_data.R") # Check slow import function enable

source("R/model_functions.R")

# Run inference  -------------------------------------------------------------

# Define parameters
kk_pick <- 0.2
daily_decline <- 0.01 # assumed daily growth rate
under_factor <- 1
run_n <- 10
rr_range <- seq(1.4,1.8,0.05)
decline_range <- seq(0.02,0.04,0.005)
imp_range <- seq(0.6,2,0.2)

# Set up grid search
parameter_list <- expand.grid(rr_range, decline_range,imp_range)
parameter_list <- data.frame(parameter_list); names(parameter_list) <- c("rr","decline","imp")

source("R/format_data_specific.R")

# Fit model
run_inference()

# Extract MLE
output_mle <- get_MLE()
output1 <- output_mle$out

# Plot outputs  -------------------------------------------------------------

source("R/plot_outputs.R")


