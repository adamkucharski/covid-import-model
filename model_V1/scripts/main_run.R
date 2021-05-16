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
run_n <- 20
rr_range <- seq(2,4,0.5)
decline_range <- seq(0.02,0.05,0.01)
imp_range <- seq(0.2,1,0.2) # seq(0.9,1.5,0.1)

# Set up grid search
parameter_list <- expand.grid(rr_range, decline_range,imp_range)
parameter_list <- data.frame(parameter_list); names(parameter_list) <- c("rr","decline","imp")

source("R/format_data_specific.R")

# Fit model
run_inference()

# Extract MLE
output_mle <- get_MLE()

# Simulate & plot outputs  -------------------------------------------------------------

# Run fitted model
add_d <- 10
mle_val <- output_mle$mle
mle_val$rr <- 1.6; mle_val$decline <- 0.03; mle_val$imp <-  1

output1 <- fit_R(mle_val$rr,run_n,add_days = add_d, import_f = mle_val$imp, daily_decline=mle_val$decline)

source("R/plot_outputs.R")


