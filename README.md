# covid-import-model

Model of variant importations and onwards community transmission. Currently implementation uses a grid search and profile likelihood approach to estimate paramaters using a continuous time overdispersed branching process model.

### Quick start guide

Data loading and main model run script (including file paths) are in `scripts/main_run.r`. This calls the following R files:

> `R/load_data.R` - Load data from `covidregionaldata`, COG-UK, PHE imports and India sequences (stored elsewhere).

> `R/format_data_specific.R` - Format data to line up dates and convert in moving averages ec.

> `R/model_functions.R` - Main model runs and likelihood for MLE fitting.

> `R/plot_outputs.R` - Plot data, parameter estimates and model simulations using MLE fits.
