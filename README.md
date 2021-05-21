# covid-import-model

Preliminary modelling analysis of variant importations and onwards community transmission. Currently implementation uses adaptive Metropolis Hastings MCMC to estimate paramaters using a deterministic approximation of a continuous time branching process model.

_Note: this is working repository, so code and data are likely to change over time_

### Quick start guide

1. Get required data files from COG-UK and India from GISAID via outbreak.data and place in a local folder:

	 * Download [cog_metadata.csv](https://cog-uk.s3.climb.ac.uk/phylogenetics/latest/cog_metadata.csv)

	 * Download `outbreakinfo_mutation_report_data_[date].tsv` from the lineage prevalence visualisation [here](https://outbreak.info/location-reports?loc=IND) and rename as `outbreakinfo_mutation_report_data.tsv`

2. Define file paths:

	 * Set GitHub directory, e.g. `setwd("~/Documents/GitHub/covid-import-model")`

	 * Specify local data path with above files stored, e.g. `data_path <- "~/Documents/COVID_data/B_617_2/"`

3. Main model script in `scripts/main_run.r`, including dependencies and data loading. This calls the following R files:

	> `R/load_data.R` - Load data from `covidregionaldata`, COG-UK and India sequences (stored in folder defined above).

	> `R/model_functions.R` - Main model simulation functions.

	> `R/model_mcmc.R` - MCMC inference functions.

	> `R/load_posteriors.R` - Load and format MCMC posteriors.

	> `R/plot_outputs.R` - Plot data, parameter estimates and model simulations.


### Archived code

An MLE framework used to implement an earlier version of this analysis is archived in `V1_code/`.

### Citation

Reference for initial version of the analysis: [report 1 citation to come]. If you plan to build on or cite this preliminary analysis for an academic publication, please ensure that you credit the underlying data sources above appropriately.
