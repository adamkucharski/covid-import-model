# Set up model-specific data and parameters ----------------------------------------------

# Imports in PHE data Feb,Mar,Apr
south_asia_imports <- 1079
total_cases_india <- all_india %>% filter(country == "India", date<as.Date("2021-04-23")) %>% select(cases_new) %>% sum()
travel_multiplier <- south_asia_imports/total_cases_india

all_india <- all_india %>% mutate(daily_imports = travel_multiplier*cases_new)

# Downweight recent imports based on incubation period (from McAloon et al. 2020)
red_list_point <- as.numeric(india_red_list-date_pick)
total_days <- length(all_india$date)
downweight_imports <- c(rep(1,red_list_point),1-plnorm(1:(total_days-red_list_point),mean=log(5.1),sd=log(1.65)))

# Format data
ma_India_variant0 <- ma(data_india$B.1.617.2,7) # Moving average of UK cases
l_IV <- length(ma_India_variant0)
final_value <- 1 #tail(ma_India_variant0,4)[1]
ma_India_variant <- c(ma_India_variant0[1:(l_IV-4)],rep(final_value,total_days-l_IV+4))
ma_India_variant[1:3] <- 0

# Estimated ndia daily B.1.617.2 imports
daily_india <- under_factor*round(downweight_imports*all_india$daily_imports*ma_India_variant)

# Imported cases
daily_india_seq <- 0*daily_india # Imports based on traveller cases
daily_india_seq[match(traveller_cases0$Date,all_india$date)] <- traveller_cases0$Number_B_1_617_2

# Moving average and scaling factors
ma_UK_cases <- ma(all_uk$cases_new,7) # Moving average of UK cases

data_fit <- data_proportion[match(all_uk$date,data_proportion$sample_date),]
t_fit <- nrow(data_fit)