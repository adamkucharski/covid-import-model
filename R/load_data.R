
# - - -
# Load covidregionaldata

if(!exists("all_countries")){
  all_india <- get_national_data("india") # Load data from covidregionaldata (slow)
  all_uk <- get_national_data("uk",source="who") # Load data from covidregionaldata (slow)
  
  # Fix missing data in WHO file
  missing_data <- c(2027,1926,1979)
  all_uk[all_uk$date==as.Date("2021-05-16"),"cases_new"] <- missing_data[1]
  all_uk[all_uk$date==as.Date("2021-05-17"),"cases_new"] <- missing_data[2]
  all_uk[all_uk$date==as.Date("2021-05-18"),"cases_new"] <- missing_data[3]
}

# - - -
# Load downloaded data

if(!exists("data_proportion")){
    
  cog0 = fread(paste0(data_path,"cog_metadata.csv"))[country == "UK" & pillar_2 == TRUE]
  cog = cog0[, .(.N,
                 B.1.617.1 = sum(lineage %like% "B\\.1\\.617\\.1"),
                 B.1.617.2 = sum(lineage %like% "B\\.1\\.617\\.2"),
                 B.1.617.3 = sum(lineage %like% "B\\.1\\.617\\.3")), keyby = sample_date]
  fwrite(cog, paste0(data_path,"COG_UK_out.csv"))

}

data_proportion <- read_csv(paste0(data_path,"COG_UK_out.csv"))

data_india <- read_tsv("data/outbreakinfo_mutation_report_data.tsv")


traveller_cases0 <- read_csv("data/voc_imports_2021_05_13.csv") # Read in from repo

# - - -
# Load variant data

india_red_list <- as.Date("2021-04-23")
date_pick <- as.Date("2021-02-01")
date_uk_fit <- as.Date("2021-04-23")

all_india <- all_india[1:nrow(all_uk),] # Avoid mismatched lengths
all_india <- all_india %>% filter(date>date_pick)
all_uk <- all_uk %>% filter(date>date_pick)

data_proportion <- head(data_proportion,-2) # Remove last 1 days
#data_proportion$long_dates <- as.Date(data_proportion$long_dates,origin="1970-01-01")

data_india$date_time <- as.Date(data_india$date_time)
data_india <- data_india %>% filter(date_time>date_pick)

# Importation data

# Imported cases
traveller_cases0 <- traveller_cases0 %>% filter(`Travel Indicator` == "Traveller")
traveller_cases <- traveller_cases0 #head(traveller_cases0,-5)

daily_india_seq <- 0*all_india$cases_new # Imports based on traveller cases
daily_india_seq[match(traveller_cases0$Date,all_india$date)] <- traveller_cases0$Number_B_1_617_2


# Scale imports
south_asia_imports <- 1000 # Assumed initial normalisation
total_cases_india <- all_india %>% filter(country == "India", date<as.Date("2021-04-23")) %>% select(cases_new) %>% sum()
travel_multiplier <- south_asia_imports/total_cases_india

all_india <- all_india %>% mutate(daily_imports = travel_multiplier*cases_new)

# Downweight recent imports based on incubation period (from McAloon et al. 2020)
red_list_point <- as.numeric(india_red_list-date_pick)
total_days <- length(all_india$date)
downweight_imports <- c(rep(1,red_list_point),1-plnorm(1:(total_days-red_list_point),mean=log(5.1),sd=0.5) )

# Format data
ma_India_variant0 <- ma(data_india$B.1.617.2,7) # Moving average of UK cases
l_IV <- length(ma_India_variant0)
final_value <- 1 #tail(ma_India_variant0,4)[1] #1 #
ma_India_variant <- c(ma_India_variant0[1:(l_IV-4)],rep(final_value,total_days-l_IV+4))
ma_India_variant[1:3] <- 0

# Moving average and scaling factors
ma_UK_cases <- ma(all_uk$cases_new,7) # Moving average of UK cases

data_fit <- data_proportion[match(all_uk$date,data_proportion$sample_date),]
t_fit <- nrow(data_fit)

# Set up for fitting
all_uk_fit <- all_uk %>% filter(date<date_uk_fit)
total_days_uk <- length(all_uk_fit$date)
ma_UK_cases_fit <- ma(all_uk_fit$cases_new,7)


# Define parameters
# serial interval from Rai et al: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7448781/
theta_f <- list(serial_mean=log(5.4),serial_sd=0.4)

#incubation_period <- EpiNow2::get_incubation_period(disease = "SARS-CoV-2", source = "lauer")

# Construct table of serial intervals
t_max <- total_days
serial_mat0 <- matrix(0,nrow=t_max,ncol=(t_max+21))
for(ii in 1:t_max){
  serial_mat0[ii,ii:(ii+20)] <- dlnorm(0:20,theta_f[["serial_mean"]],theta_f[["serial_sd"]])
}
serial_mat0 <- serial_mat0[,1:total_days]

#Estimates from PHE report (Table 9)
#r_phe_report <- c(travel = 128/(250*0.724),non_travel = 98/(287*0.805)) # Report 11
r_phe_report <- c(travel = 135/(331*0.704),non_travel = 245/(698*0.818)) # Report 12

# Variant of concern date (UK - 7th May)
voc_date <- as.Date("2021-05-07")




