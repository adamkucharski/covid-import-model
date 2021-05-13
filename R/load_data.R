# - - -
# Load country data (slow)

# all_countries <- get_national_data() 

# - - -
# Load variant data

india_red_list <- as.Date("2021-04-23")
date_pick <- as.Date("2021-02-01")
date_uk_fit <- as.Date("2021-04-23")

all_india <- all_countries %>% filter(country == "India", date>date_pick)
all_uk <- all_countries %>% filter(country == "United Kingdom",date>date_pick)

data_proportion <- read_csv(paste0(data_path,"data/B.1.617-COG-UK-20210513.csv"))
data_proportion <- head(data_proportion,-2) # Remove last 2 days
#data_proportion$long_dates <- as.Date(data_proportion$long_dates,origin="1970-01-01")

data_india <- read_tsv(paste0(data_path,"data/outbreakinfo_mutation_report_data_2021-05-13.tsv"))
data_india$date_time <- as.Date(data_india$date_time)
data_india <- data_india %>% filter(date_time>date_pick)

# Importation data

traveller_cases0 <- read_csv(paste0(data_path,"data/voc_imports.csv"))
traveller_cases0 <- traveller_cases0 %>% filter(`Travel Indicator` == "Traveller")
traveller_cases <- head(traveller_cases0,-5)