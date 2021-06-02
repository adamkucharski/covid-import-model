
# - - -
# Load covidregionaldata

if(!exists("all_uk")){
  all_india0 <- get_national_data("india") # Load data from covidregionaldata (slow)
  all_uk <- get_national_data("uk",source="who") # Load data from covidregionaldata (slow)
  
  all_uk_region <- get_regional_data("uk")
  
  # Fix missing data in WHO file
  missing_data <- c(2027,1926,1979)
  all_uk[all_uk$date==as.Date("2021-05-16"),"cases_new"] <- missing_data[1]
  all_uk[all_uk$date==as.Date("2021-05-17"),"cases_new"] <- missing_data[2]
  all_uk[all_uk$date==as.Date("2021-05-18"),"cases_new"] <- missing_data[3]
}

# - - -
# Load UK local data

if(!exists("all_uk_p2_in")){
  # all_uk_p2 <- read_csv(paste0(data_path,"merged_2021_06_01.csv"))
  # 
  # all_uk_p2$specimen_date <- lubridate::dmy(all_uk_p2$specimen_date) # convert dates
  # 
  # all_uk_p2a <- all_uk_p2 %>% filter(specimen_date>as.Date("2021-01-01"))
  # all_uk_p2a <- all_uk_p2a %>% select(NHSER_name,specimen_date,s_sgtf,v_variant,v_seq_result,v_specimen_date_sk)
  # 
  # write_csv(all_uk_p2a,paste0(data_path,"merged_2021_06_01_out.csv"))
  
  all_uk_p2_in <- read_csv(paste0(data_path,"merged_2021_06_01_out.csv"))


}

# - - -
# Load downloaded data

if(!exists("data_proportion")){
    
  cog0 = fread(paste0(data_path,"cog_metadata.csv"))[country == "UK" & is_pillar_2 == "Y"]
  
  cog1 <- cog0 %>% filter(sample_date==as.Date("2021-05-15"))
  table(cog1$lineage)
  
  cog = cog0[, .(.N,
                 B.1.617.1 = sum(lineage %like% "B\\.1\\.617\\.1"),
                 B.1.617.2 = sum(lineage %like% "B\\.1\\.617\\.2"),
                 B.1.617.3 = sum(lineage %like% "B\\.1\\.617\\.3")), keyby = sample_date]
  fwrite(cog, paste0(data_path,"COG_UK_out.csv"))

}

data_proportion <- read_csv(paste0(data_path,"COG_UK_out.csv"))

traveller_cases0 <- read_csv("data/voc_imports_2021_05_24.csv") # Read in from repo
traveller_cases_617_1 <- read_csv("data/vui_617_1_imports_2021_05_24.csv") # Read in from repo

data_india <- read_tsv("data/outbreakinfo_mutation_report_data.tsv")

# Technical report 13 data:
tr_fig1 <- read_csv(paste0(data_path,"TR13_figure1.csv"))
tr_fig4 <- read_csv(paste0(data_path,"TR13_figure4.csv"))
tr_fig10 <- read_csv(paste0(data_path,"TR13_figure10.csv"))

#all_uk_p2_in %>% filter(v_variant=="VOC-21APR-02") %>% group_by(NHSER_name) %>% tally()

# Redo file load with local data
if(local_run==TRUE){
  
  local_pick <- location_ID
  
  d1_local <- all_uk_p2_in %>% filter(NHSER_name==local_pick) # Get local data
  d1_sequenced <- d1_local %>% filter(!is.na(v_variant))
  d1_variant <- d1_sequenced %>% filter(v_variant=="VOC-21APR-02")

  # Count sequences
  data_proportion0 <- data_proportion
  data_proportion0$N <- NA; data_proportion0$B.1.617.2 <- 0
  
  d1_s_tab <- table(d1_sequenced$v_specimen_date_sk) # Tally up sequenced
  d1_s_pick <- match(data_proportion0$sample_date,as.Date(names(d1_s_tab)))
  data_proportion0$N <- as.numeric(d1_s_tab[d1_s_pick]) # EDIT

  d1_v_tab <- table(d1_variant$v_specimen_date_sk) # Tally up sequenced
  d1_v_pick <- match(data_proportion0$sample_date,as.Date(names(d1_v_tab)))
  data_proportion0$B.1.617.2 <- as.numeric(d1_v_tab[d1_v_pick]) # EDIT
  
  # Amend all data
  all_uk <- all_uk_region %>% filter(region==local_pick) # Amend to have local data
  if(local_pick=="Midlands"){
    all_uk_1 <- all_uk_region %>% filter(region=="West Midlands" | region=="East Midlands") # Amend to have local data
    all_uk <- all_uk_1 %>% group_by(date) %>% summarise(cases_new = sum(cases_new))
  }
  if(local_pick=="North East and Yorkshire"){
    all_uk <- all_uk_region %>% filter(region=="Yorkshire and The Humber" | region=="North East") # Amend to have local data
    all_uk <- all_uk_1 %>% group_by(date) %>% summarise(cases_new = sum(cases_new))
  }

  # Amend sequence data
  # tr_617_2_local <- tr_fig10 %>% filter(Region==local_pick)
  # d1_s_pick <- match(tr_617_2_local$Date,data_proportion$sample_date)
  # 
  # 
  # 
  # data_proportion[d1_s_pick,]$B.1.617.2 <- tr_617_2_local$Number_VOC # Add VOC count
  # 
  # # Add N sequenced
  # tr_N_local <- tr_fig4 %>% filter(Region==local_pick)
  # d1_n_1 <- match(data_proportion$sample_date,tr_N_local$Date) # d1_n_1 <- d1_n_1[!is.na(d1_n_1)]
  # d1_n_2 <- match(data_proportion$sample_date,all_uk$date) # d1_n_1 <- d1_n_1[!is.na(d1_n_1)]
  #       
  # data_proportion$N <- all_uk[d1_n_2,]$cases_new # Add cases
  # 
  # data_proportion$N <- round(data_proportion$N*0.01*tr_N_local[d1_n_1,]$seven_day_pr) #Add sequenced %
  
  # REMOVE NA
  data_proportion <- data_proportion0
  data_proportion <- data_proportion[!is.na(data_proportion$N),]
  all_uk <- all_uk[!is.na(all_uk$cases_new),]
  # 
  # all_uk %>% filter(date==as.Date("2021-05-16"))
  # tr_N_local %>% filter(Date==as.Date("2021-05-16"))
  
}





# - - -
# Load variant data

india_red_list <- as.Date("2021-04-23")
date_pick <- as.Date("2021-02-01")
date_uk_fit <- as.Date("2021-04-23")

fit_date <- as.Date("2021-04-01") # Specify date to fit from


all_uk <- head(all_uk,-1) # Remove final point because of bank holiday

all_india <- all_india0 %>% filter(date>=date_pick)
all_uk <- all_uk %>% filter(date>=date_pick)

all_india <- all_india[1:nrow(all_uk),] # Avoid mismatched lengths

if(local_run==T){
  data_proportion <- head(data_proportion,-2) # Remove last X days
}else{
  data_proportion <- head(data_proportion,-2) # Remove last X days
}
#data_proportion$long_dates <- as.Date(data_proportion$long_dates,origin="1970-01-01")

data_india$date_time <- as.Date(data_india$date_time)
data_india <- data_india %>% filter(date_time>date_pick)

# Importation data

# Imported cases
traveller_cases0 <- traveller_cases0 %>% filter(`Travel Indicator` == "Traveller")
traveller_cases <- traveller_cases0 #head(traveller_cases0,-5)



#daily_india_seq <- 0*all_india$cases_new # Imports based on traveller cases
#daily_india_seq[match(traveller_cases0$Date,all_india$date)] <- traveller_cases0$Number_B_1_617_2


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

# Set up for extrapolation
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

voc_n <- as.numeric(voc_date - min(all_uk$date) +1)

fit_n <- as.numeric(fit_date - min(all_uk$date) +1)





