
# Load local data

store_loc <- NULL

for(local_nn in 1:length(local_names)){
  
  data_in <- read_csv(paste0("outputs/result_",local_nn,".csv"))
  store_loc <- rbind(store_loc,data_in)
  
}

# Table 1

store_loc0 <- store_loc %>% filter(parameter=="R_non_traveller" | parameter=="R_p7" | parameter=="R_p17" | parameter=="R_non_617") 

output_A <- store_loc0 %>% pivot_wider(names_from = parameter, values_from = value)

write_csv(output_A,"outputs/Table_R_outputs.csv")

# Table 2

store_loc0 <- store_loc %>% filter(parameter=="ratio_617_2_recent" | parameter=="ratio_617_2_end") 

output_A <- store_loc0 %>% pivot_wider(names_from = parameter, values_from = value)

write_csv(output_A,"outputs/Table_R_ratio.csv")

# Table 3

store_loc0 <- store_loc %>% filter(parameter=="decline" | parameter=="decline_date") 

output_A <- store_loc0 %>% pivot_wider(names_from = parameter, values_from = value)

write_csv(output_A,"outputs/Table_decline.csv")

# Table 4

store_loc0 <- store_loc %>% filter(parameter=="doubling" | parameter=="surge_report") 

output_A <- store_loc0 %>% pivot_wider(names_from = parameter, values_from = value)

write_csv(output_A,"outputs/Table_report.csv")


# Output summary for that location
c("R_traveller","R_non_traveller","R_recent","R_non_17","ratio_617_2","ratio_617_2_recent","decline","decline_date")
