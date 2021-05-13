# - - 
# Helper functions

c.nume<-function(x){
  bp1=c(median(x),quantile(x,0.025),quantile(x,0.975))
  as.numeric(bp1)
}

# Plot data and binom CI
plot_CI <- function(dates,xx,nn,colA="black") {
  
  for(ii in 1:length(nn)){
    test_r <- binom.test(xx[ii],nn[ii])
    CI1 <- as.numeric(test_r$conf.int)[1]
    CI2 <- as.numeric(test_r$conf.int)[2]
    points(dates[ii],xx[ii]/nn[ii],col=colA,pch=19); lines(c(dates[ii],dates[ii]),c(CI1,CI2),col=colA)
  }
  
}


# Simulation code ---------------------------------------------------------

# Simulate outbreaks - serial interval from Rai et al: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7448781/

# Loop for R estimation

fit_R <- function(r_pick,run_n,add_days=25,import_est = daily_india, daily_decline){
  
  # DEBUG r_pick=1.5; run_n=20; add_days=0; import_est = daily_india
  
  theta <- list(r=r_pick,k=kk_pick,serial_mean=log(5.4),serial_sd=log(1.5))
  # plot(0:20,dlnorm(0:20,(theta[["serial_mean"]]),(theta[["serial_sd"]])))
  
  # Set up initial cases
  seed_list <- NULL
  for(ii in 1:length(import_est)){
    seed_list <- c(seed_list,rep(ii,import_est[ii]))
  }
  
  seed <- length(seed_list)
  
  # Set up dates
  runs <- run_n
  t_max <- total_days+add_days
  if(add_days>0){
    long_dates <- c(all_uk$date,max(all_uk$date)+1:add_days)
  }else{
    long_dates <- all_uk$date
  }
  

  # Extract UK fit data
  all_uk_fit <- all_uk %>% filter(date<date_uk_fit)
  total_days_uk <- length(all_uk_fit$date)
  ma_UK_cases_fit <- ma(all_uk_fit$cases_new,7)
  ma_UK_cases_2 <- c(ma_UK_cases_fit[1:(total_days_uk-4)],ma_UK_cases_fit[total_days_uk-4]*(1-daily_decline)^(1:(t_max-total_days_uk+4)) )
  
  # Store_values
  
  store_vals <- data.frame(matrix(NA, nrow = runs,ncol= t_max))
  
  for(i in 1:runs) {
    cases <- seed
    
    t <- seed_list # XX Edit seed cases based on timeseries
    
    times <- t 
    while(cases > 0) { # Stop once out of time
      secondary <- rnbinom(cases,size=theta[["k"]],mu=theta[["r"]]) 
      t.new <- numeric()
      for(j in 1:length(secondary)) {
        t.new <- c(t.new,t[j] + rlnorm(secondary[j],theta[["serial_mean"]],theta[["serial_sd"]]) )
      }
      
      t.new <- t.new[t.new<(t_max+1)] # Constrain to be before end of simulation
      cases <- length(t.new) 
      t <- t.new
      times <- c(times,t.new)
    } 
    
    # Tally times:
    tally_times <- table(round(times))
    
    tally_times
    time_ID <- as.numeric(names(tally_times))
    
    #store_times <- rep(0,total_days)
    store_times <- as.numeric(tally_times[match(1:t_max,time_ID)])
    store_times[is.na(store_times)] <- 0
    
    store_vals[i,] <- store_times
    
  }
  
  # Calculate estimation interval
  pred_interval <- apply(store_vals,2,function(x){c.nume(x)})
  mean_val <-  apply(store_vals,2,mean)
  
  # Observation model
  data_fit1 <- data_fit[!is.na(data_fit$N),]
  actual_prop <- data_fit1$B.1.617.2/data_fit1$N
  expected_prop <- mean_val/(mean_val+ma_UK_cases_2); expected_prop <- expected_prop[1:length(actual_prop)]
  
  #expected_prop <- mean_val/(ma_UK_cases_2); expected_prop <- expected_prop[1:length(actual_prop)]
  
  # Calculate likelihood on proportions
  log_L_b <- dbinom(data_fit1$B.1.617.2,data_fit1$N,expected_prop,log=T)
  log_L_b_sum <- log_L_b[!is.na(log_L_b)] %>% sum()
  
  # Calculate likelihood on cases
  log_L_p <- dpois(round(ma_UK_cases),lambda=mean_val+ma_UK_cases_2,log=T)
  log_L_p_sum <- log_L_p[!is.na(log_L_p)] %>% sum()
  
  # Total likelihood
  log_L_sum <- log_L_b_sum + log_L_p_sum 
  
  # Output fits
  list(lik = log_L_sum, traj = pred_interval,long_dates = long_dates, mov_average = ma_UK_cases_2)
  
  
}

# End loop


# Run inference -----------------------------------------------------------

fit_model <- function(){

  #rr_store <- NULL
  
  # Run grid search
  
  estimate_list <- foreach(ii = 1:nrow(parameter_list)) %dopar% {
    # for(jj in 1:decline_range){
    #   output1 <- fit_R(rr_range[ii],run_n,add_days = 0, import_est = daily_india, daily_decline = decline_range[jj])
    #   c(rr_range[ii],output1$lik)
    # }
    output1 <- fit_R(parameter_list$rr[ii],run_n,add_days = 0, import_est = daily_india, daily_decline=parameter_list$decline[ii])
    output1$lik
  }
  
  rr_store <- cbind(parameter_list,unlist(estimate_list)); names(rr_store) <- c(names(parameter_list),"lik")
  
  # for(ii in 1:length(rr_range)){
  #   output1 <- fit_R(rr_range[ii],run_n,add_days = 0, import_est = daily_india)
  #   rr_store <- rbind(rr_store,c(rr_range[ii],output1$lik))
  # }

  
  #plot(rr_store[,1],rr_store[,2],xlab="R",ylab="lik")
  
  write_csv(data.frame(rr_store),paste0("outputs/fit",kk_pick,".csv"))

}


# Extract fit and estimate MLE --------------------------------------------

get_MLE <- function(){

  # Get MLE
  rr_store <- read_csv(paste0("outputs/fit",kk_pick,".csv"))
  
  # Fit spline to estimates
  modelB.P <- gam(lik ~ s(rr,k=8) + s(decline,k=2) , data = rr_store,family = "gaussian") 
  x1_list <- seq(min(parameter_list$rr),max(parameter_list$rr),0.01)
  x2_list <- seq(min(parameter_list$decline),max(parameter_list$decline),0.001)
  
  param_long <- expand.grid(x1_list, x2_list)
  param_long <- data.frame(param_long); names(param_long) <- names(parameter_list)
    
  preds <- predict(modelB.P, newdata = list(rr=param_long$rr,decline=param_long$decline), type = "link", se.fit = TRUE)
  
  mle_val <- param_long[which(preds$fit==max(preds$fit)),]
  
  range_95 <- xx_list[which(preds$fit>=(max(preds$fit)-1.92))]
  range_95 <- c(min(range_95),max(range_95))
  
  #plot(rr_store$X1,rr_store$X2,xlim=c(min(rr_range),max(rr_range)),ylim=c(-220,-130),xlab="R",ylab="lik",main="R estimate")
  #lines(xx_list,preds$fit)
  
  # Run fitted model
  best_r <- mle_val
  add_d <- 10
  output1 <- fit_R(best_r,run_n,add_days = add_d, import_est = daily_india, daily_decline)

}


