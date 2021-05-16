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


# Simulation and likelihood code ---------------------------------------------------------

# Simulate outbreaks - 
# Loop for R estimation

fit_R <- function(r_pick,run_n,add_days=25,import_f = under_factor, daily_decline){
  
  # DEBUG r_pick=1.6; run_n=15; add_days=0; import_f = 0.2; daily_decline=0.02
  
  # Estimated India daily B.1.617.2 imports - XX DEBUGGING
  daily_india <- round(downweight_imports*all_india$daily_imports*ma_India_variant) #ma_India_variant
  #daily_india <- rpois(length(daily_india),lambda=daily_india) # Add Poisson noise
  #daily_india <- daily_india_seq
  
  
  # Define parameters
  # serial interval from Rai et al: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7448781/

  theta <- list(r=r_pick,k=kk_pick,serial_mean=log(5.4),serial_sd=log(1.5))
  #plot(0:20,dlnorm(0:20,(theta[["serial_mean"]]),(theta[["serial_sd"]])))
  
  # Set up initial cases
  seed_list <- NULL
  for(ii in 1:length(daily_india)){
    seed_list <- c(seed_list,rep(ii,daily_india[ii]))
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
  declining_cases <- ma_UK_cases_fit[total_days_uk-4]*(1-daily_decline)^(1:(t_max-total_days_uk+4)) # Add Poisson noise?
  ma_UK_cases_2 <- c(ma_UK_cases_fit[1:(total_days_uk-4)],declining_cases)
  
  # Store_values
  store_vals <- data.frame(matrix(NA, nrow = runs,ncol= t_max))
  store_generations <- data.frame(matrix(NA, nrow = runs,ncol= t_max))
  
  for(i in 1:runs) {
    cases <- seed
    
    tt <- seed_list # XX Edit seed cases based on timeseries
    
    timeT <- tt
    gen_store <- c(0,length(seed_list)); genT <- 0 # Store generation of infection
    
    while(cases > 0) { # Stop once out of time
      
      # Lower SAR in onwards generations?
      thetaR_in <- theta[["r"]]
      thetaR_in <- theta[["r"]]*(genT<1) + import_f*theta[["r"]]*(genT>=1) #
      
      secondary0 <- rnbinom(cases,size=theta[["k"]],mu=thetaR_in) 
      secondary <- secondary0[secondary0>0]; tt <- tt[secondary0>0] # Remove zeros
      
      t.new <- numeric()
      
      # Check whether any overall transmission
      if(length(secondary)>0){
        aa <- Sys.time()
        # for(j in 1:length(secondary)) {
        #   t.new <- c(t.new,tt[j] + rlnorm(secondary[j],theta[["serial_mean"]],theta[["serial_sd"]]) )
        # }
        t.new <- rep(tt,secondary) + rlnorm(sum(secondary),theta[["serial_mean"]],theta[["serial_sd"]])
        
        t.new <- t.new[t.new<=(t_max+1)] # Constrain to be before end of simulation
      }

      cases <- length(t.new) 
      tt <- t.new
      timeT <- c(timeT,t.new)
      
      genT <- genT + 1
      gen_store <- rbind(gen_store,c(genT,cases))

    } 
    
    # Tally incidence:
    tally_times <- table(round(timeT))
    
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
  
  expected_cases <- (ma_UK_cases_2); expected_cases <- expected_cases[1:length(actual_prop)]
  expected_6172 <- mean_val; expected_6172 <- expected_6172[1:length(actual_prop)]
  
  # DEBUG
  #plot(data_fit1$B.1.617.2/data_fit1$N);lines(expected_prop)
  #expected_prop <- mean_val/(ma_UK_cases_2); expected_prop <- expected_prop[1:length(actual_prop)]
  
  # XX CONVERT TO HYPERGEOMETRIC?
  # Calculate likelihood on proportions
  #log_L_b <- dhyper(x=data_fit1$B.1.617.2,m=round(expected_6172),n=round(expected_cases),k=data_fit1$N,log=T)
  
  log_L_b <- dbinom(data_fit1$B.1.617.2,data_fit1$N,expected_prop,log=T)
  log_L_b_sum <- log_L_b[!is.na(log_L_b)] %>% sum()
  
  # Calculate likelihood on cases
  log_L_p <- dpois(round(ma_UK_cases),lambda=mean_val+ma_UK_cases_2,log=T)
  log_L_p_sum <- log_L_p[!is.na(log_L_p)] %>% sum()
  
  # Total likelihood
  log_L_sum <- log_L_b_sum + log_L_p_sum 
  
  # Output fits
  list(lik = log_L_sum, traj = pred_interval,long_dates = long_dates, mov_average = ma_UK_cases_2, daily_india = daily_india)
  
  
}

# End loop




# Run inference -----------------------------------------------------------

run_inference <- function(){

  #rr_store <- NULL
  
  # Run grid search
  
  estimate_list <- foreach(ii = 1:nrow(parameter_list)) %dopar% {
    # for(jj in 1:decline_range){
    #   output1 <- fit_R(rr_range[ii],run_n,add_days = 0, import_est = daily_india, daily_decline = decline_range[jj])
    #   c(rr_range[ii],output1$lik)
    # }
    output1 <- fit_R(parameter_list$rr[ii],run_n,add_days = 0, import_f = parameter_list$imp[ii], daily_decline=parameter_list$decline[ii])
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
  rr_store <- rr_store[rr_store$lik!=-Inf,] #<- -1e5 # Remove invalid points
  # 
  # # Fit spline to estimates
  # modelB.P <- gam(lik ~ s(rr,k=8) + s(decline,k=4) + s(imp,k=4) , data = rr_store,family = "gaussian") 
  # x1_list <- seq(min(rr_store$rr),max(rr_store$rr),0.01)
  # x2_list <- seq(min(rr_store$decline),max(rr_store$decline),0.0005)
  # x3_list <- seq(min(rr_store$imp),max(rr_store$imp),0.01)
  # 
  # param_long <- expand.grid(x1_list, x2_list, x3_list)
  # param_long <- data.frame(param_long); names(param_long) <- names(parameter_list)
  # 
  # preds <- predict(modelB.P, newdata = list(rr=param_long$rr,decline=param_long$decline,imp=param_long$imp), type = "link", se.fit = TRUE)
  # 
  # # Spilne fit to surface
  # surface_est <- cbind(param_long,preds$fit); names(surface_est) <- c(names(parameter_list),"lik")
  # 
  # # DEBUG  - PLOT SOME SLICES
  # #surface_est %>% filter(imp==1)
  # rr_store %>% filter(rr==1.4,imp==1, decline==0.03)
  # 
  
  # - - -
  # Profile individual parameters:
  # R estimate
  x1 <- unique(rr_store$rr); lik1 <- sapply(x1,function(x){max(rr_store[rr_store$rr==x,"lik"])})
  data_lik <- data.frame(cbind(x1,lik1));   data_lik <- data_lik[data_lik$lik1>(max(data_lik$lik1)-100),] # Remove outliers
  #sapply(unique(rr_store$rr),function(x){max(rr_store[rr_store$rr==x,"lik"])})
  modelB.P <- gam(lik1 ~ s(x1,k=4) , data = data_lik,family = "gaussian") 
  x1_list <- seq(min(x1),max(x1),0.001); preds <- predict(modelB.P, newdata = list(x1=x1_list), type = "link", se.fit = TRUE)

  mle_val_rr <- x1_list[which(preds$fit==max(preds$fit))]; range_95 <- x1_list[which(preds$fit>=(max(preds$fit)-1.92))]
  range_95_rr <- c(min(range_95),max(range_95))
  
  # Decline estimate
  x1 <- unique(rr_store$decline); lik1 <- sapply(x1,function(x){max(rr_store[rr_store$decline==x,"lik"])})
  data_lik <- data.frame(cbind(x1,lik1));   data_lik <- data_lik[data_lik$lik1>(max(data_lik$lik1)-100),] # Remove outliers
  modelB.P <- gam(lik1 ~ s(x1,k=4) , data = data_lik,family = "gaussian") 
  x1_list <- seq(min(x1),max(x1),0.01); preds <- predict(modelB.P, newdata = list(x1=x1_list), type = "link", se.fit = TRUE)

  mle_val_dec <- x1_list[which(preds$fit==max(preds$fit))]; range_95 <- x1_list[which(preds$fit>=(max(preds$fit)-1.92))]
  range_95_decline <- c(min(range_95),max(range_95))
  
  # Import estimate
  x1 <- unique(rr_store$imp); lik1 <- sapply(x1,function(x){max(rr_store[rr_store$imp==x,"lik"])})
  data_lik <- data.frame(cbind(x1,lik1)); data_lik <- data_lik[data_lik$lik1>(max(data_lik$lik1)-100),] # Remove outliers
  modelB.P <- gam(lik1 ~ s(x1,k=4) , data = data_lik,family = "gaussian") 
  x1_list <- seq(min(x1),max(x1),0.001); preds <- predict(modelB.P, newdata = list(x1=x1_list), type = "link", se.fit = TRUE)
  
  mle_val_imp <- x1_list[which(preds$fit==max(preds$fit))]; range_95 <- x1_list[which(preds$fit>=(max(preds$fit)-1.92))]
  range_95_imp <- c(min(range_95),max(range_95))
  
  # Combine estimates
  mle_val <- matrix(c(mle_val_rr,mle_val_dec,mle_val_imp),ncol=3); mle_val <- data.frame(mle_val);  names(mle_val) <- names(parameter_list)

  
  plot(x1,lik1,ylim=c(max(lik1)-50,max(lik1))); lines(x1_list,preds$fit)
  
  #modelB.P <- gam(lik ~ s(rr,k=8) , data = rr_store,family = "gaussian") 
  
  # Estimate MLE
  #mle_val <- surface_est[which(preds$fit==max(preds$fit)),]
  #range_95 <- param_long[which(preds$fit>=mle_val$lik-1.92),]

  # Calculate MLE from grid search
  # mle_val <- rr_store[which(rr_store$lik==max(rr_store$lik)),]
  # range_95 <- rr_store[which(rr_store$lik>=(mle_val$lik-1.92)),]
  
  # Sense check MLE:
  #rr_store %>% filter(rr==1.4,decline==0.02,imp==1)
  #plot(rr_store$lik,ylim=c(max(rr_store$lik)-5,max(rr_store$lik)))
  
  # range_95_rr <- c(min(range_95$rr),max(range_95$rr))
  # range_95_decline <- c(min(range_95$decline),max(range_95$decline))
  # range_95_imp <- c(min(range_95$imp),max(range_95$imp))

  
  list(mle = mle_val, r95_rr = range_95_rr, r95_dec = range_95_decline, r95_imp = range_95_imp)

}




# Inference code - White & Pagano ---------------------------------------------------------

fit__WP <- function(r_pick,run_n,add_days=25,import_f = under_factor, daily_decline){
  
  # DEBUG r_pick=1.5; run_n=20; add_days=0; import_f = under_factor
  
  # Estimated India daily B.1.617.2 imports
  daily_india <- import_f*round(downweight_imports*all_india$daily_imports*ma_India_variant)
  
  # Define parameters
  # serial interval from Rai et al: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7448781/
  
  theta <- list(r=r_pick,k=kk_pick,serial_mean=log(5.4),serial_sd=log(2))
  #plot(0:20,dlnorm(0:20,(theta[["serial_mean"]]),(theta[["serial_sd"]])))
  
  
  # Set up dates
  runs <- run_n
  t_max <- total_days+add_days
  
  daily_india
  
  
  # To do - restructure fits to get raw values?
  
  # Extract UK fit data
  all_uk_fit <- all_uk %>% filter(date<date_uk_fit)
  total_days_uk <- length(all_uk_fit$date)
  ma_UK_cases_fit <- ma(all_uk_fit$cases_new,7)
  ma_UK_cases_2 <- c(ma_UK_cases_fit[1:(total_days_uk-4)],ma_UK_cases_fit[total_days_uk-4]*(1-daily_decline)^(1:(t_max-total_days_uk+4)) )
  
  # Store_values
  store_vals <- data.frame(matrix(NA, nrow = runs,ncol= t_max))
  store_generations <- data.frame(matrix(NA, nrow = runs,ncol= t_max))
  
  for(i in 1:runs) {
    cases <- seed
    
    t <- seed_list # XX Edit seed cases based on timeseries
    
    timeT <- t 
    gen_store <- c(0,length(seed_list)); genT <- 0 # Store generation of infection
    
    while(cases > 0) { # Stop once out of time
      secondary <- rnbinom(cases,size=theta[["k"]],mu=theta[["r"]]) 
      t.new <- numeric()
      for(j in 1:length(secondary)) {
        t.new <- c(t.new,t[j] + rlnorm(secondary[j],theta[["serial_mean"]],theta[["serial_sd"]]) )
      }
      
      t.new <- t.new[t.new<(t_max+1)] # Constrain to be before end of simulation
      cases <- length(t.new) 
      t <- t.new
      timeT <- c(timeT,t.new)
      
      genT <- genT + 1
      gen_store <- rbind(gen_store,c(genT,cases))
      
    } 
    
    # Tally incidence:
    tally_times <- table(round(timeT))
    
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
  list(lik = log_L_sum, traj = pred_interval,long_dates = long_dates, mov_average = ma_UK_cases_2, daily_india = daily_india)
  
  
}
