# - - 
# Helper functions

c.nume<-function(x){
  bp1=c(median(x),quantile(x,0.025),quantile(x,0.975))
  as.numeric(bp1)
}

doubling_est <- function(x){log((tail(x,1)-head(x,1))/length(x),2)}

doubling_est2 <- function(x){log(2)*length(x)/log(tail(x,1)/head(x,1))}

c.text <- function(x,sigF=2){
  bp1=signif(c(median(x),quantile(x,0.025),quantile(x,0.975)),sigF)
  paste(bp1[1]," (95% CrI: ",bp1[2],"-",bp1[3],")",sep="")
}

c.nume.50<-function(x){
  bp1=c(median(x),quantile(x,0.025),quantile(x,0.25),quantile(x,0.75),quantile(x,0.975))
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

decline_f <- function(total_days_uk,t_max,daily_decline,dt_decline,ma_UK_cases_fit){
  d_list <- (1:(t_max-total_days_uk+4))
  daily_declineX <- daily_decline + (dt_decline-1)*d_list/10
  declining_cases <- ma_UK_cases_fit[total_days_uk-4]*(1-daily_declineX)^d_list # Add Poisson noise?
  
  ma_UK_cases_2 <- c(ma_UK_cases_fit[1:(total_days_uk-4)],declining_cases)
  
  # Alternative version with sequence numbers
  # prop_non_b16172 <- (1-data_fit$B.1.617.2/data_fit$N)
  # prop_non_b16172 <- ma(prop_non_b16172,7); prop_non_b16172[1:3] <- 1 # Fix early points
  # 
  # na_remove <- prop_non_b16172[!is.na(prop_non_b16172)]
  # tail_val <- tail(na_remove,1)
  # prop_non_b16172 <- c(prop_non_b16172[1:length(na_remove)],rep(tail_val,t_max-length(na_remove)))
  # #ma_UK_cases_2 <- ma_UK_cases*prop_non_b16172 # Includes overall UK case data
  # 
  # na_remove_cases <- ma_UK_cases_fit[1:(total_days_uk-4)]
  # ma_UK_cases0 <- c(na_remove_cases,rep(tail(na_remove_cases,1),t_max-length(na_remove_cases)))
  # ma_UK_cases_2 <- ma_UK_cases0*prop_non_b16172 # Includes overall UK case data
  
}



# Simulation deterministic  ---------------------------------------------------------

fit_R_deterministic <- function(theta,run_n,add_days=25){
  
  # DEBUG r_pick=1.6; run_n=15; add_days=0; import_f = 0.2; daily_decline=0.02
  # theta <- c(rr=1.6,r_scale=1,r_scale_2=1,decline=0.02,dt_decline=1,imp=2,rep_vol=1); add_days=0
  
  r_pick <- theta[["rr"]]
  r_scale <- theta[["r_scale"]]
  r_scale_2 <- theta[["r_scale_2"]]
  import_f <- theta[["imp"]]
  daily_decline <- theta[["decline"]]
  dt_decline <- 1 #theta[["dt_decline"]] # Scale rate of change
  rep_scale <- theta[["rep_scale"]]
  rep_vol <- theta[["rep_vol"]]
  rep_vol_seq <- theta[["rep_vol_seq"]]
  surge_scale <- theta[["surge_scale"]]
  surge_time <- theta[["surge_time"]]
  end_scale <- theta[["end_scale"]]
  end_time <- theta[["end_time"]]
  
  # Estimated India daily B.1.617.2 imports 
  daily_india <- (import_f*downweight_imports*all_india$daily_imports*ma_India_variant) #ma_India_variant
  
  #daily_india <- daily_india_seq + import_f*downweight_imports*all_india$daily_imports*ma_India_variant
  
  daily_india <- c(daily_india,rep(0,add_days)) # Add more initial points 
  #daily_india <- rpois(length(daily_india),lambda=daily_india) # Add Poisson noise

  #plot(0:20,dlnorm(0:20,(theta[["serial_mean"]]),(theta[["serial_sd"]])))

  # Set up dates
  t_max <- total_days+add_days
  if(add_days>0){
    long_dates <- c(all_uk$date,max(all_uk$date)+1:add_days)
    # Set up matrix
    serial_mat <- matrix(0,nrow=t_max,ncol=(t_max+21))
    for(ii in 1:t_max){
      serial_mat[ii,ii:(ii+20)] <- dlnorm(0:20,theta_f[["serial_mean"]],theta_f[["serial_sd"]])
    }
    serial_mat <- serial_mat[,1:t_max]
  }else{
    long_dates <- all_uk$date
    serial_mat <- serial_mat0
  }

  # Extract UK fit data
  ma_UK_cases_2 <- decline_f(total_days_uk,t_max,daily_decline,dt_decline,ma_UK_cases_fit)
  
  # Get VOC date - or extract while fitting
  #voc_pick <- which(long_dates>voc_date) %>% min()
  voc_pick <- surge_time
  end_pick <- end_time
  
  # Store_values
  store_vals <- rep(0,t_max)

  # Add onwards transmission from travellers
  serial_mat2 <- serial_mat*daily_india*r_pick
  store_vals1 <- store_vals + colSums(serial_mat2)
  
  # Add onwards transmission from first genertion - exclude travellers here to avoid double counting
  serial_mat2 <- serial_mat*store_vals1*r_pick*r_scale
  store_vals <- store_vals + colSums(serial_mat2)
  
  # Add onwards transmission from subsequent generations - reduces twice
  for(ii in 1:t_max){
    if(ii<voc_pick){r_pick_c <- r_pick} # Include possible effect of VOC measures
    if(ii>=voc_pick){r_pick_c <- r_pick*surge_scale} # Include possible effect of VOC measures (r_scale_2 = 1)
    if(ii>=end_pick){r_pick_c <- r_pick*surge_scale*end_scale} # Include possible effect of VOC measures (r_scale_2 = 1)
    store_vals <- store_vals + store_vals[ii]*serial_mat[ii,]*r_pick_c*r_scale*r_scale_2
  }
  
  store_vals <- daily_india + store_vals1 + store_vals # Add imports + first generation to avoid double counting
  
  # Calculate estimation interval
  pred_interval <- rbind(store_vals,store_vals,store_vals)
  pred_interval[,voc_n:t_max] <- pred_interval[,voc_n:t_max]*rep_scale # Scale by surge testing
  
  #mean_val <-  store_vals
  mean_val <-  pred_interval[1,] # Include scale by surge testing

  # Observation model
  data_fit1 <- data_fit[!is.na(data_fit$N),]
  actual_617 <- data_fit1$B.1.617.2
  actual_prop <- data_fit1$B.1.617.2/data_fit1$N
  actual_non_6172 <- data_fit1$N - data_fit1$B.1.617.2
  
  expected_prop <- mean_val/(mean_val+ma_UK_cases_2); expected_prop <- expected_prop[1:length(actual_prop)]
  
  expected_cases <- (ma_UK_cases); expected_cases <- expected_cases[1:length(actual_prop)]
  expected_6172 <- mean_val; expected_6172 <- expected_6172[1:length(actual_prop)]
  
  
  # DEBUG
  #plot(data_fit1$B.1.617.2/data_fit1$N);lines(expected_prop)
  #expected_prop <- mean_val/(ma_UK_cases_2); expected_prop <- expected_prop[1:length(actual_prop)]
  
  #
  # - - 
  # Calculate likelihood on 617.2
  log_L_b <- dnbinom(actual_617,mu=(expected_6172*data_fit1$N/expected_cases),size=1/rep_vol_seq,log=T)

  log_L_b_sum <- log_L_b[!is.na(log_L_b)] %>% sum()
  
  # Calculate likelihood of overall cases
  pick_fit <- all_uk$date>fit_date # Only from April 1st 2021
  
  log_L_p <- dnbinom(all_uk$cases_new,mu=(mean_val+ma_UK_cases_2),size=1/rep_vol,log=T)

  log_L_p <- log_L_p[pick_fit]
  
  log_L_p_sum <- log_L_p[!is.na(log_L_p)] %>% sum()
  
  # Total likelihood
  log_L_sum <- log_L_b_sum + log_L_p_sum 
  
  # Output fits
  list(lik = log_L_sum, traj = pred_interval,long_dates = long_dates, mov_average = ma_UK_cases_2, daily_india = daily_india)
  
  
}

# End loop


