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
  rep_vol <- theta[["rep_vol"]]
  
  # Estimated India daily B.1.617.2 imports - XX DEBUGGING
  daily_india <- (import_f*downweight_imports*all_india$daily_imports*ma_India_variant) #ma_India_variant
  
  daily_india <- daily_india_seq + import_f*downweight_imports*all_india$daily_imports*ma_India_variant
  
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
  
  # Store_values
  store_vals <- rep(0,t_max)

  # Add onwards transmission from travellers
  serial_mat2 <- serial_mat*daily_india*r_pick
  store_vals1 <- store_vals + colSums(serial_mat2)
  
  # Add onwards transmission from first genertion - exclude travellers for now to avoid double counting
  serial_mat2 <- serial_mat*store_vals1*r_pick*r_scale
  store_vals <- store_vals + colSums(serial_mat2)
  
  
  # Add onwards transmission from subsequent generations - reduces twice
  for(ii in 1:t_max){
    store_vals <- store_vals + store_vals[ii]*serial_mat[ii,]*r_pick*r_scale*r_scale_2
  }
  
  store_vals <- daily_india + store_vals1 + store_vals # Add imports + first generation to avoid double counting
  
  # Calculate estimation interval
  pred_interval <- rbind(store_vals,store_vals,store_vals)
  mean_val <-  store_vals
  
  # Observation model
  data_fit1 <- data_fit[!is.na(data_fit$N),]
  actual_617 <- data_fit1$B.1.617.2
  actual_prop <- data_fit1$B.1.617.2/data_fit1$N
  actual_non_6172 <- data_fit1$N - data_fit1$B.1.617.2
  
  expected_prop <- mean_val/(mean_val+ma_UK_cases_2); expected_prop <- expected_prop[1:length(actual_prop)]
  
  expected_cases <- (ma_UK_cases); expected_cases <- expected_cases[1:length(actual_prop)]
  #expected_cases <- (ma_UK_cases_2); expected_cases <- expected_cases[1:length(actual_prop)]
  expected_6172 <- mean_val; expected_6172 <- expected_6172[1:length(actual_prop)]
  
  
  # DEBUG
  #plot(data_fit1$B.1.617.2/data_fit1$N);lines(expected_prop)
  #expected_prop <- mean_val/(ma_UK_cases_2); expected_prop <- expected_prop[1:length(actual_prop)]
  
  # XX CONVERT TO HYPERGEOMETRIC?
  # - - 
  # Calculate likelihood on 617.2
  #log_L_b <- dhyper(x=data_fit1$B.1.617.2,m=round(expected_6172),n=round(expected_cases),k=data_fit1$N,log=T)
  
  log_L_b <- dnbinom(actual_617,mu=(expected_6172*data_fit1$N/expected_cases),size=1/rep_vol,log=T)
  
  #log_L_b <- dbinom(data_fit1$B.1.617.2,data_fit1$N,expected_prop,log=T)
  log_L_b_sum <- log_L_b[!is.na(log_L_b)] %>% sum()
  
  # Calculate likelihood on cases
  #log_L_p <- dpois(round(ma_UK_cases),lambda=mean_val+ma_UK_cases_2,log=T)

  #log_L_bP <- dhyper(x=actual_non_6172,m=round(expected_cases),n=round(expected_6172),k=data_fit1$N,log=T)

  #log_L_p <- dpois(all_uk$cases_new,lambda=(mean_val+ma_UK_cases_2),log=T)
  pick_fit <- all_uk$date>as.Date("2021-04-01") # Only ffrom April 1st 2021
  
  log_L_p <- dnbinom(all_uk$cases_new,mu=(mean_val+ma_UK_cases_2),size=1/rep_vol,log=T)
  #log_L_p <- dpois(all_uk$cases_new,lambda=(mean_val+ma_UK_cases_2),log=T)
  log_L_p <- log_L_p[pick_fit]
  
  log_L_p_sum <- log_L_p[!is.na(log_L_p)] %>% sum()
  
  # Total likelihood
  log_L_sum <- log_L_b_sum + log_L_p_sum 
  
  # Output fits
  list(lik = log_L_sum, traj = pred_interval,long_dates = long_dates, mov_average = ma_UK_cases_2, daily_india = daily_india)
  
  
}

# End loop


