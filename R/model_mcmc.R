# - - - - - - - - - - - - - - - - - - - - - - - 
# Run MCMC
# - - - - - - - - - - - - - - - - - - - - - - - 


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# RUN MCMC Model

run_transmission_mcmc <- function(MCMC.runs = 10){

  multichain <- c(2) # run in parallel
  #iiM <- 1
  
  # DEBUG  multichain=c(4); MCMC.runs=10; iiM = multichain; prop.risk = F
  
  foreach(iiM=multichain) %dopar% {  # Loop over scenarios with parallel MCMC chains
  #for(iiM in multichain){

  # - - - - - - - - - - - 
  # Load relevant data
  thetaR_IC = read_csv("data/thetaR_IC.csv")
  thetaR_IC <- thetaR_IC[1,]
  
  # - - - - - - - - - - - 
  # Initialise ICs 

  # Global parameters
  theta = c(rr = thetaR_IC$rr,
            r_scale = thetaR_IC$r_scale,
            r_scale_2 = thetaR_IC$r_scale_2,
            decline = thetaR_IC$decline,
            dt_decline = thetaR_IC$dt_decline,
            imp = thetaR_IC$imp,
            rep_vol = thetaR_IC$rep_vol
            )
  
  
  # Covariance matrices - Add theta and thetaAll together in MCMC runs
  nparam = length(theta) 
  npc = rep(1,nparam)
  pmask = NULL # fix parameters
  
  # Define models (i.e. parameters to fit)
  if(iiM==1){pmask <- c("dt_decline")}
  if(iiM==2){pmask <- c("r_scale_2","dt_decline")}
  
  npc[match(pmask,names(theta))]=0
  cov_matrix_theta0 = diag(npc)
  
  # Quick simulation to check IC OK
  output1 <- fit_R_deterministic(theta = theta,run_n,add_days=0)
  sim_marg_lik_star <- output1$lik
  
  # - - - - - - - - - - - 
  # Set up matrices for MCMC run
  # - - - - - - - - - - -
  
  m = 1 #set initial MCMC
  thetatab=matrix(NA,nrow=(MCMC.runs+1),ncol=length(theta))
  colnames(thetatab)=names(theta)
  thetatab[1,]=theta
  
  sim_liktab=rep(-Inf,(MCMC.runs+1))
  accepttab=rep(NA,(MCMC.runs))
  max.length = t_fit  # Need to store enough values
  c_trace_tab=array(NA, dim=c(MCMC.runs+1,max.length)) # Note max length
  intro_trace_tab=array(NA, dim=c(MCMC.runs+1,max.length)) # Note max length
  
  # - - - - - - - - - - - 
  #RUN MCMC:
  # - - - - - - - - - - - 
  
  for (m in 1:MCMC.runs){
    
    if(m==1){
      epsilon0=0.001
      cov_matrix_theta=epsilon0*cov_matrix_theta0
    }else{
      epsilon0 = max(min(0.1,exp(log(epsilon0)+(accept_rate-0.234)*0.999^m)),1e-6) # Stop epsilon getting too big or small
      cov_matrix_theta=epsilon0*cov_matrix_theta0
    }
    
    # Resample global theta
    if(m>1){
      output_theta = SampleTheta(thetatab[m,],cov_matrix_theta,pmask) #sample nearby global parameter space
      theta_star=output_theta$thetaS
    }else{
      theta_star=thetatab[m,]
    }
    

    # RUN MODEL
    if(m==1){
      output1 <- fit_R_deterministic(theta = theta_star,run_n=10,add_days=0) # To avoid zero likelihood
    }else{
      output1 <- fit_R_deterministic(theta = theta_star,run_n=1,add_days=0)
    }

    sim_marg_lik_star <- output1$lik

    # Calculate probability function
    output_prob <- ComputeProbability(sim_liktab[m],sim_marg_lik_star,thetatab[m,],theta_star) 
    
    # Update parameter values
    if(runif(1,0,1) < output_prob){
      thetatab[m+1,] = theta_star
      sim_liktab[m+1] = sim_marg_lik_star
      c_trace_tab[m+1,] = as.numeric(output1$traj[1,])
      intro_trace_tab[m+1,] = as.numeric(output1$daily_india)
      accepttab[m]=1
      
    }else{
      thetatab[m+1,] = thetatab[m,]
      sim_liktab[m+1] = sim_liktab[m]
      c_trace_tab[m+1,] = c_trace_tab[m,]
      intro_trace_tab[m+1,] = intro_trace_tab[m,]
      accepttab[m]=0
    }
    
    # Upate acceptance rate
    if(m<20){
      accept_rate=0.234
    }else{
      accept_rate=sum(accepttab[1:m])/m
    }

    # Save outputs every 1000 iterations
    
    if(m %% min(MCMC.runs,1000) == 0){
      print(c(m,accept_rate,sim_liktab[m],epsilon0))
      save(sim_liktab,accepttab,c_trace_tab,intro_trace_tab,thetatab,file=paste("outputs/outputR",iiM,".RData",sep=""))
    }
    
  } # End MCMC run
  
  } # End multichains
  
}


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Compute acceptance probability

ComputeProbability<-function(sim_likelihood,sim_likelihood_star,thetatab,theta_star,pmask=NULL){
  
  # Include priors
  p_theta_star = priorScale(theta_star[["r_scale"]])*priorScale(theta_star[["r_scale_2"]])
  p_theta = priorScale(thetatab[["r_scale"]])*priorScale(thetatab[["r_scale_2"]])
  

  # Calculate acceptance probability. Adjust for log sampling
  val = exp((sim_likelihood_star-sim_likelihood + sum(log(theta_star)) - sum(log(thetatab))  ))*(p_theta_star/p_theta)
  min(val, 1)
  
  
}


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Resample theta

SampleTheta<-function(theta_in, covartheta, pmask){
  
  # sample new parameters from nearby: 
  mean_vector_theta = theta_in
  mean_vector_theta0 = mean_vector_theta 
  
  theta_star = as.numeric(exp(rmvnorm(1,log(mean_vector_theta0), covartheta)))
  names(theta_star) = names(theta_in)

  return(list(thetaS=theta_star))
  
}


