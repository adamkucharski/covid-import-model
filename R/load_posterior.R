thetatabA=NULL
c_trace_tab0=NULL
intro_trace_tab0=NULL

# Load values
load(paste("outputs/outputR",iiM,".RData",sep=""))
thetatab=cbind(data.frame(thetatab))

# Remove burn-in
mcmc.burn = 0.2
mcmc_samples=length(sim_liktab)
maxB=sum(sim_liktab!=-Inf)/mcmc_samples
minB=mcmc.burn*maxB
picks0=c(round(minB*mcmc_samples):round(maxB*mcmc_samples))

# Compile posteriors
sim_likOut = sim_liktab[picks0]
thetatabA=rbind(thetatabA,thetatab[picks0,])

c_trace_tab0 = rbind(c_trace_tab0,c_trace_tab[picks0,])
intro_trace_tab0 = rbind(intro_trace_tab0,intro_trace_tab[picks0,])

picks=c(1:length(thetatabA[,1]))

pick.max = picks[sim_likOut[picks]==max(sim_likOut[picks])][1] # Maximum likelihood

thetatab=thetatabA
c_trace_tab=c_trace_tab0
intro_trace_tab=intro_trace_tab0

# MLE
thetatab[pick.max,]

#plot(sim_likOut,type="l")

# Extract posteriors

medP <- apply(c_trace_tab,2,function(x){median(x)})
ciP1 <- apply(c_trace_tab,2,function(x){quantile(x,0.025)})
ciP2 <- apply(c_trace_tab,2,function(x){quantile(x,0.975)})

pred_interval <- rbind(medP,ciP1,ciP2)


