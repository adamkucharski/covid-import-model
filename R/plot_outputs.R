# Plot data ---------------------------------------------------------------

rr_store <- read_csv(paste0("outputs/fit",kk_pick,".csv"))

par(mfrow=c(2,4),mar=c(3,3,1,1),mgp=c(2,0.6,0),las=0)

x_range <- c(as.Date("2021-03-01"),max(all_uk$date)+25)
col2a <- rgb(1,0.2,0,0.2)
col2c <- rgb(1,0.7,0,0.2)
col2b <- rgb(0,0.2,1,0.2)
col2_grey <- rgb(0,0,0,0.2)
col2_green <- rgb(0,0.7,0,0.2)

btsp <- 100 # Bootstrap sample

# Extract parameter range
# mle = mle_val, r95_rr = range_95_rr, r95_dec = range_95_decline, r95_imp = range_95_imp

# XX NEED TO UPDATE
# mle_r <- output_mle$mle[1]
# range_95 <- output_mle$r95_rr
# rr_est <- paste0(mle_r," (95% CI:",range_95[1],"-",range_95[2],")")

# Extrace values of interest
pred_interval <- output1$traj

long_dates <- output1$long_dates
ma_UK_cases_2 <- output1$mov_average



# - - -
# India cases

plot(all_india$date,all_india$cases_new,xlim=x_range,ylim=c(0,5e5),yaxs="i",ylab="cases",xlab="",main="India cases")
lines(c(india_red_list,india_red_list),c(0,1e7),col="grey",lty=2)

letter_ii <- 1
title(main=LETTERS[letter_ii],adj=0); letter_ii <- letter_ii+1

# - - -
# India local B.1.617

plot(data_india$date_time,data_india$B.1.617.2,xlim=x_range,ylim=c(0,1.05),yaxs="i",ylab="Proportion",xlab="",main="Proportion B.1.617.2 in India")
lines(all_india$date,ma_India_variant,lwd=2)

title(main=LETTERS[letter_ii],adj=0); letter_ii <- letter_ii+1

# - - -
# India imports
daily_india <- output1$daily_india

import_CI <- apply(intro_trace_tab,2,c.nume)

plot(long_dates,-1+0*daily_india,xlim=x_range,type="l",ylim=c(0,100),yaxs="i",ylab="cases",xlab="",main="Estimated imports/clusters in UK")
lines(c(india_red_list,india_red_list),c(0,1e7),col="grey",lty=2)
 
# Estimate for imports
imp_scale <- c.nume(thetatab[,"imp"])

ratio_rep <- sum(traveller_cases$Number_B_1_617_2)/sum(daily_india)

polygon(c(long_dates,rev(long_dates)),c(import_CI[2,],rev(import_CI[3,])),lty=0,col=col2c)
lines(long_dates,import_CI[1,],col="orange")

# PHE data
points(traveller_cases$Date,traveller_cases$Number_B_1_617_2)

# Estimates outbreak
lines(long_dates,pred_interval[1,],col="red")
polygon(c(long_dates,rev(long_dates)),c(pred_interval[2,],rev(pred_interval[3,])),lty=0,col=col2a)

text(x=tail(all_india$date,1),y=90,labels="Imported+cluster cases",col="red",cex=0.8)
text(x=tail(all_india$date,1),y=30,labels="Imported cases",col="orange",cex=0.8)

title(main=LETTERS[letter_ii],adj=0); letter_ii <- letter_ii+1

# - - -
# UK cases

# Extract UK fit data - copy of code from fit_R()
all_uk_fit <- all_uk %>% filter(date<date_uk_fit)
total_days_uk <- length(all_uk_fit$date)
ma_UK_cases_fit <- ma(all_uk_fit$cases_new,7)
t_max <- length(long_dates)

decline_95 <- c.nume(thetatab[,"decline"])
dt_95 <- c.nume(thetatab[,"dt_decline"])

ma_UK_cases_f1 <- decline_f(total_days_uk,t_max,decline_95[1],1,ma_UK_cases_fit)
ma_UK_cases_f2 <- decline_f(total_days_uk,t_max,decline_95[2],1,ma_UK_cases_fit)
ma_UK_cases_f3 <- decline_f(total_days_uk,t_max,decline_95[3],1,ma_UK_cases_fit)

cvector <- matrix(NA,nrow=btsp,ncol=length(ma_UK_cases_2))

for(ii in 1:btsp){
  ma_UK_cases_2_a <- ma_UK_cases_2; ma_UK_cases_2_a[1:3] <- 0
  cvector[ii,]= sapply(ma_UK_cases_2_a+pred_interval[1,],function(x){rnbinom(1,mu=x,size=1/sample(thetatab[,"rep_vol"],1))})
}
pred_interval_1 <- apply(cvector,2,c.nume)

# Plot data
plot(all_uk$date,all_uk$cases_new,xlim=x_range,ylim=c(0,6000),yaxs="i",ylab="cases",xlab="",main="UK cases")
lines(all_uk$date,ma_UK_cases,col="black",lwd=3)
lines(c(india_red_list,india_red_list),c(0,1e7),col="grey",lty=2)

# Plot decline
lines(long_dates,ma_UK_cases_f1,col="dark green",lty=1,lwd=2)
polygon(c(long_dates,rev(long_dates)),c(ma_UK_cases_f2,rev(ma_UK_cases_f3)),lty=0,col=col2_green)

# Plot B.1.617.2 and overall cases
lines(long_dates,pred_interval[1,],col="red")
#polygon(c(long_dates,rev(long_dates)),c(pred_interval[2,],rev(pred_interval[3,])),lty=0,col=col2a)

lines(long_dates,pred_interval_1[1,],col="blue")
polygon(c(long_dates,rev(long_dates)),c(pred_interval_1[2,],rev(pred_interval_1[3,])),lty=0,col=col2b)

title(main=LETTERS[letter_ii],adj=0); letter_ii <- letter_ii+1

doubling_time2 <- doubling_est2(pred_interval[1,long_dates>=as.Date("2021-05-01")])

# Number non_B117

cvector <- matrix(NA,nrow=btsp,ncol=length(ma_UK_cases_2_a))

for(ii in 1:btsp){
  ma_UK_cases_2_a <- as.numeric(ma_UK_cases_f1*data_fit$N/ma_UK_cases); ma_UK_cases_2_a[is.na(ma_UK_cases_2_a)] <- 0
  cvector[ii,]= sapply(ma_UK_cases_2_a,function(x){rnbinom(1,mu=x,size=1/sample(thetatab[,"rep_vol"],1))})
}
pred_interval_1 <- apply(cvector,2,c.nume)


plot(all_uk$date,-100+0*ma_UK_cases,xlim=x_range,ylim=c(0,2e3),yaxs="i",ylab="Number",xlab="",main="Non-B.1.617.2 sequences in UK")

date_cut2 <- length(data_proportion$sample_date) # Subtract less reliable data
date_cut1 <- date_cut2 - 5
polygon(c(data_proportion$sample_date[date_cut1],as.Date("2021-07-01"),as.Date("2021-07-01"),data_proportion$sample_date[date_cut1]),c(0,0,1e4,1e4),lty=0,col=col2_grey)


lines(c(india_red_list,india_red_list),c(0,1e5),col="grey",lty=2)
points(data_proportion$sample_date,data_proportion$N - data_proportion$B.1.617.2)

lines(long_dates,pred_interval_1[1,],col="dark green",lty=1,lwd=2) # Put nbinom uncertainty
polygon(c(long_dates,rev(long_dates)),c(pred_interval_1[2,],rev(pred_interval_1[3,])),lty=0,col=col2_green)

text(x=tail(data_proportion$sample_date,1),y=2e3*0.95,labels="Subject to delay",col=rgb(0.4,0.4,0.4),cex=0.8,adj=0)

title(main=LETTERS[letter_ii],adj=0); letter_ii <- letter_ii+1

# seq_pick <- data_proportion[match(long_dates,data_proportion$sample_date),]
# (seq_pick$N-seq_pick$B.1.617.2)-pred_interval_1[1,]

# - - -
# Number India-linked

cvector <- matrix(NA,nrow=btsp,ncol=length(ma_UK_cases_2_a))

for(ii in 1:btsp){
  ma_UK_cases_2_a <- as.numeric(pred_interval[1,]*data_fit$N/ma_UK_cases); ma_UK_cases_2_a[is.na(ma_UK_cases_2_a)] <- 0
  cvector[ii,]= sapply(ma_UK_cases_2_a,function(x){rnbinom(1,mu=x,size=1/sample(thetatab[,"rep_vol"],1))})
}
pred_interval_1 <- apply(cvector,2,c.nume)


plot(all_uk$date,-100+0*ma_UK_cases,xlim=x_range,ylim=c(0,400),yaxs="i",ylab="Number",xlab="",main="B.1.617.2 sequences in UK")

date_cut2 <- length(data_proportion$sample_date) # Subtract less reliable data
date_cut1 <- date_cut2 - 5
polygon(c(data_proportion$sample_date[date_cut1],as.Date("2021-07-01"),as.Date("2021-07-01"),data_proportion$sample_date[date_cut1]),c(0,0,1e3,1e3),lty=0,col=col2_grey)


lines(c(india_red_list,india_red_list),c(0,1e5),col="grey",lty=2)
points(data_proportion$sample_date,data_proportion$B.1.617.2)

lines(long_dates,pred_interval_1[1,],col="red",lty=1,lwd=2) # Put nbinom uncertainty
polygon(c(long_dates,rev(long_dates)),c(pred_interval_1[2,],rev(pred_interval_1[3,])),lty=0,col=col2a)

text(x=tail(data_proportion$sample_date,1),y=400*0.95,labels="Subject to delay",col=rgb(0.4,0.4,0.4),cex=0.8,adj=0)

title(main=LETTERS[letter_ii],adj=0); letter_ii <- letter_ii+1


# - - -
# Proportion India-linked
plot(all_uk$date,-1+0*ma_UK_cases,xlim=x_range,ylim=c(0,0.7),yaxs="i",ylab="Proportion",xlab="",main="Proportion B.1.617.2 in UK")
lines(c(india_red_list,india_red_list),c(0,1e5),col="grey",lty=2)

plot_CI(data_proportion$sample_date,data_proportion$B.1.617.2,data_proportion$N)
lines(long_dates,pred_interval[1,]/(ma_UK_cases_2+pred_interval[1,]),col="blue",lty=1)

#polygon(c(long_dates,rev(long_dates)),c(pred_interval[2,]/(pred_interval[2,]+ma_UK_cases_2),rev(pred_interval[3,]/(pred_interval[3,]+ma_UK_cases_2))),lty=0,col=col2b)
#lines(all_india$date,pred_interval[1,1:total_days]/(ma_UK_cases+pred_interval[1,1:total_days]),col="blue",lwd=2)


# polygon(c(long_dates,rev(long_dates)),c(pred_interval[2,]/ma_UK_cases_2,rev(pred_interval[3,]/ma_UK_cases_2)),lty=0,col=col2b)
# lines(all_india$date,pred_interval[1,1:total_days]/ma_UK_cases,col="blue",lwd=2)
# lines(long_dates,pred_interval[1,]/ma_UK_cases_2,col="blue",lty=2)



#text(x=tail(all_india$date,1),y=0.04,labels="R=0",col="blue",cex=0.8)
#text(x=tail(data_proportion$sample_date,1)-5,y=0.08,labels=paste0("R=",rr_est),col="blue",cex=0.8,adj=0)

title(main=LETTERS[letter_ii],adj=0); letter_ii <- letter_ii+1

# - - 
# Plot R estimates
kk <- iiM
if(kk==1){xmax <- 3.3; label_x <- c(1,2,3); label_list <- c("Traveller","Contact of traveller","Onwards") }
if(kk==2){xmax <- 2.3; label_x <- c(1,2); label_list <- c("Traveller","Non-traveller") }

plot(c(1:3),-1*c(1:3),xlim=c(0.7,xmax),ylim=c(0,7),xlab="",ylab="R",xaxt="n")
grid(ny = NULL, nx=NA, col = "lightgray")

#axis(1, at=c(1,2,3,4), labels=c("Traveller","Contact of traveller","Onwards","SPI-M"))
axis(1, at=label_x, labels=label_list)

xx <- 4
#lines(c(xx,xx),c(0.8,1.1),lwd=4,lend=1) # SPI-M consensus

store_val <- NULL

for(ii in 1:3){
  if(ii==1){range_ii <- c.nume.50(thetatab[,1])}
  if(ii==2){range_ii <- c.nume.50(thetatab[,1]*thetatab[,2])}
  if(ii==3 & kk==1){range_ii <- c.nume.50(thetatab[,1]*thetatab[,2]*thetatab[,3])}
  
  
  if(ii==1){range_A <- c.text(thetatab[,1])}
  if(ii==2){range_A <- c.text(thetatab[,1]*thetatab[,2])}
  store_val <- rbind(store_val,c(ii,range_A))
  
  lines(c(ii,ii),c(range_ii[2],range_ii[5]),col=col2b)
  lines(c(ii,ii),c(range_ii[3],range_ii[4]),lwd=2.5,col="light blue")
  points(ii,range_ii[1],pch=19,cex=1,col="blue")
  
  if(ii==1){points(ii,r_phe_report[["travel"]])}
  if(ii>1){points(ii,r_phe_report[["non_travel"]])}
  
}
 title(main=LETTERS[letter_ii],adj=0); letter_ii <- letter_ii+1

# plot(rr_store$X1,rr_store$X2,xlim=c(min(rr_range),1.9),ylim=c(-220,-130),xlab="R",ylab="lik",main=paste0("R=", rr_est,", k=",kk_pick,""))
# lines(xx_list,preds$fit)
# lines(c(mle_val,mle_val),c(-1e3,0),lty=2)
# lines(c(range_95[1],range_95[1]),c(-1e3,0),lty=3)
# lines(c(range_95[2],range_95[2]),c(-1e3,0),lty=3)
# 
# title(main=LETTERS[letter_ii],adj=0); letter_ii <- letter_ii+1

# Output plots
dev.copy(png,paste0("outputs/plot_",iiM,".png"),units="cm",width=30,height=15,res=200)
dev.off()

# Plot posteriors ---------------------------------------------------------

plot_post <- function(){
  par(mfcol=c(3,2),mar=c(3,3,1,1),mgp=c(2,0.6,0),las=0)
  
  
  # Plot likelihood
  plot(sim_likOut,type="l")
  
  hist(thetatab[,1],main="R traveller")
  hist(thetatab[,1]*thetatab[,2],main="R second")
  hist(thetatab[,1]*thetatab[,3],main="R onwards")
  hist(thetatab[,4],main="decline")
  #hist(thetatab[,5]-1,main="dt_decline")
  hist(thetatab[,6],main="imports")
  
  # Output plots
  dev.copy(png,paste0("outputs/posterior_",iiM,".png"),units="cm",width=20,height=15,res=200)
  dev.off()

}

# Plot R fits on same plot  ---------------------------------------------------------

compare_R_fits <- function(){
  
  par(mfcol=c(1,2),mar=c(3,3,1,1),mgp=c(2,0.6,0),las=0)
  letter_ii <- 1
  
  for(kk in 2:1){
    
    iiM <- kk
    source("R/load_posterior.R",local=TRUE)
  
    if(kk==1){xmax <- 3.3; label_x <- c(1,2,3); label_list <- c("Traveller","Contact of traveller","Onwards") }
    if(kk==2){xmax <- 2.3; label_x <- c(1,2); label_list <- c("Traveller","Non-traveller") }
    
    plot(c(1:3),-1*c(1:3),xlim=c(0.7,xmax),ylim=c(0,7),xlab="",ylab="R",xaxt="n")
    grid(ny = NULL, nx=NA, col = "lightgray")
    
    #axis(1, at=c(1,2,3,4), labels=c("Traveller","Contact of traveller","Onwards","SPI-M"))
    axis(1, at=label_x, labels=label_list)
    
    xx <- 4
    #lines(c(xx,xx),c(0.8,1.1),lwd=4,lend=1) # SPI-M consensus
    
    
    for(ii in 1:3){
      if(ii==1){range_ii <- c.nume.50(thetatab[,1])}
      if(ii==2){range_ii <- c.nume.50(thetatab[,1]*thetatab[,2])}
      if(ii==3 & kk==1){range_ii <- c.nume.50(thetatab[,1]*thetatab[,2]*thetatab[,3])}
      
      lines(c(ii,ii),c(range_ii[2],range_ii[5]),col=col2b)
      lines(c(ii,ii),c(range_ii[3],range_ii[4]),lwd=2.5,col="light blue")
      points(ii,range_ii[1],pch=19,cex=1,col="blue")
      
      if(ii==1){points(ii,r_phe_report[["travel"]])}
      if(ii>1){points(ii,r_phe_report[["non_travel"]])}
      
    }
    
    title(main=LETTERS[letter_ii],adj=0); letter_ii <- letter_ii+1
    
  }
  
  dev.copy(png,paste0("outputs/compare_R.png"),units="cm",width=25,height=10,res=200)
  dev.off()
  
}

# Output R estimates
print(store_val)




