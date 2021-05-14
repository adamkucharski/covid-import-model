# Plot data ---------------------------------------------------------------

rr_store <- read_csv(paste0("outputs/fit",kk_pick,".csv"))

par(mfcol=c(3,2),mar=c(3,3,1,1),mgp=c(2,0.6,0),las=0)

x_range <- c(as.Date("2021-03-01"),max(all_uk$date)+25)
col2a=rgb(1,0.2,0,0.2)
col2c=rgb(1,0.7,0,0.2)
col2b=rgb(0,0.2,1,0.2)


# Extract parameter range
# mle = mle_val, r95_rr = range_95_rr, r95_dec = range_95_decline, r95_imp = range_95_imp

mle_val <- output_mle$mle[1]
range_95 <- output_mle$r95_1
rr_est <- paste0(mle_val," (95% CI:",range_95[1],"-",range_95[2],")")

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

pred_interval <- output1$traj
long_dates <- output1$long_dates
ma_UK_cases_2 <- output1$mov_average

# - - -
# India imports
daily_india <- output1$daily_india

plot(all_india$date,daily_india,xlim=x_range,type="l",ylim=c(0,100),yaxs="i",ylab="cases",xlab="",main="Estimated imported cases/clusters into UK")
lines(c(india_red_list,india_red_list),c(0,1e7),col="grey",lty=2)

points(traveller_cases$Date,traveller_cases$Number_B_1_617_2)

lines(long_dates,pred_interval[1,],col="red")
polygon(c(long_dates,rev(long_dates)),c(pred_interval[2,],rev(pred_interval[3,])),lty=0,col=col2a)

text(x=tail(all_india$date,1),y=90,labels="Imported+cluster cases",col="red",cex=0.8)
text(x=tail(all_india$date,1),y=50,labels="Imported cases",col="black",cex=0.8)

title(main=LETTERS[letter_ii],adj=0); letter_ii <- letter_ii+1

# - - -
# UK cases
plot(all_uk$date,all_uk$cases_new,xlim=x_range,ylim=c(0,1e4),yaxs="i",ylab="cases",xlab="",main="UK cases")
lines(all_uk$date,ma_UK_cases,col="black",lwd=3)
lines(long_dates,ma_UK_cases_2,col="black",lty=2)
lines(c(india_red_list,india_red_list),c(0,1e7),col="grey",lty=2)

lines(long_dates,pred_interval[1,],col="red")
polygon(c(long_dates,rev(long_dates)),c(pred_interval[2,],rev(pred_interval[3,])),lty=0,col=col2a)

lines(long_dates,ma_UK_cases_2+pred_interval[1,],col="orange")
polygon(c(long_dates,rev(long_dates)),c(ma_UK_cases_2+pred_interval[2,],rev(ma_UK_cases_2+pred_interval[3,])),lty=0,col=col2c)

title(main=LETTERS[letter_ii],adj=0); letter_ii <- letter_ii+1

# - - -
# Proportion India-linked
plot(all_uk$date,-1+0*ma_UK_cases,xlim=x_range,ylim=c(0,0.5),yaxs="i",ylab="Proportion",xlab="",main="Proportion B.1.617.2 in UK")
lines(c(india_red_list,india_red_list),c(0,1e5),col="grey",lty=2)

polygon(c(long_dates,rev(long_dates)),c(pred_interval[2,]/(pred_interval[2,]+ma_UK_cases_2),rev(pred_interval[3,]/(pred_interval[3,]+ma_UK_cases_2))),lty=0,col=col2b)
#lines(all_india$date,pred_interval[1,1:total_days]/(ma_UK_cases+pred_interval[1,1:total_days]),col="blue",lwd=2)
lines(long_dates,pred_interval[1,]/(ma_UK_cases_2+pred_interval[1,]),col="blue",lty=1)

# polygon(c(long_dates,rev(long_dates)),c(pred_interval[2,]/ma_UK_cases_2,rev(pred_interval[3,]/ma_UK_cases_2)),lty=0,col=col2b)
# lines(all_india$date,pred_interval[1,1:total_days]/ma_UK_cases,col="blue",lwd=2)
# lines(long_dates,pred_interval[1,]/ma_UK_cases_2,col="blue",lty=2)

plot_CI(data_proportion$sample_date,data_proportion$B.1.617.2,data_proportion$N)

#text(x=tail(all_india$date,1),y=0.04,labels="R=0",col="blue",cex=0.8)
text(x=tail(all_india$date,1),y=0.08,labels=paste0("R=",mle_val),col="blue",cex=0.8,adj=0)

title(main=LETTERS[letter_ii],adj=0); letter_ii <- letter_ii+1

# Plot R

# plot(rr_store$X1,rr_store$X2,xlim=c(min(rr_range),1.9),ylim=c(-220,-130),xlab="R",ylab="lik",main=paste0("R=", rr_est,", k=",kk_pick,""))
# lines(xx_list,preds$fit)
# lines(c(mle_val,mle_val),c(-1e3,0),lty=2)
# lines(c(range_95[1],range_95[1]),c(-1e3,0),lty=3)
# lines(c(range_95[2],range_95[2]),c(-1e3,0),lty=3)
# 
# title(main=LETTERS[letter_ii],adj=0); letter_ii <- letter_ii+1


# Output plots
dev.copy(png,paste0("outputs/plot_k",kk_pick,"_u",under_factor,"_d",daily_decline,".png"),units="cm",width=20,height=15,res=200)
dev.off()



# Other plots
# 
# par(mfcol=c(1,1),mar=c(3,3,1,1),mgp=c(2,0.6,0),las=0)
# plot(all_india$date,-1+0*daily_india,xlim=x_range,type="l",ylim=c(0,40),yaxs="i",ylab="cases",xlab="",main="B.1.617.2 traveller cases")
# lines(c(india_red_list,india_red_list),c(0,1e7),col="grey",lty=2)
# 
# points(traveller_cases$Date,traveller_cases$Number_B_1_617_2,pch=19)
# lines(traveller_cases$Date,traveller_cases$Number_B_1_617_2)
# 
# dev.copy(png,paste0("outputs/traveller_cases.png"),units="cm",width=15,height=10,res=200)
# dev.off()


