#library(robustbase)

files <- list.files(path=".", pattern="^CE.*\\.txt.table$");
mydata <- lapply(files, read.table, header=FALSE, col.names=c("time","voltage"));
files <- sub(".txt.table","",files);
names(mydata) <- files;
write.table(t(c("file","ChargeDensityCE")), file="outputChargeDensityCE.txt", append=TRUE, col.names=F, row.names=F);



trashfornullmessages <- lapply(files, function(x) {
	message(x);
	if(!file.exists(paste(x, ".png", sep=""))){
	len<-length(mydata[[x]]$voltage)
	endingvoltage <- mean(mydata[[x]]$voltage[(len*0.6):len])
#	fit <- lm(voltage ~ poly(time,2),data=mydata[[x]])#start=list(A=startingvoltage,B=0,C=0);
	#fitR <- nlrob(voltage ~  A + B*time + C*time^2, trace=T,data=mydata[[x]], start=list(A=startingvoltage,B=0,C=0);#A=coef(fit)[1],B=coef(fit)[2],C=coef(fit)[3]),)#,C=coef(fit)[3], D=coef(fit)[4], E=coef(fit)[5]));
	#fitR <- lmrob(voltage ~ A, trace=T, data=mydata[[x]]);
	message(endingvoltage)
	voltage2 <- mydata[[x]]$voltage - endingvoltage#, newdata = data.frame(mydata[[x]]$time))
	current <- voltage2/50
	charge <- cumsum(current)*(mydata[[x]]$time[2]-mydata[[x]]$time[1])

	charge=charge-charge[match(0,mydata[[x]]$time)]
#	print(quantile(sum))
	totalcharge=quantile(charge[charge>0],0.75)
	totalchargedensity=totalcharge/0.09
        outputChargeDensityCE <- t(c(x, abs(totalchargedensity)));
	write.table(outputChargeDensityCE, file="outputChargeDensityCE.txt", append=TRUE, col.names=F, row.names=F, quote=F);


#	png(paste(x, "-test.png", sep=""), width=1280, heigh=800)
#	par(mar=c(5,4,4,5)+.1)
#	plot(mydata[[x]]$time, current, type="l", ylab="Current (A)", xlab="Time (s)", xlim=c(-1e-6,1e-5), main=paste(x,"CE"))
#	par(new=TRUE)
#	plot(mydata[[x]]$time,charge, type="l", col="red", xaxt="n",yaxt="n",xlab="",ylab="")
#	#lines(mydata[[x]]$time,charge,col="red")
#	abline(h=0,col="red")
#	abline(h=totalcharge,col="red")
#	axis(4,col.ticks="red",col.axis="red", col="red")
#	mtext("Collected Charge (C)",side=4,line=3,col="red")
##	text(tail(mydata[[x]]$time,1)*0.9,totalcharge*0.95,labels=paste(abs(signif(totalcharge,4)), "C"),cex=2,col="red")
##	graphics.off()


	png(paste(x, ".png", sep=""), width=1280, heigh=800)
	par(mar=c(5,4,4,5)+.1)
	plot(mydata[[x]],type="l", ylab="Voltage (V)", xlab="Time (s)", main=paste(x,"CE"))
	abline(h=endingvoltage, col="green")
	par(new=TRUE)
	plot(mydata[[x]]$time,charge/0.09, type="l", col="red", xaxt="n",yaxt="n",xlab="",ylab="")
	abline(h=0,col="red")
	abline(h=totalchargedensity,col="red")
	axis(4,col.ticks="red",col.axis="red", col="red")
	mtext("Collected Charge Density (C/cm2)",side=4,line=3,col="red")
	text(tail(mydata[[x]]$time,1)*0.9,totalchargedensity*0.95,labels=paste(abs(signif(totalchargedensity,4)), "C/cm2"),cex=2,col="red")
	graphics.off()
}})

