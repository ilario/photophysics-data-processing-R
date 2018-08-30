ce <- function(cedir="ce")
{
require(minpack.lm)
library(robustbase)
print("CE: INTEGRATING")
files <- list.files(path=cedir, pattern="^CE.*\\.txt.table$");
mydata <- lapply(file.path(cedir,files), read.table, header=FALSE, col.names=c("time","voltage"));
files <- sub(".txt.table","",files);
names(mydata) <- files;
write.table(t(c("Voc","ChargeDensityCE")), file=file.path(cedir,"outputChargeDensityCE.txt"), append=FALSE, col.names=F, row.names=F);
write.table(t(c("Voc","CEmonoexpTime")), file=file.path(cedir,"outputMonoexpCE.txt"), append=FALSE, col.names=F, row.names=F);

trashfornullmessages <- lapply(files, function(x) {
	message(x);
	len<-length(mydata[[x]]$voltage)
#	endingvoltage <- mean(mydata[[x]]$voltage[(len*0.9):len]) #could be just zero, better a linear baseline from the beginning value to the final one
#	voltage2 <- mydata[[x]]$voltage - endingvoltage#, newdata = data.frame(mydata[[x]]$time))
	startVoltage <- mean(mydata[[x]]$voltage[1:600])
	endVoltage <- mean(mydata[[x]]$voltage[(len-600):len])
	baseline <- seq(startVoltage, endVoltage, length.out=len)
	voltage2 <- mydata[[x]]$voltage - baseline
	
	current <- voltage2/50
	charge <- cumsum(current)*(mydata[[x]]$time[2]-mydata[[x]]$time[1])

	charge=charge-charge[match(0,mydata[[x]]$time)]
	totalcharge=mean(charge[round(length(charge)*0.9):round(length(charge)*0.95)])
#	totalcharge=quantile(charge[charge>0],0.75)
	totalchargedensity=totalcharge/0.09

	voltagezero <- mydata[[x]]$voltage
	currentzero <- voltagezero/50
	chargezero <- cumsum(currentzero)*(mydata[[x]]$time[2]-mydata[[x]]$time[1])
	chargezero=chargezero-chargezero[match(0,mydata[[x]]$time)]

	b<-strsplit(x, "_")
	c<-unlist(b)
	c2 <- c[grepl("mV",c)]
	d<-as.numeric(sub("mV.*", "", c2))
        outputChargeDensityCE <- t(c(d, totalchargedensity));
	write.table(outputChargeDensityCE, file=file.path(cedir,"outputChargeDensityCE.txt"), append=TRUE, col.names=F, row.names=F, quote=F);

# I want the fitting data starting from the first point after the FWHM
	#decay <- mydata[[x]][mydata[[x]]$time>0,]
	peak_times = mydata[[x]]$time[mydata[[x]]$voltage > 0.5 * max(mydata[[x]]$voltage)];
	timeEndFWHM = tail(peak_times, n=1);
# sometimes the noise is way bigger than the signal peak (at low light intensity), so I can help the start decay time to be after the noise putting it after the minimum value (on our equipment corresponds to the end of the noise)
	negative_peak_time = mydata[[x]]$time[which.min(mydata[[x]]$voltage)]
	time_start_decay = max(timeEndFWHM, negative_peak_time)

	decay <- mydata[[x]][mydata[[x]]$time > time_start_decay,]

	png(file.path(cedir,paste(x, ".png", sep="")), width=1280, heigh=800)
	par(mar=c(5,4,4,5)+.1)
	plot(mydata[[x]],type="l", ylab="Voltage (V)", xlab="Time (s)", main=paste(x,"CE"), xlim=c(0,1e-5))

tryCatch({
	expfit <- nlsLM(voltage~ C*exp(-time/D), start=list(C=max(decay$voltage),D=0.01*tail(decay$time, n=1)), data=decay)
	tryCatch({
		expfit <- nlrob(voltage~ C*exp(-time/D), start=list(C=coef(expfit)["C"],D=coef(expfit)["D"]), data=decay)
	}, error=function(e) print("Failed monoexponential robust fit"))

	lines(decay$time, predict(expfit), lwd=2, col="magenta")

	outputMonoexpCE <- t(c(d, coef(expfit)["D"]));

}, error=function(e) print("Failed monoexponential fit"))

	lines(mydata[[x]]$time, baseline, col="green")
	par(new=TRUE)
	plot(mydata[[x]]$time,charge/0.09, type="l", col="red", xaxt="n",yaxt="n",xlab="",ylab="")
	abline(h=0,col="red")
	abline(h=totalchargedensity,col="red")
	axis(4,col.ticks="red",col.axis="red", col="red")
	mtext("Collected Charge Density (C/cm2)",side=4,line=3,col="red")
	text(tail(mydata[[x]]$time,1)*0.9,totalchargedensity*0.95,labels=paste(signif(totalchargedensity,4), "C/cm2"),cex=2,col="red")
	
	par(new=TRUE)
	plot(mydata[[x]]$time,chargezero/0.09, type="l", col="orange", xaxt="n",yaxt="n",xlab="",ylab="")
	
	graphics.off()
#}

write.table(outputChargeDensityCE, file=file.path(cedir,"outputChargeDensityCE.txt"), append=TRUE, col.names=F, row.names=F, quote=F);
write.table(outputMonoexpCE, file=file.path(cedir,"outputMonoexpCE.txt"), append=TRUE, col.names=F, row.names=F, quote=F);

})
}
