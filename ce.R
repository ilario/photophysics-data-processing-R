#library(robustbase)
ce <- function(cedir="ce")
{
require(minpack.lm)
print("CE: INTEGRATING")
files <- list.files(path=cedir, pattern="^CE.*\\.txt.table$");
mydata <- lapply(file.path(cedir,files), read.table, header=FALSE, col.names=c("time","voltage"));
files <- sub(".txt.table","",files);
names(mydata) <- files;
write.table(t(c("Voc","ChargeDensityCE")), file=file.path(cedir,"outputChargeDensityCE.txt"), append=FALSE, col.names=F, row.names=F);


trashfornullmessages <- lapply(files, function(x) {
	message(x);
#	if(!file.exists(paste(x, ".png", sep=""))){
	len<-length(mydata[[x]]$voltage)
	endingvoltage <- mean(mydata[[x]]$voltage[(len*0.9):len]) #could be just zero
	message(endingvoltage)
	voltage2 <- mydata[[x]]$voltage - endingvoltage#, newdata = data.frame(mydata[[x]]$time))
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
	c<-unlist(b)[length(b[[1]])]
	d<-as.numeric(gsub("mV", "", c))
        outputChargeDensityCE <- t(c(d, totalchargedensity));
	write.table(outputChargeDensityCE, file=file.path(cedir,"outputChargeDensityCE.txt"), append=TRUE, col.names=F, row.names=F, quote=F);

	decay <- mydata[[x]][mydata[[x]]$time>0,]


	png(file.path(cedir,paste(x, ".png", sep="")), width=1280, heigh=800)
	par(mar=c(5,4,4,5)+.1)
	plot(mydata[[x]],type="l", ylab="Voltage (V)", xlab="Time (s)", main=paste(x,"CE"))
tryCatch({
expfit <- nlsLM(voltage~ C*exp(D*time), start=list(C=0.5*max(decay$voltage),D=-0.1*tail(decay$time, n=1)), data=decay)
	lines(decay$time, predict(expfit), lwd=2, col="blue")
}, error=function(e) print("Failed monoexponential fit"))
	abline(h=endingvoltage, col="green")
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
})
}
