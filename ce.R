ce <- function(cedir="ce")
{
library(sfsmisc)

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


	png(file.path(cedir,paste(x, ".png", sep="")), width=image_width, height=image_height)
	op <- par(mar = c(5,7,4,8.5) + 0.1) ## default is c(5,4,4,2) + 0.1

	xlim=c(0,0.8e-5)
	plot(mydata[[x]],type="l", ylab="", xlab="", xlim=xlim, xaxt="n", yaxt="n")
	title(ylab="Voltage (V)", cex.lab=2, line=4)
	title(xlab="Time (s)", cex.lab=2, line=3.5)
	mtext(bquote("Collected Charge Density (C/cm"^"2"*")"), cex=2, side=4,line=7,col="red")
	eaxis(side=1, cex.axis=1.5)
	eaxis(side=2, cex.axis=1.5)
	
	lines(mydata[[x]]$time, baseline, col="green")
	par(new=TRUE)
	plot(mydata[[x]]$time,charge/0.09, type="l", col="red", xaxt="n",yaxt="n",xlab="",ylab="", xlim=xlim)
	abline(h=0,col="red")
	#abline(h=totalchargedensity,col="red")
	eaxis(4,col.ticks="red",col.axis="red", col="red", cex.axis=1.5)
	text(xlim[2]*0.75,totalchargedensity*0.9,labels=bquote(.(signif(totalchargedensity,3))~"C/cm"^"2"),cex=2,col="red")
	
	#par(new=TRUE)
	#plot(mydata[[x]]$time,chargezero/0.09, type="l", col="orange", xaxt="n",yaxt="n",xlab="",ylab="")
	
	graphics.off()
#reset the plotting margins
par(op)


write.table(outputChargeDensityCE, file=file.path(cedir,"outputChargeDensityCE.txt"), append=TRUE, col.names=F, row.names=F, quote=F);
})
}
