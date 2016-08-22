# Copyright (C) 2015-2016 Ilario Gelmetti <iochesonome@gmail.com>
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.

#library(robustbase)

tpc <- function(tpcdir="tpc")
{
files <- list.files(path=tpcdir, pattern="^TPC.*\\.txt.table$");
mydata <- lapply(file.path(tpcdir,files), read.table, header=FALSE, col.names=c("time","voltage"));
files <- sub(".txt.table","",files);
names(mydata) <- files;
write.table(t(c("file","ChargeDensityTPC")), file=file.path(tpcdir,"outputChargeDensityTPC.txt"), append=FALSE, col.names=F, row.names=F);

trashfornullmessages <- lapply(files, function(x) {
	message(x);
#	if(!file.exists(paste(x, ".png", sep=""))){
	len<-length(mydata[[x]]$voltage)
	endingvoltage <- mean(mydata[[x]]$voltage[(len*0.9):len])
	message(endingvoltage)
	voltage2 <- mydata[[x]]$voltage - endingvoltage
	current <- voltage2/50
	charge <- cumsum(current)*(mydata[[x]]$time[2]-mydata[[x]]$time[1])

#	charge=charge-charge[which.min(mydata[[x]]$voltage)]
	charge=charge-charge[match(0,mydata[[x]]$time)]
	totalcharge=quantile(charge,0.30)
	totalchargedensity=totalcharge/0.09
        outputChargeDensityTPC <- t(c(x, abs(totalchargedensity)));
	write.table(outputChargeDensityTPC, file=file.path(tpcdir,"outputChargeDensityTPC.txt"), append=TRUE, col.names=F, row.names=F, quote=F);


	png(file.path(tpcdir,paste(x, ".png", sep="")), width=1280, heigh=800)
	par(mar=c(5,4,4,5)+.1)
	plot(mydata[[x]],type="l", ylab="Voltage (V)", xlab="Time (s)", main=paste(x, "TPC"))
#	lines(mydata[[x]]$time, predict(fitR), col="green")
	abline(h=endingvoltage, col="green")
	par(new=TRUE)
	plot(mydata[[x]]$time,charge/0.09, type="l", col="red", xaxt="n",yaxt="n",xlab="",ylab="")
	abline(h=0,col="red")
	abline(h=totalchargedensity,col="red")
	axis(4,col.ticks="red",col.axis="red", col="red")
	mtext("Collected Charge Density (C/cm2)",side=4,line=3,col="red")
	text(tail(mydata[[x]]$time,1)*0.9,totalchargedensity*0.95,labels=paste(abs(signif(totalchargedensity,4)), "C/cm2"),cex=2,col="red")
	graphics.off()
#}
})
}
