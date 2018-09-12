# Copyright (C) 2018 Ilario Gelmetti <iochesonome@gmail.com>
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

ceWithLimits <- function(cedir="ce", dcdir=".")
{
print("CE with limits: PLOTTING")

require(minpack.lm)
require(robustbase)

files <- list.files(path=cedir, pattern="^CE.*\\.txt.table$");
mydata <- lapply(file.path(cedir,files), read.table, header=FALSE, col.names=c("time","voltage"));
files <- sub(".txt.table","",files);
names(mydata) <- files;
write.table(t(c("Voc","RCtime")), file=file.path(cedir,"outputRCtime.txt"), append=FALSE, col.names=F, row.names=F);

impedance = 50

DCcapacitance <- read.table(file.path(dcdir,"outputDCcapacitance.txt"), header=TRUE);
DCcapacitance <- DCcapacitance[with(DCcapacitance, order(Voc)), ]

tryCatch({
	expfitDC <- nls(capacitance ~ exp(B) + exp(C)*exp(D)*exp(exp(D)*Voc), start=list(B=log(max(1e-8,min(DCcapacitance$capacitance))),C=log(1e-10),D=2), data=DCcapacitance)
}, error=function(e) print("Failed restricted to positive gamma fit - first"))
tryCatch({
	expfitDC <- nlsLM(capacitance ~ exp(B) + exp(C)*exp(D)*exp(exp(D)*Voc), start=list(B=log(max(1e-8,min(DCcapacitance$capacitance))),C=log(1e-10),D=2), data=DCcapacitance)
}, error=function(e) print("Failed restricted to positive gamma fit - second"))
tryCatch({
	expfitDC <- nlrob(capacitance ~ exp(B) + exp(C)*exp(D)*exp(exp(D)*Voc), start=list(B=coef(expfitDC)[[1]],C=coef(expfitDC)[[2]],D=coef(expfitDC)[[3]]), data=DCcapacitance)
}, error=function(e) print("Failed restricted to positive gamma robust fit"))

	#importante che la variabile in new abbia lo stesso nome di quella fittata
	new <- data.frame(Voc = 0)
	geometrical_capacitance <- 0.09*predict(expfitDC, new)

#RC time is R[ohm]*C[F/cm2]*area[cm2]
discharge_func = function(t, C) exp(-t/(impedance*C))
discharge_func_fixedC = function(t) exp(-t/(impedance*geometrical_capacitance))

lapply(files, function(x) {
	message(x);
	lo = loess(mydata[[x]]$voltage~mydata[[x]]$time, span=0.01);
	loess_voltage = predict(lo)
	maxV = max(loess_voltage)
#maxV=mean(sort(mydata[[x]]$voltage, decreasing = TRUE)[1:10])

	b<-strsplit(x, "_")
	c<-unlist(b)
	c2 <- c[grepl("mV",c)]
	Voc_fromfilename<-as.numeric(sub("mV.*", "", c2))

	#importante che la variabile in new abbia lo stesso nome di quella fittata
	new <- data.frame(Voc = Voc_fromfilename)

	local_capacitance <- 0.09*predict(expfitDC, new)
	t_max = mydata[[x]]$time[which.max(loess_voltage)]

	discharge_profile = maxV * discharge_func(mydata[[x]]$time, local_capacitance)
	discharge_profile_fixedC = maxV * discharge_func_fixedC(mydata[[x]]$time)

	png(file.path(cedir,paste(x, "-limits.png", sep="")), width=800, height=600)
	par(mar=c(5,4,4,5)+.1)
	plot(mydata[[x]],type="l", ylab="Voltage (V)", xlab="Time (s)", main=paste(x,"CE"), xlim=c(-1e-6, 8e-6))
	lines(mydata[[x]]$time + t_max, discharge_profile, col="red")
	lines(mydata[[x]]$time + t_max, discharge_profile_fixedC, col="blue")
#abline(h=maxV)
#lines(mydata[[x]]$time, predict(lo), col="green", lwd=2)



# I want the fitting data starting from the first point after the FWHM
	peak_times = mydata[[x]]$time[mydata[[x]]$voltage > 0.5 * max(mydata[[x]]$voltage)];
	timeEndFWHM = tail(peak_times, n=1);
# sometimes the noise is way bigger than the signal peak (at low light intensity), so I can help the start decay time to be after the noise putting it after the minimum value (on our equipment corresponds to the end of the noise)
	negative_peak_time = mydata[[x]]$time[which.min(mydata[[x]]$voltage)]
	time_start_decay = max(timeEndFWHM, negative_peak_time)

	decay <- mydata[[x]][mydata[[x]]$time > time_start_decay,]
	
tryCatch({
	expfitCE <- nlsLM(voltage~ C*exp(-time/D), start=list(C=max(decay$voltage),D=0.01*tail(decay$time, n=1)), data=decay)
	tryCatch({
		expfitCE <- nlrob(voltage~ C*exp(-time/D), start=list(C=coef(expfitCE)["C"],D=coef(expfitCE)["D"]), data=decay)
	}, error=function(e) cat("Failed monoexponential robust fit", e$message, "\n"))

	lines(decay$time, predict(expfitCE), lwd=2, col="magenta")

	b<-strsplit(x, "_")
	c<-unlist(b)
	c2 <- c[grepl("mV",c)]
	d<-as.numeric(sub("mV.*", "", c2))
	outputMonoexpCE <- t(c(d, coef(expfitCE)["D"]));

}, error=function(e) cat("Failed monoexponential fit", e$message, "\n"))
	

	legend(x="topright",inset=0.1,c("Charge Extraction", paste("RC time local cap DC",signif(impedance*geometrical_capacitance,3)), paste("RC time geom cap DC",signif(impedance*local_capacitance,3)), paste("CE monoexp fit",signif(coef(expfitCE)["D"],3))), lty=1, lwd=6, cex=1.5, col=c("black","red","blue","magenta"))

	graphics.off()

	outputRCtime <- t(c(Voc_fromfilename, impedance*local_capacitance));
	write.table(outputRCtime, file=file.path(cedir,"outputRCtime.txt"), append=TRUE, col.names=F, row.names=F, quote=F);
	write.table(outputMonoexpCE, file=file.path(cedir,"outputMonoexpCE.txt"), append=TRUE, col.names=F, row.names=F, quote=F);
})
}
