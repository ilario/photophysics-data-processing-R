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

dcFromOutputToGraph <- function(tpvdir="tpv", tpcdir="tpc")
{
print("DC: PLOTTING")
require(minpack.lm)
require(robustbase)

write.table(t(c("Voc","capacitance")), file="outputDCcapacitance.txt", append=FALSE, col.names=F, row.names=F);
write.table(t(c("Voc","ChargeDensityDC")), file="outputDCcharge.txt", append=FALSE, col.names=F, row.names=F);
write.table(t(c("Voc","ChargeDensityDC")), file="outputDCcharge-nogeom.txt", append=FALSE, col.names=F, row.names=F);

a <- read.table(file.path(tpcdir, "outputChargeDensityTPC.txt"), header=T)
#chargeDark <- mean(a[grep("dark", a$file, ignore.case=T),]$ChargeDensityTPC)
#chargeSun <- mean(a[grep("sun", a$file, ignore.case=T),]$ChargeDensityTPC)
#chargeArray = seq(chargeDark, chargeSun, length.out=len)
# in case TPC in dark and in sun are different, the choice of what to use is arbitrary, I would use the third quartile of all the TPC measurements
charge <- quantile(a$ChargeDensityTPC, 0.75)

if(file.exists(file.path(tpvdir, "outputDeltaVmixed.txt"))){
	print("DC: using DeltaV from mixed monoexp and biexp")
	b <- read.table(file.path(tpvdir, "outputDeltaVmixed.txt"), header=T)
	names(b) <- c("file", "Voc", "deltaV")
	b <- b[with(b, order(b$Voc)), ]
}else{
	b <- read.table(file.path(tpvdir, "outputDeltaVfirstPoints.txt"), header=T)

#	b <- read.table(file.path(tpvdir, "outputDeltaV.txt"), header=T)
#	bLoess <- read.table(file.path(tpvdir, "outputDeltaVloess.txt"), header=T)
#	bMonoexp <- read.table(file.path(tpvdir, "outputDeltaVmonoexp.txt"), header=T)
#	names(b) <- c("file", "Voc", "deltaV")
#	names(bLoess) <- c("file", "Voc", "deltaV")
#	names(bMonoexp) <- c("file", "Voc", "deltaV")
#	b <- b[with(b, order(b$Voc)), ]
#	bLoess <- bLoess[with(bLoess, order(bLoess$Voc)), ]
#	bMonoexp <- bMonoexp[with(bMonoexp, order(bMonoexp$Voc)), ]
#
#	#remove lines where monoexp fit failed
#	matchIndexMonoexp <- match(bMonoexp$file, b$file)
#	b <- b[matchIndexMonoexp,]
#	bLoess <- bLoess[matchIndexMonoexp,]
#
#	#element wise maximum
#	bMonoexpLoess <- bMonoexp
#	bMonoexpLoess$deltaV <- rowMeans(cbind(bMonoexp$deltaV, bLoess$deltaV))
#	# uses the maximum between Monoexp and Loess close to 1 sun and the plain deltaV close to dark, with a linear mixing between the two
#	lenMatch <- length(matchIndexMonoexp)
#	b$deltaV <- (seq(lenMatch,1)*b$deltaV + seq(1,lenMatch)*bMonoexpLoess$deltaV) / (lenMatch+1)
}

write.table(b, file=file.path(tpvdir, "outputDeltaVprocessedForDC.txt"), append=FALSE, row.names=FALSE)

len <- length(b$deltaV)

chargeArray = rep(charge, len)

getExpFit <- function(){
	capacitance <<- chargeArray/b$deltaV
	directory <<- tail(strsplit(getwd(), "/")[[1]], n=1)

	outputDCcapacitance <<- data.frame(b$Voc, capacitance);
	names(outputDCcapacitance) <<- c("Voc","capacitance")

	tryCatch({
		expfit <<- nls(capacitance ~ exp(B) + exp(C)*exp(D)*exp(exp(D)*Voc), start=list(B=log(max(1e-8,min(outputDCcapacitance$capacitance))),C=log(1e-10),D=2), data=outputDCcapacitance)
	}, error=function(e) print("Failed restricted to positive gamma fit - first"))
	tryCatch({
		expfit <<- nlsLM(capacitance ~ exp(B) + exp(C)*exp(D)*exp(exp(D)*Voc), start=list(B=log(max(1e-8,min(outputDCcapacitance$capacitance))),C=log(1e-10),D=2), data=outputDCcapacitance)
	}, error=function(e) print("Failed restricted to positive gamma fit - second"))
	tryCatch({
		expfit <<- nlrob(capacitance ~ exp(B) + exp(C)*exp(D)*exp(exp(D)*Voc), start=list(B=coef(expfit)[[1]],C=coef(expfit)[[2]],D=coef(expfit)[[3]]), data=outputDCcapacitance)
	}, error=function(e) print("Failed restricted to positive gamma robust fit"))
}


getExpFit()

write.table(outputDCcapacitance, file="outputDCcapacitance.txt", append=TRUE, col.names=F, row.names=F, quote=F);

png(paste("DC-capacitance-", directory, ".png", sep=""), width=400, height=400)
par(mar=c(5,6,1,1))
plot(b$Voc, capacitance, ylab=bquote("Specific Capacitance (F/cm"^"2"*")"), xlab=bquote("V"["oc"]~"(V)"),cex.axis=1, cex.lab=1.4, log="y")
lines(outputDCcapacitance$Voc, predict(expfit),lwd=2, col="red")
graphics.off()

write.table(t(c("B","Ch0","gamma","GeomCh","ChemCh")), file="outputDC-fit.txt", append=FALSE, col.names=F, row.names=F);
eB=exp(coef(expfit)[[1]])
eCh0=exp(coef(expfit)[[2]])
egamma=exp(coef(expfit)[[3]])
Voc=max(outputDCcapacitance$Voc)
GeomCh=eB*Voc
ChemCh=eCh0*(exp(egamma*Voc)-1)
output <- t(c(eB, eCh0, egamma, GeomCh, ChemCh))
write.table(output, file="outputDC-fit.txt", append=TRUE, col.names=F, row.names=F)

c<- data.frame(b$Voc,capacitance)
d <- c[with(c, order(b.Voc)), ]
e <- d[1:(nrow(d)/2),]

f <- data.frame(d$b.Voc, d$capacitance)

names(f) <- c("Voc","capacitance")
g <- f
g$capacitance[g$capacitance < 0] <- 0

z <- approxfun(g$Voc, g$capacitance, method="linear", 0, 0)
integral=Vectorize(function(X)integrate(z,0,X)$value)

outputDCcharge <- data.frame(f$Voc, integral(f$Voc));
write.table(outputDCcharge, file="outputDCcharge.txt", append=TRUE, col.names=F, row.names=F, quote=F);

png(paste("DC-charge-", directory, ".png", sep=""), width=400, height=400)
par(mar=c(5,6,1,1))
plot(integral,range(f$Voc)[1], range(f$Voc)[2], ylab=bquote("Charge Density (C/cm"^"2"*")"), xlab=bquote("V"["oc"]~"(V)"),cex.axis=1, cex.lab=1.4)#, log="y")
graphics.off()

dataframe <- data.frame(Voc=g$Voc,capacitance=g$capacitance)

geometrical <- quantile(dataframe$capacitance, 0.05)
g$capacitance <- g$capacitance - geometrical

z <- approxfun(g$Voc, g$capacitance, method="linear", 0, 0)
integral=Vectorize(function(X)integrate(z,0,X)$value)

outputDCcharge <- data.frame(f$Voc, integral(f$Voc));
write.table(outputDCcharge, file="outputDCcharge-nogeom.txt", append=TRUE, col.names=F, row.names=F, quote=F);

png(paste("DC-nogeom-charge-", directory, ".png", sep=""), width=400, height=400)
par(mar=c(5,6,1,1))
plot(integral,range(f$Voc)[1], range(f$Voc)[2], ylab=bquote("Charge Density (C/cm"^"2"*")"), xlab=bquote("V"["oc"]~"(V)"),cex.axis=1, cex.lab=1.4)#, log="y")
graphics.off()
}
