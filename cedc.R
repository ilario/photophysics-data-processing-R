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

cea <- read.table("ce/outputChargeDensityCE.txt",header=T,stringsAsFactors=F)
ceb<-strsplit(cea$file, "_")
cec<-unlist(ceb)[length(ceb[[1]])*(1:length(cea$file))]
ced<-as.numeric(gsub("mV", "", cec))
directory <- tail(strsplit(getwd(), "/")[[1]], n=1)

a <- read.table("tpc/outputChargeDensityTPC.txt",header=T)
charge <- mean(a$ChargeDensityTPC)
b <- read.table("tpv/outputDeltaV.txt",header=T)
capacitance <- charge/b$deltaV

c<- data.frame(b$Voc,capacitance)
d <- c[with(c, order(b.Voc)), ]
e <- d[1:(nrow(d)/2),]

#fit <- nls(capacitance ~ cbind(1, exp(C*(b.Voc))), trace=F, alg="plinear", start=list(C=1), data=e)
#A = coef(fit)[".lin1"]
#B = coef(fit)[".lin2"]
#C = coef(fit)["C"]

#f <- data.frame(d$b.Voc, d$capacitance - A)
f <- data.frame(d$b.Voc, d$capacitance)
names(f) <- c("Voc","capacitance")
#g <- subset(f,capacitance>0)
g <- f
g$capacitance[g$capacitance < 0] <- 0

z <- approxfun(g$Voc, g$capacitance, method="linear", 0, 0)


ChargeDensityCE <- cea$ChargeDensityCE*1e9
png(paste("DC-CE-", directory, ".png", sep=""), heigh=600, width=600)
par(mar=c(5.1,5,4.1,2.1))
plot(Vectorize(function(X)integrate(z,0,X)$value*1e9),xlim=c(min(ced,range(f$Voc)[1]),max(ced,range(f$Voc)[2]))#, ylim=c(min(0,ChargeDensityCE),max(ChargeDensityCE))
     , ylab=bquote("Charge Density (nC/cm"^"2"*")"), xlab="Voltage (V)", main=paste(directory,"DC and CE"), cex.main=1.5, cex.lab=1.5, cex.axis=1.5)
points(ced, ChargeDensityCE)
legend(x="topleft",inset=0.1,c("Differential Charging", "Charge Extraction"), lty=c(1,NA), pch=c(NA,1), lwd=2, cex=1.5)
graphics.off()
