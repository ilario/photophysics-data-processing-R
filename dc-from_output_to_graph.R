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


a <- read.table("tpc/outputChargeDensityTPC.txt",header=T)
charge <- mean(a$ChargeDensityTPC)
b <- read.table("tpv/outputDeltaV.txt",header=T)
capacitance <- charge/b$deltaV
directory <- tail(strsplit(getwd(), "/")[[1]], n=1)

png(paste("DC-capacitance-", directory, ".png", sep=""), width=1000, heigh=600)
plot(b$Voc, capacitance, ylab="Specific Capacitance (F/cm2)", xlab="Voltage (V)", main=paste(directory,"DC capacitance"))
graphics.off()

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
#integrate(z, range(f$Voc)[1], range(f$Voc)[2])

png(paste("DC-charge-", directory, ".png", sep=""), width=1000, heigh=600)
plot(Vectorize(function(X)integrate(z,0,X)$value),range(f$Voc)[1], range(f$Voc)[2], ylab="Charge Density (C/cm2)", xlab="Voltage (V)", main=paste(directory,"DC charge"))
graphics.off()


#plot(f$Voc[-nrow(f)]?????,cumsum(f$capacitance[-nrow(f)]*diff(f$Voc)))


