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
a <- read.table(file.path(tpcdir, "outputChargeDensityTPC.txt"), header=T)
charge <- mean(a$ChargeDensityTPC)
b <- read.table(file.path(tpvdir, "outputDeltaV.txt"), header=T)
capacitance <- charge/b$deltaV
directory <- tail(strsplit(getwd(), "/")[[1]], n=1)

png(paste("DC-capacitance-", directory, ".png", sep=""), width=400, heigh=400)
par(mar=c(5,6,1,1))
plot(b$Voc, capacitance, ylab=bquote("Specific Capacitance (F/cm"^"2"*")"), xlab=bquote("V"["oc"]~"(V)"),cex.axis=1.4, cex.lab=1.4, ylim=c(0,max(capacitance)))
graphics.off()

c<- data.frame(b$Voc,capacitance)
d <- c[with(c, order(b.Voc)), ]
e <- d[1:(nrow(d)/2),]

f <- data.frame(d$b.Voc, d$capacitance)

names(f) <- c("Voc","capacitance")
g <- f
g$capacitance[g$capacitance < 0] <- 0

z <- approxfun(g$Voc, g$capacitance, method="linear", 0, 0)

png(paste("DC-charge-", directory, ".png", sep=""), width=400, heigh=400)
par(mar=c(5,6,1,1))
plot(Vectorize(function(X)integrate(z,0,X)$value),range(f$Voc)[1], range(f$Voc)[2], ylab=bquote("Charge Density (C/cm"^"2"*")"), xlab=bquote("V"["oc"]~"(V)"),cex.axis=1.4, cex.lab=1.4)
graphics.off()
}
