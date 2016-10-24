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

#name=""


library(RColorBrewer)
i <- 0
dirs <- list.dirs(recursive=FALSE)
dirs <- sub("./","",dirs)
colors=brewer.pal(max(length(dirs),3),"Spectral")

png(paste(name,"-DCs-charge.png",sep=""), width=640, height=640)
par(mar=c(5.1,5,4.1,2.1))
plot(NULL,xlim=c(0,1),ylim=c(1e-10,1.2e-7),cex.main=1.5,cex.lab=1.5, cex.axis=1.2, xlab="Voltage (V)",ylab=bquote("Charge Density (C/cm"^"2"*")"), main=paste(name,"DCs charge"), log="y")

lapply(dirs, function(x) {print(x);
a <- read.table(paste(x,"/tpc/outputChargeDensityTPC.txt",sep=""),header=T)
charge <- mean(a$ChargeDensityTPC)
b <- read.table(paste(x,"/tpv/outputDeltaV.txt",sep=""),header=T)
capacitance <- charge/b$deltaV

#points(b$Voc, capacitance, lwd=2, pch=i, col=colors[i+1])

c<- data.frame(b$Voc,capacitance)
d <- c[with(c, order(b.Voc)), ]
e <- d[1:(nrow(d)/2),]

f <- data.frame(d$b.Voc, d$capacitance)

names(f) <- c("Voc","capacitance")
g <- f
g$capacitance[g$capacitance < 0] <- 0

z <- approxfun(g$Voc, g$capacitance, method="linear", 0, 0)

w <- Vectorize(function(X)integrate(z,0,X)$value)

curve(w,range(f$Voc)[1], range(f$Voc)[2], lwd=2, col=colors[i+1], add=T)

ww <- function(X)integrate(z,0,X)$value
x <- c(seq(range(f$Voc)[1], range(f$Voc)[2], 0.1), range(f$Voc)[2])
www <- unlist(lapply(x, ww))
points(x, www, lwd=2, col=colors[i+1], cex=2)

i <<- i+1
})
legend(x="topleft",inset=0.05,sub("-ig..-...-.","",dirs),pch=1, lwd=4, pt.cex=2, pt.lwd=2, cex=1.5, col=colors)
graphics.off()
