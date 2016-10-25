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
#library(magicaxis)
library(sfsmisc)
library(Hmisc)


i <- 0
dirs <- list.dirs(recursive=FALSE)
dirs <- sub("./","",dirs)
colors=colorRampPalette(c("red","orange","springgreen","royalblue"))(max(length(dirs),3))
#brewer.pal(max(length(dirs),3),"Spectral")

data <- lapply(dirs, function(x) {print(x);
 a <- read.table(paste(x,"/tpc/outputChargeDensityTPC.txt",sep=""),header=T)
 charge <- mean(a$ChargeDensityTPC)
 b <- read.table(paste(x,"/tpv/outputDeltaV.txt",sep=""),header=T)
 capacitance <- charge/b$deltaV
 c<- data.frame(b$Voc,capacitance)
 d <- c[with(c, order(b.Voc)), ]
 e <- d[1:(nrow(d)/2),]
 f <- data.frame(d$b.Voc, d$capacitance)
 names(f) <- c("Voc","capacitance")
 f})
names(data) <- dirs

png(paste(name,"-DCs-charge-linlog.png",sep=""), width=640, height=640)
par(mar=c(5.1,5,2,2.1))
plot(NULL,xlim=c(0,1),ylim=c(3e-10,7e-8),cex.main=1.5,cex.lab=1.5, cex.axis=1.2, xlab="Voltage (V)",ylab=bquote("Charge Density (C/cm"^"2"*")"),  log="y", las=1, yaxt="n")#main=paste(name,"DCs charge"),
#magaxis(side=1:2, ratio=0.5, unlog=FALSE, labels=FALSE, tcl=-0.5)
eaxis(side=2,at=c(1e-10,1e-9,1e-8,1e-7,1e-6,1e-5,1e-4,1e-3,1e-2,0.1,1,10,100,1e3), cex.axis=1.2)
minor.tick(nx=10)#, tick.ratio=n)
lapply(dirs, function(x) {print(x);
g <- data[[x]]
g$capacitance[g$capacitance < 0] <- 0
z <- approxfun(g$Voc, g$capacitance, method="linear", 0, 0)
w <- Vectorize(function(X)integrate(z,0,X)$value)
curve(w,range(data[[x]]$Voc)[1], range(data[[x]]$Voc)[2], lwd=2, col=colors[i+1], add=T)
ww <- function(X)integrate(z,0,X)$value
x <- c(seq(range(data[[x]]$Voc)[1], range(data[[x]]$Voc)[2], 0.1), range(data[[x]]$Voc)[2])
www <- unlist(lapply(x, ww))
points(x, www, lwd=1, bg=colors[i+1], cex=2, pch=21+i)
i <<- i+1
})
legend(x="bottomright",inset=0.05,sub("-ig..-...-.","",dirs),pch=seq(21,25), lwd=4, pt.cex=2, pt.lwd=2, pt.bg=colors, cex=1.5, col=colors, title=paste("DC charge\n",sub("-"," - ",sub("_"," ",name))), bty="n")
graphics.off()

i<-0
png(paste(name,"-DCs-charge.png",sep=""), width=640, height=640)
par(mar=c(5.1,7,2,2.1))
plot(NULL,xlim=c(0,1),ylim=c(3e-10,7e-8),cex.main=1.5,cex.lab=1.5, cex.axis=1.2, xlab="Voltage (V)",ylab="", las=1, yaxt="n")
eaxis(side=2, cex.axis=1.2)
minor.tick(nx=10, ny=10)
title(ylab=bquote("Charge Density (C/cm"^"2"*")"), mgp=c(5,1,0), cex.lab=1.5)
lapply(dirs, function(x) {print(x);
g <- data[[x]]
g$capacitance[g$capacitance < 0] <- 0
z <- approxfun(g$Voc, g$capacitance, method="linear", 0, 0)
w <- Vectorize(function(X)integrate(z,0,X)$value)
curve(w,range(data[[x]]$Voc)[1], range(data[[x]]$Voc)[2], lwd=2, col=colors[i+1], add=T)
ww <- function(X)integrate(z,0,X)$value
x <- c(seq(range(data[[x]]$Voc)[1], range(data[[x]]$Voc)[2], 0.1), range(data[[x]]$Voc)[2])
www <- unlist(lapply(x, ww))
points(x, www, lwd=1, bg=colors[i+1], cex=2, pch=21+i)
i <<- i+1
})
legend(x="topleft",inset=0.05,sub("-ig..-...-.","",dirs),pch=seq(21,25), lwd=4, pt.cex=2, pt.lwd=2, pt.bg=colors, cex=1.5, col=colors, title=paste("DC charge\n",sub("-"," - ",sub("_"," ",name))), bty="n")
graphics.off()

