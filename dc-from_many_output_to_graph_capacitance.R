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
#files.tpc <- list.files(path=".", pattern="*-outputDeltaV.txt$")
#files.tpv <- list.files(path=".", pattern="*-outputChargeDensityTPC.txt$")
colors=brewer.pal(length(dirs),"Set1")

png(paste(name,"-DCs-capacitance.png",sep=""), width=800, height=640)

plot(NULL,xlim=c(0,1),ylim=c(0,7e-7),cex.main=1.5,xlab="Voltage (V)",ylab="Specific Capacitance (F/cm2)", main=paste(name,"DCs capacitance"))

lapply(dirs, function(x) {print(x);
a <- read.table(paste(x,"/outputChargeDensityTPC.txt",sep=""),header=T)
charge <- mean(a$ChargeDensityTPC)
b <- read.table(paste(x,"/outputDeltaV.txt",sep=""),header=T)
capacitance <- charge/b$deltaV

points(b$Voc, capacitance, lwd=2, pch=i, col=colors[i+1])
 i <<- i+1
})
legend(x="topleft",inset=0.05,dirs,pch=seq(0,10,1), col=colors)
graphics.off()




