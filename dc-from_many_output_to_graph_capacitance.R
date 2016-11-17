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
library(sfsmisc)
library(Hmisc)

ylim=c(2e-8,4e-7)
xlim=c(0,1)

i <- 0
dirs <- list.dirs(recursive=FALSE)
dirs <- sub("./","",dirs)
#files.tpc <- list.files(path=".", pattern="*-outputDeltaV.txt$")
#files.tpv <- list.files(path=".", pattern="*-outputChargeDensityTPC.txt$")
colors=colorRampPalette(c("red","orange","springgreen","royalblue"))(max(length(dirs),3))
#brewer.pal(max(length(dirs),3),"Spectral")

png(paste(name,"-DCs-capacitance.png",sep=""), width=640, height=640)
par(mar=c(5.1,5,2,2.1))
plot(NULL,xlim=xlim,ylim=ylim,cex.main=1.5,cex.axis=1.2,cex.lab=1.5,xlab="Voltage (V)",ylab=bquote("Specific Capacitance (F/cm"^"2"*")"), log="y", yaxt="n")#main=paste(name,"DCs capacitance"), )
eaxis(side=2,at=c(1e-10,1e-9,1e-8,1e-7,1e-6,1e-5,1e-4,1e-3,1e-2,0.1,1,10,100,1e3), cex.axis=1.2)
minor.tick(nx=10)

lapply(dirs, function(x) {print(x);
a <- read.table(paste(x,"/tpc/outputChargeDensityTPC.txt",sep=""),header=T)
charge <- mean(a$ChargeDensityTPC)
if(file.exists(file.path(x, "tpv", "outputDeltaVmixed.txt"))){
	        b <- read.table(file.path(x, "tpv", "outputDeltaVmixed.txt"), header=T)
}else{
	        b <- read.table(file.path(x, "tpv", "outputDeltaV.txt"), header=T)
}
capacitance <- charge/b$deltaV

points(b$Voc, capacitance, lwd=1, bg=colors[i+1], cex=2, pch=21+i)
 i <<- i+1
})
#abline(h=0)
legend(x="topleft",inset=0.1,sub("-ig..-...-.","",sub("^0","",dirs)),pch=seq(21,25), pt.bg=colors,pt.cex=2, cex=1.5, pt.lwd=2, lwd=4,col=colors, title=paste("DC capacitance\n","with geom. cap.\n",sub("-"," - ",sub("_"," ",name))), bty="n")
graphics.off()

i <- 0
png(paste(name,"-DCs-nogeom-capacitance.png",sep=""), width=640, height=640)
par(mar=c(5.1,5,2,2.1))
plot(NULL,xlim=xlim,ylim=ylim,cex.main=1.5,cex.axis=1.2,cex.lab=1.5,xlab="Voltage (V)",ylab=bquote("Specific Capacitance (F/cm"^"2"*")"), log="y", yaxt="n")#main=paste(name,"DCs capacitance"), )
eaxis(side=2,at=c(1e-10,1e-9,1e-8,1e-7,1e-6,1e-5,1e-4,1e-3,1e-2,0.1,1,10,100,1e3), cex.axis=1.2)
minor.tick(nx=10)

lapply(dirs, function(x) {print(x);
a <- read.table(paste(x,"/tpc/outputChargeDensityTPC.txt",sep=""),header=T)
charge <- mean(a$ChargeDensityTPC)
if(file.exists(file.path(x, "tpv", "outputDeltaVmixed.txt"))){
	        b <- read.table(file.path(x, "tpv", "outputDeltaVmixed.txt"), header=T)
}else{
	        b <- read.table(file.path(x, "tpv", "outputDeltaV.txt"), header=T)
}
capacitance <- charge/b$deltaV
capacitance <- capacitance - min(capacitance)
points(b$Voc, capacitance, lwd=1, bg=colors[i+1], cex=2, pch=21+i)
 i <<- i+1
})
#abline(h=0)
legend(x="topleft",inset=0.1,sub("-ig..-...-.","",sub("^0","",dirs)),pch=seq(21,25), pt.bg=colors,pt.cex=2, cex=1.5, pt.lwd=2, lwd=4,col=colors, title=paste("DC capacitance\n","no geom. cap.\n",sub("-"," - ",sub("_"," ",name))), bty="n")
graphics.off()
