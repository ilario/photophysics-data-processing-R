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

#title=gsub("_"," ",tail(unlist(strsplit(name,"-")),1))
filename=gsub(",","",gsub(":","",name))

library(RColorBrewer)
library(sfsmisc)
library(Hmisc)

ylim=lim.DCcapacitance.capacitance
xlim=lim.DCcapacitance.voltage
ylimnogeom=lim.DCcapacitance.nogeom.capacitance
xlimnogeom=lim.DCcapacitance.nogeom.voltage

output=list()
output.nogeom=list()

i <- 0
dirs <- list.dirs(recursive=FALSE)
dirs <- sub("./","",dirs)
legend=sub("-ig..-...-.","",sub("^0","",dirs))
#files.tpc <- list.files(path=".", pattern="*-outputDeltaV.txt$")
#files.tpv <- list.files(path=".", pattern="*-outputChargeDensityTPC.txt$")
colors=colorRampPalette(c("red","orange","springgreen","royalblue"))(max(length(dirs),3))
#brewer.pal(max(length(dirs),3),"Spectral")

jpeg(quality=98, paste(filename,"-DCs-capacitance.jpg",sep=""), width=image_width, height=image_height)
par(mar=c(5.1,7,2,2.1))
plot(NULL,xlim=xlim,ylim=ylim,cex.main=1.5,cex.axis=1.2,cex.lab=1.5,xlab="Voltage (V)",ylab="",# log="y", 
yaxt="n",xaxt="n")#main=paste(name,"DCs capacitance"), )
#eaxis(side=2,at=c(1e-12,1e-11,1e-10,1e-9,1e-8,1e-7,1e-6,1e-5,1e-4,1e-3,1e-2,0.1,1,10,100,1e3), cex.axis=1.2)
eaxis(side=2, cex.axis=1.2)
eaxis(side=1, cex.axis=1.2)
minor.tick(nx=10, ny=10)
title(ylab=bquote("Specific Capacitance (F/cm"^"2"*")"), mgp=c(5,1,0), cex.lab=1.5)

lapply(dirs, function(x) {print(x);
a <- read.table(paste(x,"/tpc/outputChargeDensityTPC.txt",sep=""),header=T)
charge <- mean(a$ChargeDensityTPC)
if(file.exists(file.path(x, "tpv", "outputDeltaVmixed.txt"))){
	        b <- read.table(file.path(x, "tpv", "outputDeltaVmixed.txt"), header=T)
}else{
	        b <- read.table(file.path(x, "tpv", "outputDeltaVloess.txt"), header=T)
}
output[[paste("Voc",sub("nm","",sub("-ig..-...-.","",sub("^0","",x))),sep="")]] <<- signif(b$Voc,5)
capacitance <- charge/b$deltaV
output[[sub("-ig..-...-.","",sub("^0","",x))]] <<- signif(capacitance,5)

points(b$Voc, capacitance, lwd=1, bg=colors[i+1], cex=2, pch=21+(i%%5))
 i <<- i+1
})
#abline(h=0)
legend(x="topleft",inset=0.05,legend,pch=seq(21,25), pt.bg=colors,pt.cex=2, cex=1.5, pt.lwd=2, lwd=4,col=colors, title=title, bg="gray90"# bty="n"
)
graphics.off()

maxlength = max(sapply(output,length))
output = lapply(output, function(x){length(x)=maxlength; print(x)})
output = as.data.frame(output,check.names=FALSE)
write.table(output, file=paste(filename,"-DCs-capacitance.csv",sep=""), row.names=FALSE, na="", sep=",")

i <- 0
jpeg(quality=98, paste(filename,"-DCs-nogeom-capacitance.jpg",sep=""), width=image_width, height=image_height)
par(mar=c(5.1,7,2,2.1))
plot(NULL,xlim=xlimnogeom,ylim=ylimnogeom,cex.main=1.5,cex.axis=1.2,cex.lab=1.5,xlab="Voltage (V)",ylab=bquote("Specific Capacitance (F/cm"^"2"*")"), #log="y", 
yaxt="n",xaxis="n")#main=paste(name,"DCs capacitance"), )
#eaxis(side=2,at=c(1e-12,1e-11,1e-10,1e-9,1e-8,1e-7,1e-6,1e-5,1e-4,1e-3,1e-2,0.1,1,10,100,1e3), cex.axis=1.2)
eaxis(side=2, cex.axis=1.2)
eaxis(side=1, cex.axis=1.2)
minor.tick(nx=10, ny=10)

lapply(dirs, function(x) {print(x);
a <- read.table(paste(x,"/tpc/outputChargeDensityTPC.txt",sep=""),header=T)
charge <- mean(a$ChargeDensityTPC)
if(file.exists(file.path(x, "tpv", "outputDeltaVmixed.txt"))){
	        b <- read.table(file.path(x, "tpv", "outputDeltaVmixed.txt"), header=T)
}else{
	        b <- read.table(file.path(x, "tpv", "outputDeltaVloess.txt"), header=T)
}
output.nogeom[[paste("Voc",sub("nm","",sub("-ig..-...-.","",sub("^0","",x))),sep="")]] <<- signif(b$Voc,5)
capacitance <- charge/b$deltaV
dataframe <- data.frame(Voc=b$Voc,capacitance=capacitance)

dataframe1 <- subset(dataframe, Voc < max(Voc)/2)
dataframe2 <- dataframe[1:2,]
dataframe <- unique(rbind(dataframe1, dataframe2))

geometrical <- mean(min(dataframe$capacitance),median(dataframe$capacitance))
capacitance <- capacitance - geometrical
#capacitance <- capacitance - min(capacitance)#mean(sort(capacitance)[1:3])

output.nogeom[[sub("-ig..-...-.","",sub("^0","",x))]] <<- signif(capacitance,5)
points(b$Voc, capacitance, lwd=1, bg=colors[i+1], cex=2, pch=21+(i%%5))
 i <<- i+1
})
#abline(h=0)
legend(x="topleft",inset=0.05,legend,pch=seq(21,25), pt.bg=colors,pt.cex=2, cex=1.5, pt.lwd=2, lwd=4,col=colors, title=#paste("DC capacitance\n","no geom. cap.\n",
       title, bg="gray90"#bty="n"
)
graphics.off()

output.nogeom = lapply(output.nogeom, function(x){length(x)=maxlength; print(x)})
output.nogeom = as.data.frame(output.nogeom,check.names=FALSE)
write.table(output.nogeom, file=paste(filename,"-DCs-nogeom-capacitance.csv",sep=""), row.names=FALSE, na="", sep=",")
