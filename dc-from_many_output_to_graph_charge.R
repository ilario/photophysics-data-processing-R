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
#library(magicaxis)
library(sfsmisc)
library(Hmisc)

ylim=lim.DCcharge.charge
xlim=lim.DCcharge.voltage
ylimnogeom=lim.DCcharge.nogeom.charge
xlimnogeom=lim.DCcharge.nogeom.voltage

output=list()
output.nogeom=list()

i <- 0
dirs <- list.dirs(recursive=FALSE)
dirs <- sub("./","",dirs)
legend=sub("_.*","",sub("^0","",dirs))

# try to obtain the color from the file name
mycolors=gsub(".*-col_","",dirs[grepl("-col_", dirs)])
# if the color is not set, use the default one
if(!length(mycolors)){mycolors=brewer.pal(8,"Dark2")}

data <- lapply(dirs, function(x) {print(x);
 subdirs <- list.dirs(path=x, recursive=F)
 subdirs.tpc <- subdirs[grep("tpc", subdirs, ignore.case=T)]
 subdirs.tpv <- subdirs[grep("tpv", subdirs, ignore.case=T)]
 a <- read.table(paste(subdirs.tpc,"/outputChargeDensityTPC.txt",sep=""),header=T)
# in case TPC in dark and in sun are different, the choice of what to use is arbitrary, I would use the third quartile of all the TPC measurements
 charge <- quantile(a$ChargeDensityTPC, 0.75)
 b <- read.table(file.path(subdirs.tpv, "outputDeltaVprocessedForDC.txt"), header=T)

capacitance <- charge/b$deltaV
 c<- data.frame(b$Voc,capacitance)
 d <- c[with(c, order(b.Voc)), ]
 e <- d[1:(nrow(d)/2),]
 f <- data.frame(d$b.Voc, d$capacitance)
 names(f) <- c("Voc","capacitance")
 f})
names(data) <- dirs

jpeg(quality=98, paste(filename,"-DCs-charge-linlog.jpg",sep=""), width=image_width, height=image_height)
par(mar=c(5.1,7,2,2.1))
plot(NULL,xlim=xlim,ylim=ylim,cex.main=1.5,cex.lab=1.5, cex.axis=1.2, xlab="Voltage (V)",ylab=bquote("Charge Density (C/cm"^"2"*")"),  log="y", las=1, yaxt="n")#main=paste(name,"DCs charge"),
#magaxis(side=1:2, ratio=0.5, unlog=FALSE, labels=FALSE, tcl=-0.5)
eaxis(side=2,at=c(1e-12,1e-11,1e-10,1e-9,1e-8,1e-7,1e-6,1e-5,1e-4,1e-3,1e-2,0.1,1,10,100,1e3), cex.axis=1.2)
minor.tick(nx=10)#, tick.ratio=n)
lapply(dirs, function(x) {print(x);
g <- data[[x]]
g$capacitance[g$capacitance < 0] <- 0
z <- approxfun(g$Voc, g$capacitance, method="linear", 0, 0)
w <- Vectorize(function(X)integrate(z,0,X)$value)
curve(w,range(data[[x]]$Voc)[1], range(data[[x]]$Voc)[2], lwd=2, col=mycolors[i+1], add=T)
ww <- function(X)integrate(z,0,X)$value
xx <- c(seq(range(data[[x]]$Voc)[1], range(data[[x]]$Voc)[2], 0.1), range(data[[x]]$Voc)[2])
output[[paste("Voc",sub("nm","",sub("_.*","",sub("^0","",x))),sep="")]] <<- signif(xx,5)
www <- unlist(lapply(xx, ww))
output[[sub("_.*","",sub("^0","",x))]] <<- signif(www,5)
points(xx, www, lwd=1, bg=mycolors[i+1], cex=2, pch=21+(i%%5))
i <<- i+1
})
legend(x="bottomright",inset=0.05,legend,pch=seq(21,25), lwd=4, pt.cex=2, pt.lwd=2, pt.bg=mycolors, cex=1.5, col=mycolors, title=#paste("DC charge\n",
title,bg="gray90"#), bty="n"
)
graphics.off()

maxlength = max(sapply(output,length))
output = lapply(output, function(x){length(x)=maxlength; print(x)})
output = as.data.frame(output,check.names=FALSE)
write.table(output, file=paste(filename,"-DCs-charge.csv",sep=""), row.names=FALSE, na="", sep=",")


i<-0
jpeg(quality=98, paste(filename,"-DCs-charge.jpg",sep=""), width=image_width, height=image_height)
par(mar=c(5.1,7,2,2.1))
plot(NULL,xlim=xlim,ylim=ylim,cex.main=1.5,cex.lab=1.5, cex.axis=1.2, xlab="Voltage (V)",ylab="", las=1, yaxt="n", xaxt="n")
eaxis(side=2, cex.axis=1.2)
eaxis(side=1, cex.axis=1.2)
minor.tick(nx=10, ny=10)
title(ylab=bquote("Charge Density (C/cm"^"2"*")"), mgp=c(5,1,0), cex.lab=1.5)
lapply(dirs, function(x) {print(x);
g <- data[[x]]
g$capacitance[g$capacitance < 0] <- 0
z <- approxfun(g$Voc, g$capacitance, method="linear", 0, 0)
w <- Vectorize(function(X)integrate(z,0,X)$value)
curve(w,range(data[[x]]$Voc)[1], range(data[[x]]$Voc)[2], lwd=2, col=mycolors[i+1], add=T)
ww <- function(X)integrate(z,0,X)$value
xx <- c(seq(range(data[[x]]$Voc)[1], range(data[[x]]$Voc)[2], 0.1), range(data[[x]]$Voc)[2])
www <- unlist(lapply(xx, ww))
points(xx, www, lwd=1, bg=mycolors[i+1], cex=2, pch=21+(i%%5))
i <<- i+1
})
legend(x="topleft",inset=0.05,legend,pch=seq(21,25), lwd=4, pt.cex=2, pt.lwd=2, pt.bg=mycolors, cex=1.5, col=mycolors, title=#paste("DC charge\n","with geom. cap.\n",
title,bg="gray90"#), bty="n"
)
graphics.off()

i<-0
jpeg(quality=98, paste(filename,"-DCs-nogeom-charge-linlog.jpg",sep=""), width=image_width, height=image_height)
par(mar=c(5.1,7,2,2.1))
plot(NULL,xlim=xlimnogeom,ylim=ylimnogeom,cex.main=1.5,cex.lab=1.5, cex.axis=1.2, xlab="Voltage (V)",ylab=bquote("Charge Density (C/cm"^"2"*")"),  log="y", las=1, yaxt="n")#main=paste(name,"DCs charge"),
#magaxis(side=1:2, ratio=0.5, unlog=FALSE, labels=FALSE, tcl=-0.5)
eaxis(side=2,at=c(1e-12,1e-11,1e-10,1e-9,1e-8,1e-7,1e-6,1e-5,1e-4,1e-3,1e-2,0.1,1,10,100,1e3), cex.axis=1.2)
minor.tick(nx=10)#, tick.ratio=n)
lapply(dirs, function(x) {print(x);
g <- data[[x]]
g$capacitance[g$capacitance < 0] <- 0
dataframe <- data.frame(Voc=g$Voc,capacitance=g$capacitance)

geometrical <- quantile(dataframe$capacitance, 0.05)
g$capacitance <- g$capacitance - geometrical
#g$capacitance <- g$capacitance - min(g$capacitance)#mean(sort(g$capacitance)[1:3])
z <- approxfun(g$Voc, g$capacitance, method="linear", 0, 0)
w <- Vectorize(function(X)integrate(z,0,X)$value)
curve(w,range(data[[x]]$Voc)[1], range(data[[x]]$Voc)[2], lwd=2, col=mycolors[i+1], add=T)
ww <- function(X)integrate(z,0,X)$value
xx <- c(seq(range(data[[x]]$Voc)[1], range(data[[x]]$Voc)[2], 0.1), range(data[[x]]$Voc)[2])
output.nogeom[[paste("Voc",sub("nm","",sub("_.*","",sub("^0","",x))),sep="")]] <<- signif(xx,5)
www <- unlist(lapply(xx, ww))
output.nogeom[[sub("_.*","",sub("^0","",x))]] <<- signif(www,5)
points(xx, www, lwd=1, bg=mycolors[i+1], cex=2, pch=21+(i%%5))
i <<- i+1
})
legend(x="bottomright",inset=0.05,legend,pch=seq(21,25), lwd=4, pt.cex=2, pt.lwd=2, pt.bg=mycolors, cex=1.5, col=mycolors, title=#paste("DC charge\n","no geom. cap.\n",
       title,bg="gray90"#, bty="n"
)
graphics.off()

output.nogeom = lapply(output.nogeom, function(x){length(x)=maxlength; print(x)})
output.nogeom = as.data.frame(output.nogeom,check.names=FALSE)
write.table(output.nogeom, file=paste(filename,"-DCs-nogeom-charge.csv",sep=""), row.names=FALSE, na="", sep=",")


i<-0
jpeg(quality=98, paste(filename,"-DCs-nogeom-charge.jpg",sep=""), width=image_width, height=image_height)
par(mar=c(5.1,7,2,2.1))
plot(NULL,xlim=xlimnogeom,ylim=ylimnogeom,cex.main=1.5,cex.lab=1.5, cex.axis=1.2, xlab="Voltage (V)",ylab="", las=1, yaxt="n", xaxt="n")
eaxis(side=2, cex.axis=1.2)
eaxis(side=1, cex.axis=1.2)
minor.tick(nx=10, ny=10)
title(ylab=bquote("Charge Density (C/cm"^"2"*")"), mgp=c(5,1,0), cex.lab=1.5)
lapply(dirs, function(x) {print(x);
g <- data[[x]]
g$capacitance[g$capacitance < 0] <- 0
dataframe <- data.frame(Voc=g$Voc,capacitance=g$capacitance)

geometrical <- quantile(dataframe$capacitance,0.05)
g$capacitance <- g$capacitance - geometrical
#g$capacitance <- g$capacitance - min(g$capacitance) #mean(sort(g$capacitance)[1:3])
z <- approxfun(g$Voc, g$capacitance, method="linear", 0, 0)
w <- Vectorize(function(X)integrate(z,0,X)$value)
curve(w,range(data[[x]]$Voc)[1], range(data[[x]]$Voc)[2], lwd=2, col=mycolors[i+1], add=T)
ww <- function(X)integrate(z,0,X)$value
xx <- c(seq(range(data[[x]]$Voc)[1], range(data[[x]]$Voc)[2], 0.1), range(data[[x]]$Voc)[2])
www <- unlist(lapply(xx, ww))
points(xx, www, lwd=1, bg=mycolors[i+1], cex=2, pch=21+(i%%5))
i <<- i+1
})
legend(x="topleft",inset=0.05,legend,pch=seq(21,25), lwd=4, pt.cex=2, pt.lwd=2, pt.bg=mycolors, cex=1.5, col=mycolors, title=#paste("DC charge\n","no geom. cap.\n",
title,bg="gray90"#, bty="n"
)
graphics.off()
