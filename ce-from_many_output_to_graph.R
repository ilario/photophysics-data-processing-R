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
library(robustbase)
library(sfsmisc)
library(Hmisc)

ylim=c(1e-9,1.1e-7)
xlim=c(0,0.95)

i <- 0
dirs <- list.dirs(recursive=FALSE)
dirs <- sub("./","",dirs)
colors=colorRampPalette(c("red","orange","springgreen","royalblue"))(max(length(dirs),3))
data <- lapply(dirs, function(x) {print(x);
 a <- read.table(paste(x,"/ce/outputChargeDensityCE.txt",sep=""),header=T,stringsAsFactors=F)
 b<-strsplit(a$file, "_")
 c<-unlist(b)[length(b[[1]])*(1:length(a$file))]
 d<-as.numeric(gsub("mV", "", c))
 a$d <- d
 a <- a[with(a, order(a$d)),]
 exp <- nlrob(ChargeDensityCE~ A+C*exp(D*d), start=list(A=0,C=2e-9,D=9), data=a)
 g <- predict(exp,a$d)
 a$g <- g
 a})
names(data) <- dirs

png(paste(name,"-CEs-linlog.png",sep=""), width=640, height=640)
par(mar=c(5.1,5,2,2.1))
plot(NULL,xlim=xlim,ylim=ylim,cex.main=1.5,xlab="Voltage (V)",ylab=bquote("Extracted Charge Density (C/cm"^"2"*")"),  cex.lab=1.5, cex.axis=1.2, log="y", yaxt="n");#main=paste(name,"CEs"),
eaxis(side=2,at=c(1e-10,1e-9,1e-8,1e-7,1e-6,1e-5,1e-4,1e-3,1e-2,0.1,1,10,100,1e3), cex.axis=1.2)
minor.tick(nx=10)
lapply(dirs, function(x) {print(x);
 lines(data[[x]]$d, data[[x]]$g, col=colors[i+1],lwd=2)
 points(data[[x]]$d, data[[x]]$ChargeDensityCE, lwd=1, bg=colors[i+1], pch=21+i, cex=2)
# mtext(bquote(.(gsub("-outputChargeDensityCE.txt","",x))~": n" == .(signif(exp$coefficients["A"],3)) + 
#	      .(signif(exp$coefficients["C"],3)) ~ "e" ^ {.(signif(exp$coefficients["D"],3))~V}),side=3,line=-(i*2+4),cex=1.5,col=colors[i+1])
# mtext(bquote(.(x)~": n" == .(signif(exp$coefficients["A"],3)) + 
#	      .(signif(exp$coefficients["C"],3)) ~ "e" ^ {.(signif(exp$coefficients["D"],3))~V}),side=3,line=-(i*2+4),cex=1.5,col=colors[i+1])
 i <<- i+1
})
legend(x="bottomright",inset=0.05,sub("-ig..-...-.","",sub("^0","",dirs)), pch=seq(21,25), pt.bg=colors, col=colors, pt.cex=2, cex=1.5, pt.lwd=2, lwd=4, title=paste("Charge Extraction\n",sub("-"," - ",sub("_"," ",name))), bty="n")
graphics.off()

i<-0
png(paste(name,"-CEs.png",sep=""), width=640, height=640)
par(mar=c(5.1,7,2,2.1))
plot(NULL,xlim=xlim,ylim=ylim,cex.main=1.5,xlab="Voltage (V)", ylab="", cex.lab=1.5, cex.axis=1.2, yaxt="n");#main=paste(name,"CEs"),
eaxis(side=2, cex.axis=1.2)
minor.tick(nx=10, ny=10)
title(ylab=bquote("Extracted Charge Density (C/cm"^"2"*")"), mgp=c(5,1,0), cex.lab=1.5)
lapply(dirs, function(x) {print(x);
 lines(data[[x]]$d, data[[x]]$g, col=colors[i+1],lwd=2)
 points(data[[x]]$d, data[[x]]$ChargeDensityCE, lwd=1, bg=colors[i+1], pch=21+i, cex=2)
# mtext(bquote(.(gsub("-outputChargeDensityCE.txt","",x))~": n" == .(signif(exp$coefficients["A"],3)) + 
#	      .(signif(exp$coefficients["C"],3)) ~ "e" ^ {.(signif(exp$coefficients["D"],3))~V}),side=3,line=-(i*2+4),cex=1.5,col=colors[i+1])
# mtext(bquote(.(x)~": n" == .(signif(exp$coefficients["A"],3)) + 
#	      .(signif(exp$coefficients["C"],3)) ~ "e" ^ {.(signif(exp$coefficients["D"],3))~V}),side=3,line=-(i*2+4),cex=1.5,col=colors[i+1])
 i <<- i+1
})
legend(x="topleft",inset=0.05,sub("-ig..-...-.","",sub("^0","",dirs)), pch=seq(21,25), pt.bg=colors, col=colors, pt.cex=2, cex=1.5, pt.lwd=2, lwd=4, title=paste("Charge Extraction\n",sub("-"," - ",sub("_"," ",name))), bty="n")
graphics.off()

