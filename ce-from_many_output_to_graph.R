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

#title=gsub("-","\n\n",gsub("_"," ",name))
title=gsub("_"," ",tail(unlist(strsplit(name,"-")),1))
filename=gsub(",","",gsub(":","",name))

library(RColorBrewer)
library(robustbase)
library(sfsmisc)
library(Hmisc)

ylim=limcecharge
xlim=limvoltage

output=list()

i <- 0
dirs <- list.dirs(recursive=FALSE)
dirs <- sub("./","",dirs)
legend=sub("-ig..-...-.","",sub("^0","",dirs))
colors=colorRampPalette(c("red","orange","springgreen","royalblue"))(max(length(dirs),3))
data <- lapply(dirs, function(x) {print(x);
 a <- read.table(paste(x,"/ce/outputChargeDensityCE.txt",sep=""),header=T,stringsAsFactors=F)
 output[[paste("Voc",sub("nm","",sub("-ig..-...-.","",sub("^0","",x))),sep="")]] <<- a$Voc
 a <- a[with(a, order(a$Voc)),]
 exp <- nlrob(ChargeDensityCE~ A+C*exp(D*Voc), start=list(A=0,C=1e-10,D=9), data=a)
 g <- predict(exp,a$Voc)
 a$g <- g
 output[[sub("-ig..-...-.","",sub("^0","",x))]] <<- signif(g,5)
 a})
names(data) <- dirs

maxlength = max(sapply(output,length))
output = lapply(output, function(x){length(x)=maxlength; print(x)})
output = as.data.frame(output,check.names=FALSE)
write.table(output, file=paste(filename,"-CEs.csv",sep=""), row.names=FALSE, na="", sep=",")

jpeg(quality=98, paste(filename,"-CEs-linlog.jpg",sep=""), width=640, height=480)
par(mar=c(5.1,5,2,2.1))
plot(NULL,xlim=xlim,ylim=ylim,cex.main=1.5,xlab="Voltage (V)",ylab=bquote("Extracted Charge Density (C/cm"^"2"*")"),  cex.lab=1.5, cex.axis=1.2, log="y", yaxt="n");#main=paste(name,"CEs"),
eaxis(side=2,at=c(1e-10,1e-9,1e-8,1e-7,1e-6,1e-5,1e-4,1e-3,1e-2,0.1,1,10,100,1e3), cex.axis=1.2)
minor.tick(nx=10)
lapply(dirs, function(x) {print(x);
 lines(data[[x]]$Voc, data[[x]]$g, col=colors[i+1],lwd=2)
 points(data[[x]]$Voc, data[[x]]$ChargeDensityCE, lwd=1, bg=colors[i+1], pch=21+i, cex=2)
# mtext(bquote(.(gsub("-outputChargeDensityCE.txt","",x))~": n" == .(signif(exp$coefficients["A"],3)) + 
#	      .(signif(exp$coefficients["C"],3)) ~ "e" ^ {.(signif(exp$coefficients["D"],3))~V}),side=3,line=-(i*2+4),cex=1.5,col=colors[i+1])
# mtext(bquote(.(x)~": n" == .(signif(exp$coefficients["A"],3)) + 
#	      .(signif(exp$coefficients["C"],3)) ~ "e" ^ {.(signif(exp$coefficients["D"],3))~V}),side=3,line=-(i*2+4),cex=1.5,col=colors[i+1])
 i <<- i+1
})
legend(x="bottomright",inset=0.05,legend, pch=seq(21,25), pt.bg=colors, col=colors, pt.cex=2, cex=1.5, pt.lwd=2, lwd=4, title=title, bg="gray90"#, bty="n"
)
graphics.off()

i<-0
jpeg(quality=98, paste(filename,"-CEs.jpg",sep=""), width=640, height=480)
par(mar=c(5.1,7,2,2.1))
plot(NULL,xlim=xlim,ylim=ylim,cex.main=1.5,xlab="Voltage (V)", ylab="", cex.lab=1.5, cex.axis=1.2, yaxt="n");#main=paste(name,"CEs"),
eaxis(side=2, cex.axis=1.2)
minor.tick(nx=10, ny=10)
title(ylab=bquote("Extracted Charge Density (C/cm"^"2"*")"), mgp=c(5,1,0), cex.lab=1.5)
lapply(dirs, function(x) {print(x);
 lines(data[[x]]$Voc, data[[x]]$g, col=colors[i+1],lwd=2)
 points(data[[x]]$Voc, data[[x]]$ChargeDensityCE, lwd=1, bg=colors[i+1], pch=21+i, cex=2)
# mtext(bquote(.(gsub("-outputChargeDensityCE.txt","",x))~": n" == .(signif(exp$coefficients["A"],3)) + 
#	      .(signif(exp$coefficients["C"],3)) ~ "e" ^ {.(signif(exp$coefficients["D"],3))~V}),side=3,line=-(i*2+4),cex=1.5,col=colors[i+1])
# mtext(bquote(.(x)~": n" == .(signif(exp$coefficients["A"],3)) + 
#	      .(signif(exp$coefficients["C"],3)) ~ "e" ^ {.(signif(exp$coefficients["D"],3))~V}),side=3,line=-(i*2+4),cex=1.5,col=colors[i+1])
 i <<- i+1
})
legend(x="topleft",inset=0.05,legend, pch=seq(21,25), pt.bg=colors, col=colors, pt.cex=2, cex=1.5, pt.lwd=2, lwd=4, title=title, bg="gray90"#bty="n"
)
graphics.off()

