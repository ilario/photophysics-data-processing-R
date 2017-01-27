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

ylim=limlifetime
xlim=limtpvdccharge

output=list()
output.nogeom=list()

library(robustbase)
#library(RColorBrewer)
library(minpack.lm)
library(sfsmisc)
library(Hmisc)

#i <- 0
dirs <- list.dirs(recursive=FALSE)
colors=colorRampPalette(c("red","orange","springgreen","royalblue"))(max(length(dirs),3))
##brewer.pal(max(length(dirs),3),"Spectral")
dirs <- sub("./","",dirs)
legend=sub("-ig..-...-.","",sub("^0","",dirs))
#lapply(dirs, function(x) {print(x);
# a <- read.table(paste(x,"/outputDCcharge.txt",sep=""),header=T,stringsAsFactors=F)
#  lo <- loess(a$ChargeDensityDC~a$Voc, span=0.9)
#       	exp <- nlrob(ChargeDensityDC~ B*Voc + C*(exp(D*Voc)-1), start=list(B=1e-8,C=2e-9,D=9), data=a)
# expend <- nlsLM(ChargeDensityDC~ C*(exp(D*Voc)-A), start=list(A=1,C=coef(exp)["C"],D=coef(exp)["D"]), data=a[round(length(a$Voc)/2):length(a$Voc),])
# jpeg(quality=98, paste(x,"-DCfitting.jpg",sep=""), width=640, height=480)
# plot(NULL,xlim=c(0,1),ylim=c(0,2e-7),cex.main=1.5,xlab="Voltage (V)",ylab="Charge Density (C/cm2)", main=paste(x,"DC fitted"));
# points(a, lwd=2, pch=1, col=colors[i+1])
# lines(a$Voc,predict(exp))
# lines(a$Voc,predict(lo),col="green")
# lines(a$Voc[round(length(a$Voc)/2):length(a$Voc)],predict(expend),col="red")
# graphics.off()
# i <<- i+1
#})

i <- 0
jpeg(quality=98, paste(filename,"-TPVDCs.jpg",sep=""), width=640, height=480)
par(mar=c(5.1,7,2,2.1))
plot(1,xlim=xlim,ylim=ylim,cex.main=1.5,xlab=bquote("Charge Density (C/cm"^"2"*")"), ylab="Life-time (s)",cex.lab=1.5,cex.axis=1.2,log="y", yaxt="n", xaxt="n");#, main=paste(name,"TPV decay vs Charge from DC")
eaxis(side=2,at=c(1e-10,1e-9,1e-8,1e-7,1e-6,1e-5,1e-4,1e-3,1e-2,0.1,1,10,100,1e3), cex.axis=1.2)
eaxis(side=1, cex.axis=1.2)
minor.tick(nx=10)

lapply(dirs, function(x) {print(x);
 a <- read.table(paste(x,"/outputDCcharge.txt",sep=""),header=T,stringsAsFactors=F)
 lo <- loess(a$ChargeDensityDC~a$Voc,span=0.9)
 expend <- nlsLM(ChargeDensityDC~ A+C*exp(D*Voc), start=list(A=-1e-10,C=1e-10,D=8), data=a[round(length(a$Voc)/2):length(a$Voc),])
tryCatch({
 exp <- nlrob(ChargeDensityDC~ C*(exp(D*Voc)-1), start=list(C=coef(expend)["C"],D=coef(expend)["D"]), data=a)
}, error=function(e) {print("FAILED ZEROTH FIT")});
tryCatch({
 exp <- nlrob(ChargeDensityDC~ A+C*exp(D*Voc), start=list(A=coef(expend)["A"],C=coef(expend)["C"],D=coef(expend)["D"]), data=a)
}, error=function(e) {print("FAILED FIRST FIT")});
tryCatch({
 exp <- nlrob(ChargeDensityDC~ B*Voc+C*(exp(D*Voc)-1), start=list(B=1e-9,C=coef(expend)["C"],D=coef(expend)["D"]), data=a)
}, error=function(e) {print("FAILED SECOND FIT")});
tryCatch({
 exp <- nlrob(ChargeDensityDC~ B*Voc+C*(exp(D*Voc)-1), start=list(B=1e-9,C=1e-10,D=8), data=a)
}, error=function(e) {print("FAILED THIRD FIT")});
tryCatch({
 exp <- nlrob(ChargeDensityDC~ B*Voc+C*(exp(D*Voc)-1), start=list(B=1e-8,C=1e-9,D=1), data=a)
}, error=function(e) {print("FAILED THIRD FIT")});
tryCatch({
 exp <- nlrob(ChargeDensityDC~ A+B*Voc+C*exp(D*Voc), start=list(A=0,B=1e-9,C=coef(expend)["C"],D=coef(expend)["D"]), data=a)
}, error=function(e) {print("FAILED FOURTH FIT")});

filex <- file.path(x, "tpv", "output-monoexp.txt")

fulloutput <- read.table(filex, header=TRUE);
n<-tail(grep("file",fulloutput[,1]),n=1)
tpv <- read.table(filex, header=TRUE, skip=ifelse(length(n),n,0));
#importante che la variabile in new abbia lo stesso nome di quella fittata
new <- data.frame(Voc = tpv$Voc)
charge <- (predict(lo,tpv$Voc)+predict(exp,new))/2
new2 <- data.frame(Voc = tpv$Voc[is.na(charge)])
charge[is.na(charge)] <- (predict(exp,new2) + predict(expend,new2))/2
output[[paste("Charge",sub("nm","",sub("-ig..-...-.","",sub("^0","",x))),sep="")]] <<- signif(charge,5)
output[[sub("-ig..-...-.","",sub("^0","",x))]] <<- signif(tpv$T,5)
lo<-loess(tpv$T~charge,span=0.3)
lines(charge, predict(lo), lwd=2, col=colors[i+1])
points(charge, tpv$T, lwd=1, bg=colors[i+1], cex=2, pch=21+i);
 i <<- i+1
})
legend(x="topright",inset=0.05,legend,pch=seq(21,25), pt.bg=colors, lwd=4, pt.lwd=2, pt.cex=2, col=colors,cex=1.5, title=#paste("TPV vs DC\n","with geom. cap.\n",
title,bg="gray90"#), bty="n"
)
graphics.off()

maxlength = max(sapply(output,length))
output = lapply(output, function(x){length(x)=maxlength; print(x)})
output = as.data.frame(output,check.names=FALSE)
write.table(output, file=paste(filename,"-TPVDCs.csv",sep=""), row.names=FALSE, na="", sep=",")


i <- 0
jpeg(quality=98, paste(filename,"-TPVDCs-nogeom.jpg",sep=""), width=640, height=480)
par(mar=c(5.1,7,2,2.1))
plot(1,xlim=xlim,ylim=ylim,cex.main=1.5,xlab=bquote("Charge Density (C/cm"^"2"*")"), ylab="Life-time (s)",cex.lab=1.5,cex.axis=1.2,log="y", yaxt="n", xaxt="n");#, main=paste(name,"TPV decay vs Charge from DC")
eaxis(side=2,at=c(1e-10,1e-9,1e-8,1e-7,1e-6,1e-5,1e-4,1e-3,1e-2,0.1,1,10,100,1e3), cex.axis=1.2)
eaxis(side=1, cex.axis=1.2)
minor.tick(nx=10)

lapply(dirs, function(x) {print(x);
 a <- read.table(paste(x,"/outputDCcharge-nogeom.txt",sep=""),header=T,stringsAsFactors=F)
 lo <- loess(a$ChargeDensityDC~a$Voc,span=0.9)
 expend <- nlsLM(ChargeDensityDC~ A+C*exp(D*Voc), start=list(A=-1e-10,C=1e-10,D=8), data=a[round(length(a$Voc)/2):length(a$Voc),])
tryCatch({
 exp <- nlrob(ChargeDensityDC~ C*(exp(D*Voc)-1), start=list(C=coef(expend)["C"],D=coef(expend)["D"]), data=a)
}, error=function(e) {print("FAILED ZEROTH nogeom FIT")});
tryCatch({
 exp <- nlrob(ChargeDensityDC~ A+C*exp(D*Voc), start=list(A=coef(expend)["A"],C=coef(expend)["C"],D=coef(expend)["D"]), data=a)
}, error=function(e) {print("FAILED FIRST nogeom FIT")});
tryCatch({
 exp <- nlrob(ChargeDensityDC~ A+C*exp(D*Voc), start=list(A=0,C=1e-10,D=8), data=a)
}, error=function(e) {print("FAILED SECOND nogeom FIT")});

filex <- file.path(x, "tpv", "output-monoexp.txt")
fulloutput <- read.table(filex, header=TRUE);
n<-tail(grep("file",fulloutput[,1]),n=1)
tpv <- read.table(filex, header=TRUE, skip=ifelse(length(n),n,0));
#importante che la variabile in new abbia lo stesso nome di quella fittata
new <- data.frame(Voc = tpv$Voc)
charge <- (predict(lo,tpv$Voc)+predict(exp,new))/2
new2 <- data.frame(Voc = tpv$Voc[is.na(charge)])
charge[is.na(charge)] <- (predict(exp,new2) + predict(expend,new2))/2
output.nogeom[[paste("Charge",sub("nm","",sub("-ig..-...-.","",sub("^0","",x))),sep="")]] <<- signif(charge,5)
output.nogeom[[sub("-ig..-...-.","",sub("^0","",x))]] <<- signif(tpv$T,5)
lo<-loess(tpv$T~charge,span=0.3)
lines(charge, predict(lo), lwd=2, col=colors[i+1])
points(charge, tpv$T, lwd=1, bg=colors[i+1], cex=2, pch=21+i);
i <<- i+1
})
legend(x="topright",inset=0.05,legend,pch=seq(21,25), pt.bg=colors, lwd=4, pt.lwd=2, pt.cex=2, col=colors,cex=1.5, title=#paste("TPV vs DC\n","no geom. cap.\n",
title, bg="gray90"#, bty="n"
)
graphics.off()

output.nogeom = lapply(output.nogeom, function(x){length(x)=maxlength; print(x)})
output.nogeom = as.data.frame(output.nogeom,check.names=FALSE)
write.table(output.nogeom, file=paste(filename,"-TPVDCs-nogeom.csv",sep=""), row.names=FALSE, na="", sep=",")
