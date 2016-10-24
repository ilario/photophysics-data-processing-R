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

ylim=c(3e-7,1e-3)
xlim=c(0,2.5e-7)

library(robustbase)
#library(RColorBrewer)
library(minpack.lm)

i <- 0
dirs <- list.dirs(recursive=FALSE)
colors=colorRampPalette(c("red","orange","springgreen","royalblue"))(max(length(dirs),3))
#brewer.pal(max(length(dirs),3),"Spectral")
dirs <- sub("./","",dirs)

#lapply(dirs, function(x) {print(x);
# a <- read.table(paste(x,"/outputDCcharge.txt",sep=""),header=T,stringsAsFactors=F)
#  lo <- loess(a$chargeDC~a$Voc, span=0.9)
#       	exp <- nlrob(chargeDC~ A+C*exp(D*Voc), start=list(A=0,C=2e-9,D=9), data=a)
# expend <- nlsLM(chargeDC~ A+C*exp(D*Voc), start=list(A=coef(exp)["A"],C=coef(exp)["C"],D=coef(exp)["D"]), data=a[round(length(a$Voc)/2):length(a$Voc),])
# png(paste(x,"-CEs.png",sep=""), width=640, height=640)
# plot(NULL,xlim=c(0,1),ylim=c(0,2e-7),cex.main=1.5,xlab="Voltage (V)",ylab="Charge Density (C/cm2)", main=paste(x,"DC fitted"));
# points(a, lwd=2, pch=1, col=colors[i+1])
# lines(a$Voc,predict(exp))
# lines(a$Voc,predict(lo),col="green")
# lines(a$Voc[round(length(a$Voc)/2):length(a$Voc)],predict(expend),col="red")
# graphics.off()
# i <<- i+1
#})

i <- 0
png(paste(name,"-TPVDCs.png",sep=""), width=640, height=640)
par(mar=c(5.1,5,4.1,2.1))
plot(1,xlim=xlim,ylim=ylim,cex.main=1.5,xlab=bquote("Charge Density (C/cm"^"2"*")"), ylab="Life-time (s)", main=paste(name,"TPV decay vs Charge from DC"),cex.lab=1.5,cex.axis=1.2,log="y");
lapply(dirs, function(x) {print(x);
 a <- read.table(paste(x,"/outputDCcharge.txt",sep=""),header=T,stringsAsFactors=F)
 lo <- loess(a$chargeDC~a$Voc,span=0.9)
 exp <- nlrob(chargeDC~ A+C*exp(D*Voc), start=list(A=0,C=2e-9,D=9), data=a)
 expend <- nlsLM(chargeDC~ A+C*exp(D*Voc), start=list(A=coef(exp)["A"],C=coef(exp)["C"],D=coef(exp)["D"]), data=a[round(length(a$Voc)/2):length(a$Voc),])

fulloutput <- read.table(paste(x,"/tpv/output-monoexp.txt",sep=""), header=TRUE);
n<-tail(grep("file",fulloutput[,1]),n=1)
tpv <- read.table(paste(x,"/tpv/output-monoexp.txt",sep=""), header=TRUE, skip=ifelse(length(n),n,0));
#importante che la variabile in new abbia lo stesso nome di quella fittata
new <- data.frame(Voc = tpv$Voc)
charge <- (predict(lo,tpv$Voc)+predict(exp,new))/2
new2 <- data.frame(Voc = tpv$Voc[is.na(charge)])
charge[is.na(charge)] <- (predict(exp,new2) + predict(expend,new2))/2
lo<-loess(tpv$T~charge,span=0.5)
lines(charge, predict(lo), lwd=2, col=colors[i+1])
points(charge, tpv$T, lwd=2, col=colors[i+1]);
 i <<- i+1
})
legend(x="topright",inset=0.05,sub("-ig..-...-.","",dirs),pch=1, lwd=4, pt.lwd=2, pt.cex=2, col=colors,cex=1.5)
graphics.off()