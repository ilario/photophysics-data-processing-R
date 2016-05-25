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

ylim=c(3e-7,6e-4)
xlim=c(0,2.5e-7)

library(robustbase)
library(RColorBrewer)
library(minpack.lm)

i <- 0
dirs <- list.dirs(recursive=FALSE)
colors=brewer.pal(max(length(dirs),3),"Set1")
dirs <- sub("./","",dirs)

lapply(dirs, function(x) {print(x);
 a <- read.table(paste(x,"-outputChargeDensityCE.txt",sep=""),header=T,stringsAsFactors=F)
 b<-strsplit(a$file, "_")
 c<-unlist(b)[length(b[[1]])*(1:length(a$file))]
 d<-as.numeric(gsub("mV", "", c))
 lo <- loess(a$ChargeDensityCE~d,span=0.9)
 a$d <- d
#tryCatch({
       	exp <- nlrob(ChargeDensityCE~ A+C*exp(D*d), start=list(A=0,C=2e-9,D=9), data=a)
 expend <- nlsLM(ChargeDensityCE~ A+C*exp(D*d), start=list(A=coef(exp)["A"],C=coef(exp)["C"],D=coef(exp)["D"]), data=a[round(length(a$file)/2):length(a$file),])
#}, error=function(e) {
#	library(minipack.lm)
#	exp <- nlsLM(ChargeDensityCE~ A+C*exp(D*d), start=list(A=0,C=2e-9,D=5), data=a)
 #sinexp <- nls(ChargeDensityCE~ cbind(1,d,sin(D*d),exp(F*d)), start=c(D=-7,F=6), data=a, alg="plinear")
# sinexp2 <- nlrob(ChargeDensityCE~ A+B*d+C*sin(D*d)+E*exp(F*d), start=list(A=coef(sinexp)[".lin1"],B=coef(sinexp)[".lin.d"],C=coef(sinexp)[".lin3"],D=coef(sinexp)["D"],E=coef(sinexp)[".lin4"],F=coef(sinexp)["F"]), data=a)
# exp <- nls(ChargeDensityCE~ C*exp(F*d), start=c(C=1,F=6), data=a)
# powexp <- nls(ChargeDensityCE~ cbind(poly(d,2),exp(F*d)), start=c(F=coef(exp)["F"]), data=a, alg="plinear")
# sinexp2 <- nlrob(ChargeDensityCE~ C*sin(D*d)+E*exp(F*d), start=list(C=coef(sinexp)[".lin1"],D=coef(sinexp)["D"],E=coef(sinexp)[".lin2"],F=coef(sinexp)["F"]), data=a)
 png(paste(x,"-CEs.png",sep=""), width=800, height=640)
 plot(NULL,xlim=c(0,1),ylim=c(0,2e-7),cex.main=1.5,xlab="Voltage (V)",ylab="Extracted Charge Density (C/cm2)", main=paste(x,"CEs"));
 points(d, a$ChargeDensityCE, lwd=2, pch=i, col=colors[i+1])
# lines(d,predict(lo))
# lines(d,predict(exp), col="red")
# lines(d,(predict(lo)+predict(exp))/2, col="green")
 lines(d,predict(exp))
 lines(d,predict(lo),col="green")
 lines(d[round(length(a$file)/2):length(a$file)],predict(expend),col="red")

 graphics.off()
 i <<- i+1
})

i <- 0
png(paste(name,"-TPVCEs.png",sep=""), width=800, height=800)
plot(1,xlim=xlim,ylim=ylim,cex.main=1.5,xlab="Extracted Charge Density (C/cm2)", ylab="Life-time (s)", main=paste(name,"TPV decay vs Charge from CE"),cex.lab=1.5,cex.axis=1.5,log="y");
lapply(dirs, function(x) {print(x);
 a <- read.table(paste(x,"-outputChargeDensityCE.txt",sep=""),header=T,stringsAsFactors=F)
 b<-strsplit(a$file, "_")
 c<-unlist(b)[length(b[[1]])*(1:length(a$file))]
 d<-as.numeric(gsub("mV", "", c))
 lo <- loess(a$ChargeDensityCE~d,span=0.9)
 a$d <- d
 exp <- nlrob(ChargeDensityCE~ A+C*exp(D*d), start=list(A=0,C=2e-9,D=9), data=a)
 expend <- nlsLM(ChargeDensityCE~ A+C*exp(D*d), start=list(A=coef(exp)["A"],C=coef(exp)["C"],D=coef(exp)["D"]), data=a[round(length(a$file)/2):length(a$file),])
 # sinexp <- nls(ChargeDensityCE~ cbind(1,d,sin(D*d),exp(F*d)), start=c(D=-7,F=6), data=a, alg="plinear")
# sinexp2 <- nlrob(ChargeDensityCE~ A+B*d+C*sin(D*d)+E*exp(F*d), start=list(A=coef(sinexp)[".lin1"],B=coef(sinexp)[".lin.d"],C=coef(sinexp)[".lin3"],D=coef(sinexp)["D"],E=coef(sinexp)[".lin4"],F=coef(sinexp)["F"]), data=a)

fulloutput <- read.table(paste(x,"/output-monoexp.txt",sep=""), header=TRUE);
n<-tail(grep("file",fulloutput[,1]),n=1)
tpv <- read.table(paste(x,"/output-monoexp.txt",sep=""), header=TRUE, skip=ifelse(length(n),n,0));
#importante che la variabile in new abbia lo stesso nome di quella fittata
new <- data.frame(d = tpv$Voc)
charge <- (predict(lo,tpv$Voc)+predict(exp,new))/2
new2 <- data.frame(d = tpv$Voc[is.na(charge)])
charge[is.na(charge)] <- (predict(exp,new2) + predict(expend,new2))/2
points(charge, tpv$T, lwd=2, pch=i, col=colors[i+1]);
 i <<- i+1
})
legend(x="topright",inset=0.05,dirs,pch=seq(0,10,1), col=colors,cex=1.5)
graphics.off()
