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

library(robustbase)
library(RColorBrewer)
library(minpack.lm)

 a <- read.table("ce/outputChargeDensityCE.txt",header=T,stringsAsFactors=F)
 b<-strsplit(a$file, "_")
 c<-unlist(b)[length(b[[1]])*(1:length(a$file))]
 d<-as.numeric(gsub("mV", "", c))
 lo <- loess(a$ChargeDensityCE~d,span=0.9)
 a$d <- d

 tryCatch({
	  exp <- nlrob(ChargeDensityCE~ C*exp(D*d), start=list(C=2e-8,D=5), data=a)
	   }, error=function(e) {print("FAILED ZEROTH FIT")});
 
 tryCatch({
 exp <- nlrob(ChargeDensityCE~ A+C*exp(D*d), start=list(A=0,C=2e-8,D=0.1), data=a)
 }, error=function(e) {print("FAILED FIRST FIT")});

tryCatch({
	exp <- nlrob(ChargeDensityCE~ A+C*exp(D*d), start=list(A=0,C=2e-10,D=1), data=a)
	 }, error=function(e) {print("FAILED SECOND FIT")});

tryCatch({
	 exp <- nlrob(ChargeDensityCE~ A+C*exp(D*d), start=list(A=0,C=2e-9,D=9), data=a)
         }, error=function(e) {print("FAILED THIRD FIT")});


 expend <- nlsLM(ChargeDensityCE~ A+C*exp(D*d), start=list(A=0,C=coef(exp)["C"],D=coef(exp)["D"]), data=a[round(length(a$file)/2):length(a$file),])
 png(paste(name,"-CE-fit.png",sep=""), width=800, height=640)
 plot(d, a$ChargeDensityCE, cex.main=1.5,xlab="Voltage (V)",ylab="Extracted Charge Density (C/cm2)", main=paste(name,"CEs"), lwd=2)
 lines(d,predict(exp))
 lines(d,predict(lo),col="green")
 lines(d[round(length(a$file)/2):length(a$file)],predict(expend),col="red")
 graphics.off()

# a <- read.table("outputChargeDensityCE.txt",header=T,stringsAsFactors=F)
# b<-strsplit(a$file, "_")
# c<-unlist(b)[length(b[[1]])*(1:length(a$file))]
# d<-as.numeric(gsub("mV", "", c))
# lo <- loess(a$ChargeDensityCE~d,span=0.9)
# a$d <- d
# exp <- nlrob(ChargeDensityCE~ A+C*exp(D*d), start=list(A=0,C=2e-9,D=9), data=a)
# expend <- nlsLM(ChargeDensityCE~ A+C*exp(D*d), start=list(A=coef(exp)["A"],C=coef(exp)["C"],D=coef(exp)["D"]), data=a[round(length(a$file)/2):length(a$file),])

fulloutput <- read.table("tpv/output-monoexp.txt", header=TRUE);
n<-tail(grep("file",fulloutput[,1]),n=1)
tpv <- read.table("tpv/output-monoexp.txt", header=TRUE, skip=ifelse(length(n),n,0));
#importante che la variabile in new abbia lo stesso nome di quella fittata
new <- data.frame(d = tpv$Voc)
charge <- (predict(lo,tpv$Voc)+predict(exp,new))/2
new2 <- data.frame(d = tpv$Voc[is.na(charge)])
charge[is.na(charge)] <- (predict(exp,new2) + predict(expend,new2))/2
png(paste(name,"-TPVCEs.png",sep=""), width=800, height=800)
plot(charge, tpv$T,cex.main=1.5,xlab="Extracted Charge Density (C/cm2)", ylab="Life-time (s)", main=paste(name,"TPV decay vs Charge from CE"),cex.lab=1.5,cex.axis=1.5,log="y", lwd=2);
#legend(x="topright",inset=0.05,dirs,pch=seq(0,10,1), col=colors,cex=1.5)
graphics.off()
