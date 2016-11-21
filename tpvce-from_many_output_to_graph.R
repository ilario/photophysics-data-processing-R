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

title=gsub("-","\n\n",gsub("_"," ",name))
filename=gsub(",","",gsub(":","",name))

ylim=limlifetime
xlim=limcecharge

library(robustbase)
library(RColorBrewer)
library(minpack.lm)
library(sfsmisc)
library(Hmisc)

i <- 0
dirs <- list.dirs(recursive=FALSE)
colors=colorRampPalette(c("red","orange","springgreen","royalblue"))(max(length(dirs),3))
#brewer.pal(max(length(dirs),3),"Spectral")
dirs <- sub("./","",dirs)

#lapply(dirs, function(x) {print(x);
# a <- read.table(paste(x,"/outputChargeDensityCE.txt",sep=""),header=T,stringsAsFactors=F)
# b<-strsplit(a$file, "_")
# c<-unlist(b)[length(b[[1]])*(1:length(a$file))]
# d<-as.numeric(gsub("mV", "", c))
# lo <- loess(a$ChargeDensityCE~d,span=0.9)
# a$d <- d
#       	exp <- nlrob(ChargeDensityCE~ A+C*exp(D*d), start=list(A=0,C=2e-9,D=9), data=a)
# expend <- nlsLM(ChargeDensityCE~ A+C*exp(D*d), start=list(A=coef(exp)["A"],C=coef(exp)["C"],D=coef(exp)["D"]), data=a[round(length(a$file)/2):length(a$file),])
# png(paste(x,"-CEs.png",sep=""), width=640, height=640)
# plot(NULL,xlim=c(0,1),ylim=c(0,2e-7),cex.main=1.5,xlab="Voltage (V)",ylab="Extracted Charge Density (C/cm2)", main=paste(x,"CE fitted"));
# points(d, a$ChargeDensityCE, lwd=2, pch=i, col=colors[i+1])
# lines(d,predict(exp))
# lines(d,predict(lo),col="green")
# lines(d[round(length(a$file)/2):length(a$file)],predict(expend),col="red")
# graphics.off()
# i <<- i+1
#})

i <- 0
png(paste(filename,"-TPVCEs.png",sep=""), width=640, height=640)
par(mar=c(5.1,5,2,2.1))
plot(1,xlim=xlim,ylim=ylim,cex.main=1.5,xlab=bquote("Extracted Charge Density (C/cm"^"2"*")"), ylab="Life-time (s)",cex.lab=1.5,cex.axis=1.2,log="y", yaxt="n", xaxt="n")#, main=paste(name,"TPV decay vs Charge from CE");
eaxis(side=2,at=c(1e-10,1e-9,1e-8,1e-7,1e-6,1e-5,1e-4,1e-3,1e-2,0.1,1,10,100,1e3), cex.axis=1.2)
eaxis(side=1, cex.axis=1.2)
minor.tick(nx=10)

lapply(dirs, function(x) {print(x);
 a <- read.table(paste(x,"/ce/outputChargeDensityCE.txt",sep=""),header=T,stringsAsFactors=F)
 b<-strsplit(a$file, "_")
 c<-unlist(b)[length(b[[1]])*(1:length(a$file))]
 d<-as.numeric(gsub("mV", "", c))
 lo <- loess(a$ChargeDensityCE~d,span=0.9)
 a$d <- d
#exp <- nlrob(ChargeDensityCE~ A+C*exp(D*d), start=list(A=0,C=1e-10,D=9), data=a)
#expend <- nlsLM(ChargeDensityCE~ A+C*exp(D*d), start=list(A=coef(exp)["A"],C=coef(exp)["C"],D=coef(exp)["D"]), data=a[round(length(a$file)/2):length(a$file),])
 expend <- nlsLM(ChargeDensityCE~ C*(exp(D*d)-A), start=list(A=1,C=1e-10,D=9), data=a[round(length(a$file)/2):length(a$file),])
 tryCatch({
	  exp <- nlrob(ChargeDensityCE~ C*(exp(D*d)-A), start=list(A=coef(expend)["A"],C=coef(expend)["C"],D=coef(expend)["D"]), data=a)
 }, error=function(e) {print("FAILED ZEROTH FIT")});
 tryCatch({
	  exp <- nlrob(ChargeDensityCE~ B*d+C*(exp(D*d)-1), start=list(B=1e-9,C=coef(expend)["C"],D=coef(expend)["D"]), data=a)
 }, error=function(e) {print("FAILED FIRST FIT")});
 tryCatch({
	  exp <- nlrob(ChargeDensityCE~ B*d+C*(exp(D*d)-1), start=list(B=2e-8,C=1e-8,D=0.1), data=a)
 }, error=function(e) {print("FAILED SECOND FIT")});
 tryCatch({
	  exp <- nlrob(ChargeDensityCE~ B*d+C*(exp(D*d)-1), start=list(B=1e-9,C=1e-10,D=8), data=a)
 }, error=function(e) {print("FAILED THIRD FIT")});
 tryCatch({
	  exp <- nlrob(ChargeDensityCE~ B*d+C*(exp(D*d)-A), start=list(A=coef(expend)["A"],B=1e-9,C=coef(expend)["C"],D=coef(expend)["D"]), data=a)
 }, error=function(e) {print("FAILED FOURTH FIT")});
 tryCatch({
	  exp <- nlrob(ChargeDensityCE~ B*d+C*(exp(D*d)-A), start=list(A=1,B=1e-9,C=1e-10,D=8), data=a)
 }, error=function(e) {print("FAILED FIFTH FIT")});

fulloutput <- read.table(paste(x,"/tpv/output-monoexp.txt",sep=""), header=TRUE);
n<-tail(grep("file",fulloutput[,1]),n=1)
tpv <- read.table(paste(x,"/tpv/output-monoexp.txt",sep=""), header=TRUE, skip=ifelse(length(n),n,0));
#importante che la variabile in new abbia lo stesso nome di quella fittata
new <- data.frame(d = tpv$Voc)
charge <- (predict(lo,tpv$Voc)+predict(exp,new))/2
new2 <- data.frame(d = tpv$Voc[is.na(charge)])
charge[is.na(charge)] <- (predict(exp,new2) + predict(expend,new2))/2
lo<-loess(tpv$T~charge,span=0.5)
lines(charge, predict(lo), lwd=2, col=colors[i+1])
points(charge, tpv$T, lwd=1, bg=colors[i+1], cex=2, pch=21+i);
 i <<- i+1
})
legend(x="topright",inset=0.1,sub("-ig..-...-.","",sub("^0","",dirs)),pch=seq(21,25), pt.bg=colors, lwd=4, pt.lwd=2, pt.cex=2, col=colors,cex=1.5, title=#paste("TPV vs CE\n",
title,bg="gray90"#, bty="n")
	)
graphics.off()
