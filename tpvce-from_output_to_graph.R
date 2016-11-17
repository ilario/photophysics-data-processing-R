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

tpvCeFromOutputToGraph <- function(cedir="ce", tpvdir="tpv", printcefit=FALSE)
{
#name=""
directory <- tail(strsplit(getwd(), "/")[[1]], n=1)
name<- directory

library(robustbase)
library(RColorBrewer)
library(minpack.lm)
write.table(t(c("ChargeDensityCE", "T")), file="outputTPVCE-monoexp.txt", append=FALSE, col.names=F, row.names=F);
write.table(t(c("ChargeDensityCE", "T1", "T2")), file="outputTPVCE-biexp.txt", append=FALSE, col.names=F, row.names=F);
if(file.exists(file.path(tpvdir, "output-mixedbimono.txt")))
{
	write.table(t(c("ChargeDensityCE", "T1", "T2")), file="outputTPVCE-mixedbimono.txt", append=FALSE, col.names=F, row.names=F);
}

 a <- read.table(file.path(cedir, "outputChargeDensityCE.txt"), header=T,stringsAsFactors=F)
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

tryCatch({
	 exp <- nlrob(ChargeDensityCE~ A+C*exp(D*d), start=list(A=0,C=1e-10,D=9), data=a)
         }, error=function(e) {print("FAILED FOURTH FIT")});

 expend <- nlsLM(ChargeDensityCE~ A+C*exp(D*d), start=list(A=0,C=coef(exp)["C"],D=coef(exp)["D"]), data=a[round(length(a$file)/2):length(a$file),])

if(printcefit)
{
 png(paste(name,"-CE-fit.png",sep=""), width=800, height=640)
 plot(d, a$ChargeDensityCE, cex.main=1.5,xlab="Voltage (V)",ylab="Extracted Charge Density (C/cm2)", main=paste(name,"CEs"), lwd=2)
 lines(d,predict(exp))
 lines(d,predict(lo),col="green")
 lines(d[round(length(a$file)/2):length(a$file)],predict(expend),col="red")
 graphics.off()
}

fulloutput <- read.table(file.path(tpvdir, "output-monoexp.txt"), header=TRUE);
n<-tail(grep("file",fulloutput[,1]),n=1)
tpv <- read.table(file.path(tpvdir, "output-monoexp.txt"), header=TRUE, skip=ifelse(length(n),n,0));
#importante che la variabile in new abbia lo stesso nome di quella fittata
new <- data.frame(d = tpv$Voc)
charge <- (predict(lo,tpv$Voc)+predict(exp,new))/2
new2 <- data.frame(d = tpv$Voc[is.na(charge)])
charge[is.na(charge)] <- (predict(exp,new2) + predict(expend,new2))/2

outputTPVCEmonoexp <- data.frame(charge, tpv$T)
write.table(outputTPVCEmonoexp, file="outputTPVCE-monoexp.txt", append=TRUE, col.names=F, row.names=F, quote=F);

png(paste(name,"-tpvce-monoexp.png",sep=""), width=600, height=600)
plot(charge, tpv$T,cex.main=1.5,xlab="Extracted Charge Density (C/cm2)", ylab="Life-time (s)", main=paste(name,"TPV decay vs Charge from CE"),cex.lab=1.5,cex.axis=1.5,log="y", lwd=2);
#legend(x="topright",inset=0.05,dirs,pch=seq(0,10,1), col=colors,cex=1.5)
graphics.off()


fulloutput <- read.table(file.path(tpvdir, "output-biexp.txt"), header=TRUE)#, fill = TRUE);
n<-tail(grep("file",fulloutput[,1]),n=1)
#output <- fulloutput[(ifelse(length(n),n,0)+1):nrow(fulloutput),]
tpv <- read.table(file.path(tpvdir, "output-biexp.txt"), header=TRUE, skip=ifelse(length(n),n,0))#, fill=TRUE);
new <- data.frame(d = tpv$Voc)
charge <- (predict(lo,tpv$Voc)+predict(exp,new))/2
new2 <- data.frame(d = tpv$Voc[is.na(charge)])
charge[is.na(charge)] <- (predict(exp,new2) + predict(expend,new2))/2

outputTPVCEbiexp <- data.frame(charge, tpv$T1, tpv$T2)
write.table(outputTPVCEbiexp, file="outputTPVCE-biexp.txt", append=TRUE, col.names=F, row.names=F, quote=F);

charge2 <- charge*1e9
T1 <- tpv$T1*1e6
T2 <- tpv$T2*1e6

png(paste(name, "-tpvce-biexp.png",sep=""), width=600, height=600);
par(mar=c(5.1,4.5,4.1,2.1))
plot(0, ylim=c(min(T1,T2),max(T2,T1)), xlim=c(min(charge2),max(charge2)), log="y", xlab=bquote("Extracted Charge Density (nC/cm"^"2"*")"), ylab=bquote("Life-time ("*mu*"s)"), main=paste(name, "TPV biexp and monoexp decay vs Charge from CE"), cex.axis=1.5, cex.lab=1.5);
#errbar(tpv$Voc, tpv$T1, tpv$T1+tpv$T1.error, tpv$T1-tpv$T1.error,log="y",add=TRUE, pch="")
#errbar(tpv$Voc, tpv$T2, tpv$T2+tpv$T2.error, tpv$T2-tpv$T2.error,log="y",add=TRUE, pch="")
points(charge2, T1, pch=20, col="red")#, cex=4*tpv$A1/(tpv$A1+tpv$A2));
points(charge2, T2, pch=20, col="black")#, cex=4*tpv$A2/(tpv$A1+tpv$A2));
graphics.off()

if(file.exists(file.path(tpvdir, "output-mixedbimono.txt")))
{
fulloutput <- read.table(file.path(tpvdir, "output-mixedbimono.txt"), header=TRUE, fill = TRUE);
n<-tail(grep("file",fulloutput[,1]),n=1)
#output <- fulloutput[(ifelse(length(n),n,0)+1):nrow(fulloutput),]
tpv <- read.table(file.path(tpvdir, "output-mixedbimono.txt"), header=TRUE, skip=ifelse(length(n),n,0), fill=TRUE);
new <- data.frame(d = tpv$Voc)
charge <- (predict(lo,tpv$Voc)+predict(exp,new))/2
new2 <- data.frame(d = tpv$Voc[is.na(charge)])
charge[is.na(charge)] <- (predict(exp,new2) + predict(expend,new2))/2

outputTPVCEmixedbimono <- data.frame(charge, tpv$T1, tpv$T2)
write.table(outputTPVCEmixedbimono, file="outputTPVCE-mixedbimono.txt", append=TRUE, col.names=F, row.names=F, quote=F, na="");

charge2 <- charge*1e9
T1 <- tpv$T1*1e6
T2 <- tpv$T2*1e6
png(paste(name, "-tpvce-mixedbimono.png",sep=""), width=600, height=600);
par(mar=c(5.1,4.5,4.1,2.1))
plot(0, ylim=c(min(T1,T2[!is.na(tpv$T2)]),max(T2[!is.na(tpv$T2)],T1)), xlim=c(min(charge2),max(charge2)), log="y", xlab=bquote("Extracted Charge Density (nC/cm"^"2"*")"), ylab=bquote("Life-time ("*mu*"s)"), main=paste(name, "TPV biexp and monoexp decay vs Charge from CE"), cex.axis=1.5, cex.lab=1.5);
#errbar(tpv$Voc, tpv$T1, tpv$T1+tpv$T1.error, tpv$T1-tpv$T1.error,log="y",add=TRUE, pch="")
#errbar(tpv$Voc, tpv$T2, tpv$T2+tpv$T2.error, tpv$T2-tpv$T2.error,log="y",add=TRUE, pch="")
points(charge2, T1, pch=20, col="red")#, cex=4*tpv$A1/(tpv$A1+tpv$A2));
points(charge2[!is.na(tpv$T2)], T2[!is.na(tpv$T2)], pch=20, col="black")#, cex=4*tpv$A2/(tpv$A1+tpv$A2));
graphics.off()
}
}
