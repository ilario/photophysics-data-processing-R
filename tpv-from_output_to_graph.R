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

#library(Hmisc)
print("errori derivano da avere molti header in output, bisogna pulirlo")
directory <- tail(strsplit(getwd(), "/")[[1]], n=2)
name<- directory[1]

print("biexp")
fulloutput <- read.table("output-biexp.txt", header=TRUE);#,stringsAsFactors=F);
n<-tail(grep("file",fulloutput[,1]),n=1)
#output <- fulloutput[(ifelse(length(n),n,0)+1):nrow(fulloutput),]
tpv <- read.table("output-biexp.txt", header=TRUE, skip=ifelse(length(n),n,0)); 
png(paste(name, "-tpv-biexp.png",sep=""), width=600, height=600);
plot(1, ylim=c(min(tpv$T1,tpv$T2),max(tpv$T2,tpv$T1)), xlim=c(min(tpv$Voc),max(tpv$Voc)), log="y", xlab="Voc (V)", ylab="Life-time (s)", main=paste(name, "TPV biexp"), cex.axis=1.5, cex.lab=1.5);
#errbar(tpv$Voc, tpv$T1, tpv$T1+tpv$T1.error, tpv$T1-tpv$T1.error,log="y",add=TRUE, pch="")
#errbar(tpv$Voc, tpv$T2, tpv$T2+tpv$T2.error, tpv$T2-tpv$T2.error,log="y",add=TRUE, pch="")
points(tpv$Voc, tpv$T1, pch=20, col="red")#, cex=4*tpv$A1/(tpv$A1+tpv$A2));
points(tpv$Voc, tpv$T2, pch=20, col="black")#, cex=4*tpv$A2/(tpv$A1+tpv$A2));
graphics.off()


#print("robust biexp")
#
#fulloutput <- read.table("output-robustbiexp.txt", header=TRUE);
#n<-tail(grep("file",fulloutput[,1]),n=1)
#tpv <- read.table("output-robustbiexp.txt", header=TRUE, skip=ifelse(length(n),n,0)); 
#png(paste(name, "-tpv-robustbiexp.png",sep=""), width=1280, height=800);
#plot(1, ylim=c(min(tpv$T1,tpv$T2),max(tpv$T2,tpv$T1)), xlim=c(min(tpv$Voc),max(tpv$Voc)), log="y", xlab="Voc (V)", ylab="Life-time (s)", main=paste(name, "TPV robust biexp"));
#errbar(tpv$Voc, tpv$T1, tpv$T1+tpv$T1.error, tpv$T1-tpv$T1.error,log="y",add=TRUE)
#errbar(tpv$Voc, tpv$T2, tpv$T2+tpv$T2.error, tpv$T2-tpv$T2.error,log="y",add=TRUE)
#points(tpv$Voc, tpv$T1, pch=21, col="red", bg="red");
#points(tpv$Voc, tpv$T2, col="red");
#graphics.off()
#
print("monoexp")

fulloutput <- read.table("output-monoexp.txt", header=TRUE);
n<-tail(grep("file",fulloutput[,1]),n=1)
tpv <- read.table("output-monoexp.txt", header=TRUE, skip=ifelse(length(n),n,0)); 
lo<-loess(tpv$T~tpv$Voc,span=0.5)

png(paste(name, "-tpv-monoexp.png",sep=""), width=400, height=400);
par(mar=c(5.1,5,4.1,2.1))
plot(tpv$Voc, tpv$T, col="black", log="y", xlab=bquote("V"["oc"]~"(V)"), ylab="Life-time (s)", cex.axis=1.4, cex.lab=1.4)#, main=paste(name, "TPV monoexp"));
#errbar(tpv$Voc, tpv$T, tpv$T+tpv$T.error, tpv$T-tpv$T.error,log="y",add=TRUE)
lines(tpv$Voc, predict(lo), lwd=1,col="red")
graphics.off()
#
#print("robust monoexp")
#
#fulloutput <- read.table("output-robustmonoexp.txt", header=TRUE);
#n<-tail(grep("file",fulloutput[,1]),n=1)
#tpv <- read.table("output-robustmonoexp.txt", header=TRUE, skip=ifelse(length(n),n,0)); 
#png(paste(name, "-tpv-robustmonoexp.png",sep=""), width=1280, height=800);
#plot(tpv$Voc, tpv$T, col="red", log="y", xlab="Voc (V)", ylab="Life-time (s)", main=paste(name, "TPV robust monoexp"));
#errbar(tpv$Voc, tpv$T, tpv$T+tpv$T.error, tpv$T-tpv$T.error,log="y",add=TRUE)
#graphics.off()


print("mixed biexp monoexp")
fulloutput <- read.table("output-mixedbimono.txt", header=TRUE, fill = TRUE);#,stringsAsFactors=F);
n<-tail(grep("file",fulloutput[,1]),n=1)
#tpv <- fulloutput[(ifelse(length(n),n,0)+1):nrow(fulloutput),]
tpv <- read.table("output-mixedbimono.txt", header=TRUE, skip=ifelse(length(n),n,0), fill=TRUE); 
png(paste(name, "-tpv-mixedbimono.png",sep=""), width=600, height=600);
plot(1, ylim=c(min(tpv$T1,tpv$T2[!is.na(tpv$T2)]),max(tpv$T2[!is.na(tpv$T2)],tpv$T1)), xlim=c(min(tpv$Voc),max(tpv$Voc)), log="y", xlab="Voc (V)", ylab="Life-time (s)", main=paste(name, "TPV biexp and monoexp"), cex.axis=1.5, cex.lab=1.5);
#errbar(tpv$Voc, tpv$T1, tpv$T1+tpv$T1.error, tpv$T1-tpv$T1.error,log="y",add=TRUE, pch="")
#errbar(tpv$Voc, tpv$T2, tpv$T2+tpv$T2.error, tpv$T2-tpv$T2.error,log="y",add=TRUE, pch="")
points(tpv$Voc, tpv$T1, pch=20, col="red")#, cex=4*tpv$A1/(tpv$A1+tpv$A2));
points(tpv$Voc, tpv$T2, pch=20, col="black")#, cex=4*tpv$A2/(tpv$A1+tpv$A2));
graphics.off()

