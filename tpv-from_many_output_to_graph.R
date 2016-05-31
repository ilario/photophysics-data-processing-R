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

#for i in `ls`; do cp ../$i/tpv/output* $i; done

#name=""

ylim=c(3e-7,6e-4)
xlim=c(0,1)

library(Hmisc)
library(RColorBrewer)

dirs <- list.dirs(recursive=FALSE)
dirs <- sub("./","",dirs)
colors=brewer.pal(max(length(dirs),3),"Set1")

print("errori derivano da avere molti header in output, bisogna pulirlo")
#directory <- tail(strsplit(getwd(), "/")[[1]], n=2)
#name<- directory[1]

print("biexp")

i <- 0
png(paste(name, "-TPVs-biexp.png",sep=""), width=800, height=800);
plot(1, ylim=ylim, xlim=xlim, cex.main=1.5,cex.axis=1.5,cex.lab=1.5, log="y", xlab="Voc (V)", ylab="Life-time (s)", main=paste(name, "TPV biexp"));
lapply(dirs, function(x) {print(x);
fulloutput <- read.table(paste(x,"/output-biexp.txt",sep=""), header=TRUE);#,stringsAsFactors=F);
n<-tail(grep("file",fulloutput[,1]),n=1)
output <- read.table(paste(x,"/output-biexp.txt",sep=""), header=TRUE, skip=ifelse(length(n),n,0)); 
#errbar(output$Voc, output$T1, output$T1+output$T1.error, output$T1-output$T1.error,log="y",add=TRUE)
#errbar(output$Voc, output$T2, output$T2+output$T2.error, output$T2-output$T2.error,log="y",add=TRUE)
points(output$Voc, output$T1, lwd=2, pch=i+15, col=colors[i+1], bg=colors[i+1]);
points(output$Voc, output$T2, lwd=2, pch=i, col=colors[i+1]);
i <<- i+1
})
legend(x="bottomleft",inset=0.05,dirs,pt.cex=2,cex=1.5,pch=seq(0,10,1), col=colors)
graphics.off()


#print("robust biexp")
#i <- 0
#png(paste(name, "-TPVs-robustbiexp.png",sep=""), width=800, height=800);
#plot(1, ylim=ylim, xlim=xlim, cex.main=1.5,cex.axis=1.5,cex.lab=1.5, log="y", xlab="Voc (V)", ylab="Life-time (s)", main=paste(name, "TPV robust biexp"));
#lapply(dirs, function(x) {print(x);
#fulloutput <- read.table(paste(x,"/output-robustbiexp.txt",sep=""), header=TRUE);
#n<-tail(grep("file",fulloutput[,1]),n=1)
#output <- read.table(paste(x,"/output-robustbiexp.txt",sep=""), header=TRUE, skip=ifelse(length(n),n,0)); 
##errbar(output$Voc, output$T1, output$T1+output$T1.error, output$T1-output$T1.error,log="y",add=TRUE)
##errbar(output$Voc, output$T2, output$T2+output$T2.error, output$T2-output$T2.error,log="y",add=TRUE)
#points(output$Voc, output$T1, pch=21, col=colors[i+1], bg=colors[i+1]);
#points(output$Voc, output$T2, col=colors[i+1]);
#i <<- i+1
#})
#legend(x="bottomleft",inset=0.05,dirs,pt.cex=2,cex=1.5,pch=seq(0,10,1), col=colors)
#graphics.off()

print("monoexp")
i <- 0
png(paste(name, "-TPVs-monoexp.png",sep=""), width=800, height=800);
plot(1, ylim=ylim, xlim=xlim, cex.main=1.5,cex.axis=1.5,cex.lab=1.5, log="y", xlab="Voc (V)", ylab="Life-time (s)", main=paste(name, "TPV monoexp"));
lapply(dirs, function(x) {print(x);
fulloutput <- read.table(paste(x,"/output-monoexp.txt",sep=""), header=TRUE);
n<-tail(grep("file",fulloutput[,1]),n=1)
output <- read.table(paste(x,"/output-monoexp.txt",sep=""), header=TRUE, skip=ifelse(length(n),n,0)); 
#errbar(output$Voc, output$T, output$T+output$T.error, output$T-output$T.error,log="y",add=TRUE)
points(output$Voc, output$T, lwd=2, pch=i, col=colors[i+1]);
i <<- i+1
})
legend(x="bottomleft",inset=0.05,dirs,pt.cex=2,cex=1.5,pch=seq(0,10,1), col=colors)

graphics.off()

#print("robust monoexp")
#i <- 0
#png(paste(name, "-TPVs-robustmonoexp.png",sep=""), width=800, height=800);
#plot(1, ylim=ylim, xlim=xlim, cex.main=1.5,cex.axis=1.5,cex.lab=1.5, log="y", xlab="Voc (V)", ylab="Life-time (s)", main=paste(name, "TPV robust monoexp"));
#lapply(dirs, function(x) {print(x);
#fulloutput <- read.table(paste(x,"/output-robustmonoexp.txt",sep=""), header=TRUE);
#n<-tail(grep("file",fulloutput[,1]),n=1)
#output <- read.table(paste(x,"/output-robustmonoexp.txt",sep=""), header=TRUE, skip=ifelse(length(n),n,0)); 
##errbar(output$Voc, output$T, output$T+output$T.error, output$T-output$T.error,log="y",add=TRUE)
#points(output$Voc, output$T, lwd=2, pch=i, col=colors[i+1]);
#i <<- i+1
#})
#legend(x="bottomleft",inset=0.05,dirs,pt.cex=2,cex=1.5,pch=seq(0,10,1), col=colors)
#graphics.off()
