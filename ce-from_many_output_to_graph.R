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
library(minpack.lm)
i <- 0
#files <- list.files(path=".", pattern="*-outputChargeDensityCE.txt$")
dirs <- list.dirs(recursive=FALSE)
dirs <- sub("./","",dirs)
colors=brewer.pal(max(length(dirs),3),"Set1")
#colors=brewer.pal(length(files),"Set1")

# e<-1.6021766208e-19
png(paste(name,"-CEs.png",sep=""), width=800, height=640)
par(mar=c(5.1,5,4.1,2.1))
plot(NULL,xlim=c(0,1),ylim=c(0,2e-7),cex.main=1.5,xlab="Voltage (V)",ylab="Extracted Charge Density (C/cm2)", main=paste(name,"CEs"), cex.lab=1.5, cex.axis=1.5);

lapply(dirs, function(x) {print(x);
 a <- read.table(paste(x,"/outputChargeDensityCE.txt",sep=""),header=T,stringsAsFactors=F)
 b<-strsplit(a$file, "_")
 c<-unlist(b)[length(b[[1]])*(1:length(a$file))]
 d<-as.numeric(gsub("mV", "", c))


a$d <- d
#print(a)
#lin <- lmrob(ChargeDensityCE~ d, data=a)
#print(list(A=0,C=(2e-9/e),D=3))
exp <- nlrob(ChargeDensityCE~ A+C*exp(D*d), start=list(A=0,#A=coef(lin)[1],
#B=coef(lin)[2],
C=2e-9,D=9), data=a)
 
 points(d, a$ChargeDensityCE, lwd=2, pch=i, col=colors[i+1])
f <- data.frame(d = sort(d))
 lines(f$d,predict(exp,f),col=colors[i+1])
# mtext(bquote(.(gsub("-outputChargeDensityCE.txt","",x))~": n" == .(signif(exp$coefficients["A"],3)) + 
#	      .(signif(exp$coefficients["C"],3)) ~ "e" ^ {.(signif(exp$coefficients["D"],3))~V}),side=3,line=-(i*2+4),cex=1.5,col=colors[i+1])
 mtext(bquote(.(x)~": n" == .(signif(exp$coefficients["A"],3)) + 
	      .(signif(exp$coefficients["C"],3)) ~ "e" ^ {.(signif(exp$coefficients["D"],3))~V}),side=3,line=-(i*2+4),cex=1.5,col=colors[i+1])
 i <<- i+1
})

legend(x="topleft",inset=0.05,dirs,pch=seq(0,10,1), col=colors, pt.cex=2, cex=1.5)

graphics.off()
