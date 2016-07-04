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

library(robustbase)
a <- read.table("outputChargeDensityCE.txt",header=T,stringsAsFactors=F)

#b<-list.files(path=".", pattern="^CE.*\\.txt.table$");
b<-strsplit(a$file, "_")
c<-unlist(b)[length(b[[1]])*(1:length(a$file))]
d<-as.numeric(gsub("mV", "", c))

a$d <- d
exp <- nlrob(ChargeDensityCE~ A+C*exp(D*d), start=list(A=0,C=2e-9,D=9), data=a)
f <- data.frame(d = sort(d))

directory <- tail(strsplit(getwd(), "/")[[1]], n=2)
png(paste("charge_extraction-", directory[1], ".png", sep=""), width=800, heigh=800)
plot(d, a$ChargeDensityCE, ylab="Charge Density (C/cm2)", xlab="Voltage (V)",cex.lab=1.4, cex.axis=1.4)#, main=paste(directory,collapse=" "))
lines(f$d,predict(exp,f), lwd=1, col="red")
graphics.off()
