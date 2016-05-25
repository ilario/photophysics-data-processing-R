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

a <- read.table("outputChargeDensityCE.txt",header=T,stringsAsFactors=F)

#b<-list.files(path=".", pattern="^CE.*\\.txt.table$");
b<-strsplit(a$file, "_")
c<-unlist(b)[length(b[[1]])*(1:length(a$file))]
d<-as.numeric(gsub("mV", "", c))

directory <- tail(strsplit(getwd(), "/")[[1]], n=2)
png(paste("charge_extraction-", directory[1], ".png", sep=""), width=1000, heigh=600)
plot(d, a$ChargeDensityCE, ylab="Extracted Charge Density (C/cm2)", xlab="Voltage (V)", main=paste(directory,collapse=" "))
graphics.off()
