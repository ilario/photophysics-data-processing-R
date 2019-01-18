# Copyright (C) 2017 Ilario Gelmetti <iochesonome@gmail.com>
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

tpvFromOutputToGraphWithLimits <- function(tpvdir="tpv", dcdir=".")
{
  print("TPV with limits: PLOTTING")
  require(minpack.lm)
  require(robustbase)
  directory <- tail(strsplit(getwd(), "/")[[1]], n=2)
  name<- directory[1]
  impedance = 1e6
  
  fulloutput <- read.table(file.path(tpvdir,"output-biexp.txt"), header=TRUE);#,stringsAsFactors=F);
  n<-tail(grep("file",fulloutput[,1]),n=1)
  tpv <- read.table(file.path(tpvdir,"output-biexp.txt"), header=TRUE, skip=ifelse(length(n),n,0)); 
  
  DCcapacitance <- read.table(file.path(dcdir,"outputDCcapacitance.txt"), header=TRUE);
  DCcapacitance <- DCcapacitance[with(DCcapacitance, order(Voc)), ]
  RCtime <- impedance*DCcapacitance$capacitance*cellArea
  
  png(file.path(tpvdir, paste(name, "-tpv-biexp-limits.png",sep="")), width=600, height=600);
  plot(1, ylim=c(min(tpv$T1,tpv$T2), max(tpv$T2, tpv$T1, RCtime)), xlim=c(min(tpv$Voc),max(tpv$Voc)), log="y", xlab=bquote("V"["OC"]*"(V)"), ylab="Life-time (s)", main=paste(name, "TPV biexp"), cex.axis=1.5, cex.lab=1.5);
  points(tpv$Voc, tpv$T1, pch=20, col="red")#, cex=4*tpv$A1/(tpv$A1+tpv$A2));
  points(tpv$Voc, tpv$T2, pch=20, col="black")#, cex=4*tpv$A2/(tpv$A1+tpv$A2));
  
  abline(h=1e-7)
  
  lines(DCcapacitance$Voc, RCtime)
  
  graphics.off()
  
  
  if(file.exists(file.path(tpvdir, "output-mixedbimono.txt")))
  {
    fulloutput <- read.table(file.path(tpvdir,"output-mixedbimono.txt"), header=TRUE, fill = TRUE);
    n<-tail(grep("file",fulloutput[,1]),n=1)
    tpv <- read.table(file.path(tpvdir,"output-mixedbimono.txt"), header=TRUE, skip=ifelse(length(n),n,0), fill=TRUE); 
    png(file.path(tpvdir, paste(name, "-tpv-mixedbimono-limits.png",sep="")), width=600, height=600);
    par(mar=c(5.1, 4.1+1, 4.1, 2.1))
    plot(1, ylim=c(min(tpv$T1,tpv$T2[!is.na(tpv$T2)]),max(tpv$T2[!is.na(tpv$T2)],tpv$T1, RCtime)), xlim=c(min(tpv$Voc),max(tpv$Voc)), log="y", xlab=bquote("V"["OC"]*~"(V)"), ylab="Life-time (s)", cex.axis=2, cex.lab=2);#main=paste(name, "TPV biexp and monoexp"), 
    points(tpv$Voc, tpv$T1, pch=20, col="blue", cex=3)
    points(tpv$Voc, tpv$T2, pch=1, col="blue", cex=3)
    
    abline(h=1e-7)
    
    lines(DCcapacitance$Voc, RCtime)
    
    graphics.off()
  }
}
