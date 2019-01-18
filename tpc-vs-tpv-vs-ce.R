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

tpcVsTpvVsCe <- function(cedir="ce", tpvdir="tpv", tpcdir="tpc")
{
  require(sfsmisc)
  
  print("TPV vs CE vs TPC dark vs TPC sun: PLOTTING")
  ce.files <- list.files(path=cedir, pattern="^CE.*\\.txt.table$");
  tpv.files  <- list.files(path=tpvdir, pattern="^TPV.*\\.txt.table$");
  tpc.files <- list.files(path=tpcdir, pattern="^TPC.*\\.txt.table$");
  tpc.sun.files <- tpc.files[grep("sun", tpc.files, ignore.case=T)]
  tpc.dark.files <- tpc.files[grep("dark", tpc.files, ignore.case=T)]
  
  ce.file <- tail(ce.files, n=1)
  tpv.file <- tail(tpv.files, n=1)
  tpc.sun.file <- tail(tpc.sun.files, n=1)
  tpc.dark.file <- tail(tpc.dark.files, n=1)
  
  print(paste("CE file:", ce.file))
  print(paste("TPV file:", tpv.file))
  print(paste("TPC sun file:", tpc.sun.file))
  print(paste("TPC dark file:", tpc.dark.file))
  
  ce <- read.table(file.path(cedir, ce.file),header=F)
  tpv <- read.table(file.path(tpvdir, tpv.file),header=F)
  tpc.sun <- read.table(file.path(tpcdir, tpc.sun.file),header=F)
  tpc.dark <- read.table(file.path(tpcdir, tpc.dark.file),header=F)
  tpc.sun$V2 <- -tpc.sun$V2
  tpc.dark$V2 <- -tpc.dark$V2
  
  ce.max <- max(ce$V2)
  tpc.sun.max <- max(tpc.sun$V2)
  tpc.dark.max <- max(tpc.dark$V2)
  tpc.sun.min <- mean(tpc.sun$V2[1:50])
  tpc.dark.min <- mean(tpc.dark$V2[1:50])
  tpv.max <- max(tpv$V2)
  ce.min <- mean(ce$V2[1:50])
  tpv.min <- mean(tpv$V2[1:50])
  #main <- gsub(".txt.table","",gsub("CE_","",ce.file))
  #main <- strsplit(main, "_")[[1]][1]
  main <- basename(dirname(normalizePath(cedir)))
  png(paste("tpc_vs_tpv_vs_ce-", main, ".png", sep=""), width=600, height=600)
  plot(tpv$V1, (tpv$V2-tpv.min)/(tpv.max-tpv.min), xlim=c(-5e-7,0.5e-5), ylim=c(-0.1,1), type="l", ylab="Normalized Voltage", xlab="Time (s)", col="blue", yaxt='n', cex.axis=1.4, cex.lab=1.4)#main=paste(main, "TPC dark vs TPC 1sun vs TPV vs CE"),)
  lines(ce$V1, (ce$V2-ce.min)/(ce.max-ce.min), lwd=2, col="green")
  lines(tpc.sun$V1, (tpc.sun$V2-tpc.sun.min)/(max(tpc.sun.max-tpc.sun.min, tpc.dark.max-tpc.dark.min)), lwd=1, col="red")
  lines(tpc.dark$V1, (tpc.dark$V2-tpc.dark.min)/(max(tpc.sun.max-tpc.sun.min, tpc.dark.max-tpc.dark.min)), lwd=1, col="black", lty=1)
  legend(x="topright",inset=0.05, c("PI-TPC 1sun","PI-TPC dark", "PI-TPV 1sun", "PICE 1sun"), lwd=4, col=c("red","black","blue","green"), cex=2)
  graphics.off()
  
  png(paste("tpc_vs_tpc-", main, ".png", sep=""), width=600, height=600)
  plot(NULL, xlim=c(-5e-7,0.5e-5), ylim=c(-0.1,1), type="l", ylab="Voltage", xlab="Time (s)", yaxt='n', xaxt='n', cex.axis=1.4, cex.lab=1.4)# main=paste(main, "TPC dark vs TPC 1sun"),
  eaxis(side=1, cex.axis=1.2)
  lines(tpc.sun$V1, (tpc.sun$V2-tpc.sun.min)/(max(tpc.sun.max-tpc.sun.min, tpc.dark.max-tpc.dark.min)), lwd=1, col="red")
  lines(tpc.dark$V1, (tpc.dark$V2-tpc.dark.min)/(max(tpc.sun.max-tpc.sun.min, tpc.dark.max-tpc.dark.min)), lwd=1, col="black", lty=1)
  legend(x="topright",inset=0.05, c("PI-TPC 1sun","PI-TPC dark"), lwd=4, col=c("red","black"), cex=2)
  graphics.off()
  
  png(paste("tpv_vs_ce-", main, ".png", sep=""), width=600, height=600)
  plot(tpv$V1, (tpv$V2-tpv.min)/(tpv.max-tpv.min), xlim=c(-5e-7,0.5e-5), ylim=c(-0.1,1), type="l", ylab="Normalized Voltage", xlab="Time (s)", col="blue", yaxt='n', xaxt='n', cex.axis=1.4, cex.lab=1.4)#, main=paste(main, "TPV vs CE")
  eaxis(side=1, cex.axis=1.2)
  lines(ce$V1, (ce$V2-ce.min)/(ce.max-ce.min), lwd=2, col="green")
  legend(x="topright",inset=0.05, c("PI-TPV 1sun", "PICE 1sun"), lwd=4, col=c("blue","green"), cex=2)
  graphics.off()
  
}

