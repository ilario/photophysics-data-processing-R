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
  library(RColorBrewer)

xlim_time = c(-5e-7,0.5e-5)
mycolors=brewer.pal(8,"Dark2")

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
  tpc.sun.inv <- -tpc.sun$V2
  tpc.dark.inv <- -tpc.dark$V2
  
  ce.max <- max(ce$V2)
  tpc.sun.max <- max(tpc.sun.inv)
  tpc.dark.max <- max(tpc.dark.inv)
  tpc.sun.min <- mean(tpc.sun.inv[1:50])
  tpc.dark.min <- mean(tpc.dark.inv[1:50])
  tpv.max <- max(tpv$V2)
  ce.min <- mean(ce$V2[1:50])
  tpv.min <- mean(tpv$V2[1:50])
  #main <- gsub(".txt.table","",gsub("CE_","",ce.file))
  #main <- strsplit(main, "_")[[1]][1]
  main <- basename(dirname(normalizePath(cedir)))
  png(paste("tpc_vs_tpv_vs_ce-", main, ".png", sep=""), width=600, height=600)
  plot(tpv$V1, (tpv$V2-tpv.min)/(tpv.max-tpv.min), xlim=xlim_time, ylim=c(-0.1,1), type="l", ylab="Normalized Voltage", xlab="Time (s)", col="blue", yaxt='n', cex.axis=1.4, cex.lab=1.4)#main=paste(main, "TPC dark vs TPC 1sun vs TPV vs CE"),)
  lines(ce$V1, (ce$V2-ce.min)/(ce.max-ce.min), lwd=2, col="green")
  lines(tpc.sun$V1, (tpc.sun.inv-tpc.sun.min)/(max(tpc.sun.max-tpc.sun.min, tpc.dark.max-tpc.dark.min)), lwd=1, col="red")
  lines(tpc.dark$V1, (tpc.dark.inv-tpc.dark.min)/(max(tpc.sun.max-tpc.sun.min, tpc.dark.max-tpc.dark.min)), lwd=1, col="black", lty=1)
  legend(x="topright",inset=0.05, c("TPC 1 sun","TPC dark", "TPV 1 sun", "CE 1 sun"), lwd=4, col=c("red","black","blue","green"), cex=2)
  graphics.off()
  

if(output_pdf){
  pdf(paste("tpc_vs_tpc-", main, "-norm.pdf", sep=""), width=image_smallpdf_width, height=image_smallpdf_height, pointsize=7)
}else{
  png(paste("tpc_vs_tpc-", main, "-norm.png", sep=""), width=image_width, height=image_height)
}
  plot(NULL, xlim=xlim_time, ylim=c(-0.1,1), ylab="Voltage", xlab="Time (s)", yaxt='n', xaxt='n', cex.axis=1.4, cex.lab=1.7)
  eaxis(side=1, cex.axis=1.4)
  lines(tpc.sun$V1, (tpc.sun.inv-tpc.sun.min)/(max(tpc.sun.max-tpc.sun.min, tpc.dark.max-tpc.dark.min)), lwd=2, col=mycolors[1])
  lines(tpc.dark$V1, (tpc.dark.inv-tpc.dark.min)/(max(tpc.sun.max-tpc.sun.min, tpc.dark.max-tpc.dark.min)), lwd=2, col=mycolors[2])
  legend(x="topright",inset=0.05, c("1 sun","dark"), lwd=3, col=mycolors[c(1,2)], cex=1.5, bty="n")
  graphics.off()


  if(output_pdf){
  pdf(paste("tpc_vs_tpc-", main, ".pdf", sep=""), width=image_smallpdf_width, height=image_smallpdf_height, pointsize=7)
}else{
  png(paste("tpc_vs_tpc-", main, ".png", sep=""), width=image_width, height=image_height)
}
op <- par(mar = c(5,7,1,6) + 0.1) ## default is c(5,4,4,2) + 0.1 

ylim_tpc = range(c(tpc.sun$V2,tpc.dark$V2))
  plot(NULL, xlim=xlim_time, ylim=ylim_tpc, ylab="", xlab="Time (s)", yaxt='n', cex.axis=1.4, cex.lab=1.7, xaxt='n')
title(ylab = "Voltage (V)", cex.lab = 1.7, line = 5)
  eaxis(side=1, cex.axis=1.4, at=c(0,2e-6,4e-6))
  eaxis(side=2, cex.axis=1.4)

  lines(tpc.dark$V1, tpc.dark$V2, lwd=2, col=mycolors[1])
  lines(tpc.sun$V1, tpc.sun$V2, lwd=2, col=mycolors[2])
  legend(x="bottomright",inset=0.05, c("dark","1 sun"), lwd=3, col=mycolors[c(1,2)], cex=1.5, bty="n")

      par(new=TRUE)
      plot(NULL, xaxt="n",yaxt="n",xlab="",ylab="", xlim=xlim_time, ylim=1000*ylim_tpc/50)
      mtext("Current (mA)", cex=1.7, side=4, line=4.5)
  eaxis(side=4, cex.axis=1.4)
  graphics.off()
#reset the plotting margins
par(op)

  png(paste("tpv_vs_ce-", main, ".png", sep=""), width=600, height=600)
  plot(tpv$V1, (tpv$V2-tpv.min)/(tpv.max-tpv.min), xlim=xlim_time, ylim=c(-0.1,1), type="l", ylab="Normalized Voltage", xlab="Time (s)", col="blue", yaxt='n', xaxt='n', cex.axis=1.4, cex.lab=1.4)#, main=paste(main, "TPV vs CE")
  eaxis(side=1, cex.axis=1.2)
  lines(ce$V1, (ce$V2-ce.min)/(ce.max-ce.min), lwd=2, col="green")
  legend(x="topright",inset=0.05, c("TPV 1 sun", "CE 1 sun"), lwd=4, col=c("blue","green"), cex=2)
  graphics.off()



  
}

