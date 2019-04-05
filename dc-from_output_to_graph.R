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

dcFromOutputToGraph <- function(tpvdir="tpv", tpcdir="tpc")
{
  print("DC: PLOTTING")
  require(minpack.lm)
  require(robustbase)
  require(sfsmisc)
  require(Hmisc)
  library(RColorBrewer)
  
  mycolors=brewer.pal(8,"Dark2")
  
  write.table(t(c("Voc","capacitance")), file="outputDCcapacitance.txt", append=FALSE, col.names=F, row.names=F);
  write.table(t(c("Voc","ChargeDensityDC")), file="outputDCcharge.txt", append=FALSE, col.names=F, row.names=F);
  write.table(t(c("Voc","ChargeDensityDC")), file="outputDCcharge-nogeom.txt", append=FALSE, col.names=F, row.names=F);
  
  a <- read.table(file.path(tpcdir, "outputChargeDensityTPC.txt"), header=T)
  
  # in case TPC in dark and in sun are different, the choice of what to use is arbitrary, I would use the third quartile of all the TPC measurements
  charge <- quantile(a$ChargeDensityTPC, 0.75)
  
#    b <- read.table(file.path(tpvdir, "outputDeltaVbiexp.txt"), header=T)
    b <- read.table(file.path(tpvdir, "outputDeltaVfirstPoints.txt"), header=T)
#    b <- read.table(file.path(tpvdir, "outputDeltaVloess.txt"), header=T)
#    b <- read.table(file.path(tpvdir, "outputDeltaVmonoexp.txt"), header=T)
#    b <- read.table(file.path(tpvdir, "outputDeltaV.txt"), header=T)
 
  write.table(b, file=file.path(tpvdir, "outputDeltaVprocessedForDC.txt"), append=FALSE, row.names=FALSE)
  
  len <- length(b$deltaV)
  
  chargeArray = rep(charge, len)
  
  getExpFit <- function(){
    capacitance <<- chargeArray/b$deltaV
    directory <<- tail(strsplit(getwd(), "/")[[1]], n=1)
    
    outputDCcapacitance <<- data.frame(b$Voc, capacitance);
    names(outputDCcapacitance) <<- c("Voc","capacitance")
    
    tryCatch({
      expfit <<- nls(capacitance ~ exp(B) + exp(C)*exp(D)*exp(exp(D)*Voc), start=list(B=log(max(1e-8,min(outputDCcapacitance$capacitance))),C=log(1e-10),D=2), data=outputDCcapacitance)
    }, error=function(e) print("Failed restricted to positive gamma fit - first"))
    tryCatch({
      expfit <<- nlsLM(capacitance ~ exp(B) + exp(C)*exp(D)*exp(exp(D)*Voc), start=list(B=log(max(1e-8,min(outputDCcapacitance$capacitance))),C=log(1e-10),D=2), data=outputDCcapacitance)
    }, error=function(e) print("Failed restricted to positive gamma fit - second"))
    tryCatch({
      expfit <<- nlrob(capacitance ~ exp(B) + exp(C)*exp(D)*exp(exp(D)*Voc), start=list(B=coef(expfit)[[1]],C=coef(expfit)[[2]],D=coef(expfit)[[3]]), data=outputDCcapacitance)
    }, error=function(e) print("Failed restricted to positive gamma robust fit"))
  }
  
  
  getExpFit()
  
  write.table(outputDCcapacitance, file="outputDCcapacitance.txt", append=TRUE, col.names=F, row.names=F, quote=F);
  
  if(output_pdf){
    pdf(paste("DC-capacitance-", directory, ".pdf", sep=""), width=image_smallpdf_width, height=image_smallpdf_height, pointsize=7);
  }else{
    png(paste("DC-capacitance-", directory, ".png", sep=""), width=image_width, height=image_height)
  }
  par(mar=c(5,9,1,1))
  plot(b$Voc, capacitance, xlim=c(0, max(b$Voc)), ylim=c(0, max(capacitance)), ylab="", xlab="", xaxt="n", yaxt="n", cex.lab=1.7, panel.first=c(abline(h=0, col="gray80"), abline(v=0, col="gray80")))
  eaxis(side=2, cex.axis=1.4)
  eaxis(side=1, cex.axis=1.4)
  minor.tick(ny=10, nx=10)
  title(xlab="Light bias (V)", cex.lab=1.7, line=3.5)
  title(ylab=bquote("Specific Capacitance (F/cm"^"2"*")"), cex.lab=1.7, line=7)
  lines(outputDCcapacitance$Voc, predict(expfit), lwd=2, col=mycolors[1])
  graphics.off()
  
  write.table(t(c("B","Ch0","gamma","GeomCh","ChemCh")), file="outputDC-fit.txt", append=FALSE, col.names=F, row.names=F);
  eB=exp(coef(expfit)[[1]])
  eCh0=exp(coef(expfit)[[2]])
  egamma=exp(coef(expfit)[[3]])
  Voc=max(outputDCcapacitance$Voc)
  GeomCh=eB*Voc
  ChemCh=eCh0*(exp(egamma*Voc)-1)
  output <- t(c(eB, eCh0, egamma, GeomCh, ChemCh))
  write.table(output, file="outputDC-fit.txt", append=TRUE, col.names=F, row.names=F)
  
  c<- data.frame(b$Voc,capacitance)
  d <- c[with(c, order(b.Voc)), ]
  e <- d[1:(nrow(d)/2),]
  
  f <- data.frame(d$b.Voc, d$capacitance)
  
  names(f) <- c("Voc","capacitance")
  g <- f
  g$capacitance[g$capacitance < 0] <- 0
  
  z <- approxfun(g$Voc, g$capacitance, method="linear", 0, 0)
  integral=Vectorize(function(X)integrate(z,0,X)$value)
  
  outputDCcharge <- data.frame(f$Voc, integral(f$Voc));
  write.table(outputDCcharge, file="outputDCcharge.txt", append=TRUE, col.names=F, row.names=F, quote=F);
  
  
  if(output_pdf){
    pdf(paste("DC-charge-", directory, ".pdf", sep=""), width=image_smallpdf_width, height=image_smallpdf_height, pointsize=7);
  }else{
    png(paste("DC-charge-", directory, ".png", sep=""), width=image_width, height=image_height)
  }
  par(mar=c(5,9,1,1))
  plot(integral, range(f$Voc)[1], range(f$Voc)[2], ylab="", xlab="", xaxt="n", yaxt="n", cex.lab=1.7, panel.first=c(abline(h=0, col="gray80"), abline(v=0, col="gray80")))#, log="y")
  eaxis(side=2, cex.axis=1.4)
  eaxis(side=1, cex.axis=1.4)
  minor.tick(ny=10, nx=10)
  title(xlab="Light bias (V)", cex.lab=1.7, line=3.5)
  title(ylab=bquote("Charge Density (C/cm"^"2"*")"), cex.lab=1.7, line=7)
  graphics.off()
  
  dataframe <- data.frame(Voc=g$Voc,capacitance=g$capacitance)
  
  geometrical <- quantile(dataframe$capacitance, 0.05)
  g$capacitance <- g$capacitance - geometrical
  
  z <- approxfun(g$Voc, g$capacitance, method="linear", 0, 0)
  integral=Vectorize(function(X)integrate(z,0,X)$value)
  
  outputDCcharge <- data.frame(f$Voc, integral(f$Voc));
  write.table(outputDCcharge, file="outputDCcharge-nogeom.txt", append=TRUE, col.names=F, row.names=F, quote=F);
  
  
  if(output_pdf){
    pdf(paste("DC-nogeom-charge-", directory, ".pdf", sep=""), width=image_smallpdf_width, height=image_smallpdf_height, pointsize=7);
  }else{
    png(paste("DC-nogeom-charge-", directory, ".png", sep=""), width=image_width, height=image_height)
  }
  par(mar=c(5,9,1,1))
  plot(integral, range(f$Voc)[1], range(f$Voc)[2], ylab="", xlab="", xaxt="n", yaxt="n", cex.lab=1.7, panel.first=c(abline(h=0, col="gray80"), abline(v=0, col="gray80")))#, log="y")
  eaxis(side=2, cex.axis=1.4)
  eaxis(side=1, cex.axis=1.4)
  minor.tick(ny=10, nx=10)
  title(xlab="Light bias (V)", cex.lab=1.7, line=3.5)
  title(ylab=bquote("Charge Density (C/cm"^"2"*")"), cex.lab=1.7, line=7)
  graphics.off()
}
