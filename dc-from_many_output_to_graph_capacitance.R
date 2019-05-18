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

#title=gsub("_"," ",tail(unlist(strsplit(name,"-")),1))
filename=gsub(",","",gsub(":","",name))

library(RColorBrewer)
library(sfsmisc)
library(Hmisc)
library(robustbase)
library(minpack.lm)

ylim=lim.DCcapacitance.capacitance
xlim=lim.DCcapacitance.voltage
ylimnogeom=lim.DCcapacitance.nogeom.capacitance
xlimnogeom=lim.DCcapacitance.nogeom.voltage

add.alpha <- function(col, alpha=1){
  if(missing(col))
    stop("Please provide a vector of colours.")
  apply(sapply(col, col2rgb)/255, 2,
        function(x)
          rgb(x[1], x[2], x[3], alpha=alpha))
}
change.lightness <- function(col, lightness=1){
  if(missing(col))
    stop("Please provide a vector of colours.")
  apply(sapply(col, col2rgb, alpha=TRUE)/255, 2,
        function(x)
          rgb(x[1]*lightness, x[2]*lightness, x[3]*lightness, alpha=x[4]))
}

output=list()
output.nogeom=list()

i <- 0
dirs <- list.dirs(recursive=FALSE)
dirs <- sub("./","",dirs)
#remove everything after the last underscore, then convert the remaining underscores to spaces
legendlist=sub("^0","",gsub("_"," ",sub("_((?!_).)+$","",dirs, perl=TRUE)))

# try to obtain the color from the file name
mycolors=gsub(".*-col_","",dirs[grepl("-col_", dirs)])
# if the color is not set, use the default one
if(!length(mycolors)){mycolors=brewer.pal(8,"Dark2")}

if(output_pdf){
  pdf(paste(filename,"-DCs-capacitance.pdf",sep=""), width=image_smallpdf_width, height=image_smallpdf_height, pointsize=7)
}else{
  png(paste(filename,"-DCs-capacitance.png",sep=""), width=image_width, height=image_height)
}
op <- par(mar = c(5,8,1,1) + 0.1) ## default is c(5,4,4,2) + 0.1 
plot(NULL,xlim=xlim,ylim=ylim,cex.lab=1.7,xlab="Light bias (V)",ylab="", yaxt="n",xaxt="n", panel.first=c(abline(h=0, col="gray80"), abline(v=0, col="gray80")))
#eaxis(side=2,at=c(1e-12,1e-11,1e-10,1e-9,1e-8,1e-7,1e-6,1e-5,1e-4,1e-3,1e-2,0.1,1,10,100,1e3), cex.axis=1.2)
eaxis(side=2, cex.axis=1.4)
eaxis(side=1, cex.axis=1.4)
minor.tick(nx=10, ny=10)
title(ylab=bquote("Specific capacitance (F/cm"^"2"*")"), mgp=c(6,1,0), cex.lab=1.7)

getExpFit <- function(Voc, capacitance){
  outputDCcapacitance <- data.frame(Voc, capacitance);
  names(outputDCcapacitance) <- c("Voc","capacitance")
  tryCatch({
    expfit <- nls(capacitance ~ exp(B) + exp(C)*exp(D)*exp(exp(D)*Voc), start=list(B=log(max(1e-8,min(outputDCcapacitance$capacitance))),C=log(1e-10),D=2), data=outputDCcapacitance)
  }, error=function(e) print("Failed restricted to positive gamma fit - first"))
  tryCatch({
    expfit <- nlsLM(capacitance ~ exp(B) + exp(C)*exp(D)*exp(exp(D)*Voc), start=list(B=log(max(1e-8,min(outputDCcapacitance$capacitance))),C=log(1e-10),D=2), data=outputDCcapacitance)
  }, error=function(e) print("Failed restricted to positive gamma fit - second"))
  tryCatch({
    expfit <- nlrob(capacitance ~ exp(B) + exp(C)*exp(D)*exp(exp(D)*Voc), start=list(B=coef(expfit)[[1]],C=coef(expfit)[[2]],D=coef(expfit)[[3]]), data=outputDCcapacitance)
  }, error=function(e) print("Failed restricted to positive gamma robust fit"))
  return(expfit)
}

# preallocate geometric capacitance array
geometric = numeric(length(dirs))

lapply(dirs, function(x) {print(x);
  subdirs <- list.dirs(path=x, recursive=F)
  subdirs.tpc <- subdirs[grep("tpc", subdirs, ignore.case=T)]
  subdirs.tpv <- subdirs[grep("tpv", subdirs, ignore.case=T)]
  a <- read.table(file.path(subdirs.tpc,"outputChargeDensityTPC.txt"),header=T)
  # in case TPC in dark and in sun are different, the choice of what to use is arbitrary, I would use the first quartile of all the TPC measurements
  charge <- quantile(a$ChargeDensityTPC, 0.25)
  b <- read.table(file.path(subdirs.tpv, "outputDeltaVprocessedForDC.txt"), header=T)
  
  output[[paste("Voc",sub("nm","",sub("_.*","",sub("^0","",x))),sep="")]] <<- signif(b$Voc,5)
  capacitance <- charge/b$deltaV
  
  expfit = getExpFit(Voc=b$Voc, capacitance=capacitance)
  
#  geometric_quantile <- quantile(capacitance,0.05)
  geometric[i+1] <<- exp(coef(expfit)[[1]])
  
  print(paste("Geometric capacitance for", x, "is", geometric[i+1], "F/cm2"))
  output[[sub("_.*","",sub("^0","",x))]] <<- signif(capacitance,5)
#  abline(h=geometric_quantile, col=add.alpha(change.lightness(mycolors[i+1],0.8),0.6), lwd=2)

  lines(seq(0,max(b$Voc),0.01), predict(expfit, newdata=data.frame(Voc=seq(0,max(b$Voc),0.01))), lwd=2, col=add.alpha(change.lightness(mycolors[i+1],0.8),0.8))
  
  points(b$Voc, capacitance, col=add.alpha(change.lightness(mycolors[i+1],0.5),0.9), bg=add.alpha(mycolors[i+1],0.5), pch=21+(i%%5), cex=1.5)
  i <<- i+1
})
#abline(h=0)
legend(x="topleft",inset=0.05,legendlist,pch=seq(21,25), pt.bg=mycolors,pt.cex=2, cex=1.5, pt.lwd=1.5,col=change.lightness(mycolors,0.5), lwd=3, title=title, bty="n")
graphics.off()

maxlength = max(sapply(output,length))
output = lapply(output, function(x){length(x)=maxlength; print(x)})
output = as.data.frame(output,check.names=FALSE)
write.table(output, file=paste(filename,"-DCs-capacitance.csv",sep=""), row.names=FALSE, na="", sep=",")

i <- 0
if(output_pdf){
  pdf(paste(filename,"-DCs-nogeom-capacitance.pdf",sep=""), width=image_smallpdf_width, height=image_smallpdf_height, pointsize=7)
}else{
  png(paste(filename,"-DCs-nogeom-capacitance.png",sep=""), width=image_width, height=image_height)
}
op <- par(mar = c(5,8,1,1) + 0.1) ## default is c(5,4,4,2) + 0.1 
plot(NULL,xlim=xlimnogeom,ylim=ylimnogeom,cex.lab=1.7,xlab="Light bias (V)",ylab="", yaxt="n",xaxt="n", panel.first=c(abline(h=0, col="gray80"), abline(v=0, col="gray80")))
#eaxis(side=2,at=c(1e-12,1e-11,1e-10,1e-9,1e-8,1e-7,1e-6,1e-5,1e-4,1e-3,1e-2,0.1,1,10,100,1e3), cex.axis=1.2)
eaxis(side=2, cex.axis=1.4)
eaxis(side=1, cex.axis=1.4)
minor.tick(nx=10, ny=10)
title(ylab=bquote("Specific capacitance (F/cm"^"2"*")"), mgp=c(6,1,0), cex.lab=1.7)

lapply(dirs, function(x) {print(x);
  subdirs <- list.dirs(path=x, recursive=F)
  subdirs.tpc <- subdirs[grep("tpc", subdirs, ignore.case=T)]
  subdirs.tpv <- subdirs[grep("tpv", subdirs, ignore.case=T)]
  a <- read.table(file.path(subdirs.tpc,"outputChargeDensityTPC.txt"),header=T)
  # in case TPC in dark and in sun are different, the choice of what to use is arbitrary, I would use the first quartile of all the TPC measurements
  charge <- quantile(a$ChargeDensityTPC, 0.25)
  b <- read.table(file.path(subdirs.tpv, "outputDeltaVprocessedForDC.txt"), header=T)
  output.nogeom[[paste("Voc",sub("nm","",sub("_.*","",sub("^0","",x))),sep="")]] <<- signif(b$Voc,5)
  capacitance <- charge/b$deltaV
#  dataframe <- data.frame(Voc=b$Voc,capacitance=capacitance)
  
  capacitance <- capacitance - geometric[i+1]
  
  output.nogeom[[sub("_.*","",sub("^0","",x))]] <<- signif(capacitance,5)
  points(b$Voc, capacitance, col=add.alpha(change.lightness(mycolors[i+1],0.5),0.6), bg=add.alpha(mycolors[i+1],0.5), pch=21+(i%%5), cex=1.5)
  i <<- i+1
})
#abline(h=0)
legend(x="topleft",inset=0.05,legendlist,pch=seq(21,25), pt.bg=mycolors,pt.cex=2, cex=1.5, pt.lwd=1.5,col=change.lightness(mycolors,0.5), title=title, bty="n")
graphics.off()

output.nogeom = lapply(output.nogeom, function(x){length(x)=maxlength; print(x)})
output.nogeom = as.data.frame(output.nogeom,check.names=FALSE)
write.table(output.nogeom, file=paste(filename,"-DCs-nogeom-capacitance.csv",sep=""), row.names=FALSE, na="", sep=",")
#reset the plotting margins
par(op)
