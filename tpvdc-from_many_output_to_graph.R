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

ylim=lim.TPVDC.lifetime
xlim=lim.TPVDC.charge
ylimnogeom=lim.TPVDC.nogeom.lifetime
xlimnogeom=lim.TPVDC.nogeom.charge

output=list()
output.nogeom=list()

library(robustbase)
library(RColorBrewer)
library(minpack.lm)
library(sfsmisc)
library(Hmisc)

dirs <- list.dirs(recursive=FALSE)
dirs <- sub("./","",dirs)
legend=sub("_.*","",sub("^0","",dirs))

# try to obtain the color from the file name
mycolors=gsub(".*-col_","",dirs[grepl("-col_", dirs)])
# if the color is not set, use the default one
if(!length(mycolors)){mycolors=brewer.pal(8,"Dark2")}

i <- 0
png(paste(filename,"-TPVDCs.png",sep=""), width=image_width, height=image_height)
par(mar=c(5.1,7,2,2.1))
plot(1,xlim=xlim,ylim=ylim,cex.main=1.5,xlab=bquote("Charge Density (C/cm"^"2"*")"), ylab="Life-time (s)",cex.lab=1.5,cex.axis=1.2,log="y", yaxt="n", xaxt="n");#, main=paste(name,"TPV decay vs Charge from DC")
eaxis(side=2,at=c(1e-10,1e-9,1e-8,1e-7,1e-6,1e-5,1e-4,1e-3,1e-2,0.1,1,10,100,1e3), cex.axis=1.2)
eaxis(side=1, cex.axis=1.2)
minor.tick(nx=10)

lapply(dirs, function(x) {print(x);
  subdirs <- list.dirs(path=x, recursive=F)
  subdirs.tpv <- subdirs[grep("tpv", subdirs, ignore.case=T)]
  a <- read.table(file.path(x,"outputDCcharge.txt"),header=T,stringsAsFactors=F)
  lo <- loess(a$ChargeDensityDC~a$Voc,span=0.9)
  #just for the slope, no intercept
  linearfit <- lm(ChargeDensityDC ~ 0 + Voc, data=a[1:round(length(a$Voc)/2),])
  startlist=list(B=coef(linearfit),C=1e-10,D=8)
  tryCatch({
    expfit <- nlsLM(ChargeDensityDC~ B*Voc+C*(exp(D*Voc)-1), start=startlist, data=a)
    tryCatch({
      expfit <- nlrob(ChargeDensityDC~ B*Voc+C*(exp(D*Voc)-1), start=list(B=coef(expfit)[1], C=coef(expfit)[2], D=coef(expfit)[3]), data=a)
    }, error=function(e) {print("FAILED ROBUST FIT")});
  }, error=function(e) {print("FAILED non-robust FIT")});
  
  # if things go bad, the linear fit should suffice
  if(!exists("expfit")){
    print("USING LINEAR FIT FOR CHARGE EXTRACTION DATA")
    expfit <- linearfit
  }
  
  filex <- file.path(subdirs.tpv, "output-monoexp.txt")
  
  fulloutput <- read.table(filex, header=TRUE);
  n<-tail(grep("file",fulloutput[,1]),n=1)
  tpv <- read.table(filex, header=TRUE, skip=ifelse(length(n),n,0));
  #importante che la variabile in new abbia lo stesso nome di quella fittata
  new <- data.frame(Voc = tpv$Voc)
  charge <- (predict(lo,tpv$Voc)+predict(expfit,new))/2
  new2 <- data.frame(Voc = tpv$Voc[is.na(charge)])
  charge[is.na(charge)] <- predict(expfit,new2)
  output[[paste("Charge",sub("nm","",sub("_.*","",sub("^0","",x))),sep="")]] <<- signif(charge,5)
  output[[sub("_.*","",sub("^0","",x))]] <<- signif(tpv$T,5)
  points(charge, tpv$T, lwd=1, bg=mycolors[i+1], cex=2, pch=21+(i%%5));
  i <<- i+1
})
legend(x="topright",inset=0.05,legend,pch=seq(21,25), pt.bg=mycolors, lwd=4, pt.lwd=2, pt.cex=2, col=mycolors,cex=1.5, title=#paste("TPV vs DC\n","with geom. cap.\n",
         title,bg="gray90"#), bty="n"
)
graphics.off()

maxlength = max(sapply(output,length))
output = lapply(output, function(x){length(x)=maxlength; print(x)})
output = as.data.frame(output,check.names=FALSE)
write.table(output, file=paste(filename,"-TPVDCs.csv",sep=""), row.names=FALSE, na="", sep=",")


i <- 0
png(paste(filename,"-TPVDCs-nogeom.png",sep=""), width=image_width, height=image_height)
par(mar=c(5.1,7,2,2.1))
plot(1,xlim=xlimnogeom,ylim=ylimnogeom,cex.main=1.5,xlab=bquote("Charge Density (C/cm"^"2"*")"), ylab="Life-time (s)",cex.lab=1.5,cex.axis=1.2,log="y", yaxt="n", xaxt="n");#, main=paste(name,"TPV decay vs Charge from DC")
eaxis(side=2,at=c(1e-10,1e-9,1e-8,1e-7,1e-6,1e-5,1e-4,1e-3,1e-2,0.1,1,10,100,1e3), cex.axis=1.2)
eaxis(side=1, cex.axis=1.2)
minor.tick(nx=10)

lapply(dirs, function(x) {print(x);
  subdirs <- list.dirs(path=x, recursive=F)
  subdirs.tpv <- subdirs[grep("tpv", subdirs, ignore.case=T)]
  a <- read.table(file.path(x,"outputDCcharge-nogeom.txt"),header=T,stringsAsFactors=F)
  lo <- loess(a$ChargeDensityDC~a$Voc,span=0.9)
  startlist=list(C=1e-10,D=8)
  tryCatch({
    expfit <- nlsLM(ChargeDensityDC~ C*(exp(D*Voc)-1), start=startlist, data=a)
    tryCatch({
      expfit <- nlrob(ChargeDensityDC~ C*(exp(D*Voc)-1), start=list(C=coef(expfit)[1],D=coef(expfit)[2]), data=a)
    }, error=function(e) {print("FAILED nogeom FIT ROBUST")});
  }, error=function(e) {print("FAILED nogeom FIT non-robust")});
  
  filex <- file.path(subdirs.tpv, "output-monoexp.txt")
  fulloutput <- read.table(filex, header=TRUE);
  n<-tail(grep("file",fulloutput[,1]),n=1)
  tpv <- read.table(filex, header=TRUE, skip=ifelse(length(n),n,0));
  #importante che la variabile in new abbia lo stesso nome di quella fittata
  new <- data.frame(Voc = tpv$Voc)
  charge <- (predict(lo,tpv$Voc)+predict(expfit,new))/2
  new2 <- data.frame(Voc = tpv$Voc[is.na(charge)])
  charge[is.na(charge)] <- predict(expfit,new2)
  output.nogeom[[paste("Charge",sub("nm","",sub("_.*","",sub("^0","",x))),sep="")]] <<- signif(charge,5)
  output.nogeom[[sub("_.*","",sub("^0","",x))]] <<- signif(tpv$T,5)
  points(charge, tpv$T, lwd=1, bg=mycolors[i+1], cex=2, pch=21+(i%%5));
  i <<- i+1
})
legend(x="topright",inset=0.05,legend,pch=seq(21,25), pt.bg=mycolors, lwd=4, pt.lwd=2, pt.cex=2, col=mycolors,cex=1.5, title=#paste("TPV vs DC\n","no geom. cap.\n",
         title, bg="gray90"#, bty="n"
)
graphics.off()

output.nogeom = lapply(output.nogeom, function(x){length(x)=maxlength; print(x)})
output.nogeom = as.data.frame(output.nogeom,check.names=FALSE)
write.table(output.nogeom, file=paste(filename,"-TPVDCs-nogeom.csv",sep=""), row.names=FALSE, na="", sep=",")
