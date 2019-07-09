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

#title=gsub("_"," ",tail(unlist(strsplit(name,"-")),1))
filename=gsub(",","",gsub(":","",name))

ylim=lim.TPV.lifetime
xlim=lim.TPV.voltage

output.biexp=list()
output.monoexp=list()
output.robustmonoexp=list()

library(RColorBrewer)
library(sfsmisc)
library(Hmisc)

change.lightness <- function(col, lightness=1){
  if(missing(col))
    stop("Please provide a vector of colours.")
  apply(sapply(col, col2rgb, alpha=TRUE)/255, 2,
        function(x)
          rgb(x[1]*lightness, x[2]*lightness, x[3]*lightness, alpha=x[4]))
}

dirs <- list.dirs(recursive=FALSE)
dirs <- sub("./","",dirs)
legend=sub("_.*","",sub("^0","",dirs))

# try to obtain the color from the file name
mycolors=gsub(".*-col_","",dirs[grepl("-col_", dirs)])
# if the color is not set, use the default one
if(!length(mycolors)){mycolors=brewer.pal(8,"Dark2")}

print("errori derivano da avere molti header in output, bisogna pulirlo")

tryCatch({
  print("biexp")
  
  i <- 0
  if(output_pdf){
    pdf(paste(filename, "-TPVs-biexp.pdf",sep=""), width=image_smallpdf_width, height=image_smallpdf_height, pointsize=7);
  }else{
    png(paste(filename, "-TPVs-biexp.png",sep=""), width=image_width, height=image_height);
  }
  op <- par(mar = c(5,7,1,1) + 0.1) ## default is c(5,4,4,2) + 0.1
  plot(NULL, ylim=ylim, xlim=xlim, cex.axis=1.4, log="y", xlab="", ylab="", yaxt="n")
  title(ylab = "Small perturbation lifetime (s)", cex.lab = 1.7, line = 4)
  title(xlab = "Light bias (V)", cex.lab = 1.7, line = 3)
  
  eaxis(side=2,at=c(1e-10,1e-9,1e-8,1e-7,1e-6,1e-5,1e-4,1e-3,1e-2,0.1,1,10,100,1e3), cex.axis=1.4)
  minor.tick(nx=10)
  
  lapply(dirs, function(x) {print(x);
    subdirs <- list.dirs(path=x, recursive=F)
    subdirs.tpv <- subdirs[grep("tpv", subdirs, ignore.case=T)]
    fulloutput <- read.table(file.path(subdirs.tpv,"output-biexp.txt"), header=TRUE);
    n<-tail(grep("file",fulloutput[,1]),n=1)
    output <- read.table(file.path(subdirs.tpv,"output-biexp.txt"), header=TRUE, skip=ifelse(length(n),n,0)); 
    output.biexp[[paste("Voc",sub("_.*","",sub("^0","",x)),sep="")]] <<- signif(output$Voc,5)
    output.biexp[[paste(sub("_.*","",sub("^0","",x)),"T1",sep="")]] <<- signif(output$T1,5)
    output.biexp[[paste(sub("_.*","",sub("^0","",x)),"T2",sep="")]] <<- signif(output$T2,5)
    points(output$Voc, output$T1, pch=21+(i%%5), col=change.lightness(mycolors[i+1],0.5), bg=mycolors[i+1], cex=1.5);
    points(output$Voc, output$T2, pch=21+(i%%5), col=change.lightness(mycolors[i+1],0.5), cex=1.5);
    i <<- i+1
  })
  legend(x="topright",inset=0.05,legend,pt.cex=2, pt.lwd=1.5, cex=1.5, pch=seq(21,25), pt.bg=mycolors,title=title, bty="n", col=change.lightness(mycolors,0.5))
  graphics.off()
  #reset the plotting margins
  par(op)
  
  maxlength.biexp = max(sapply(output.biexp,length))
  output.biexp = lapply(output.biexp, function(x){length(x)=maxlength.biexp; print(x)})
  output.biexp = as.data.frame(output.biexp,check.names=FALSE)
  write.table(output.biexp, file=paste(filename,"-TPVs-biexp.csv",sep=""), row.names=FALSE, na="", sep=",")
}, error=function(e){cat("error in biexp plotting - ", e$message, "\n")})

tryCatch({
  print("monoexp")
  i <- 0
  if(output_pdf){
    pdf(paste(filename, "-TPVs-monoexp.pdf",sep=""), width=image_smallpdf_width, height=image_smallpdf_height, pointsize=7);
  }else{
    png(paste(filename, "-TPVs-monoexp.png",sep=""), width=image_width, height=image_height);
  }
  op <- par(mar = c(5,7,1,1) + 0.1) ## default is c(5,4,4,2) + 0.1
  plot(NULL, ylim=ylim, xlim=xlim,cex.axis=1.4, log="y", xlab="", ylab="", yaxt="n")
  title(ylab = "Small perturbation lifetime (s)", cex.lab = 1.7, line = 4)
  title(xlab = "Light bias (V)", cex.lab = 1.7, line = 3)
  
  eaxis(side=2,at=c(1e-10,1e-9,1e-8,1e-7,1e-6,1e-5,1e-4,1e-3,1e-2,0.1,1,10,100,1e3), cex.axis=1.4)
  minor.tick(nx=10)
  
  lapply(dirs, function(x) {print(x);
    subdirs <- list.dirs(path=x, recursive=F)
    subdirs.tpv <- subdirs[grep("tpv", subdirs, ignore.case=T)]
    fulloutput <- read.table(file.path(subdirs.tpv,"output-monoexp.txt"), header=TRUE);
    n<-tail(grep("file",fulloutput[,1]),n=1)
    output <- read.table(file.path(subdirs.tpv,"output-monoexp.txt"), header=TRUE, skip=ifelse(length(n),n,0)); 
    output.monoexp[[paste("Voc",sub("_.*","",sub("^0","",x)),sep="")]] <<- signif(output$Voc,5)
    output.monoexp[[sub("_.*","",sub("^0","",x))]] <<- signif(output$T,5)
    points(output$Voc, output$T, lwd=1, pch=21+(i%%5), col=change.lightness(mycolors[i+1],0.5), bg=mycolors[i+1], cex=1.5);
    i <<- i+1
  })
  legend(x="topright",inset=0.05,legend,pt.cex=2, pt.lwd=1.5,cex=1.5, pch=seq(21,25), pt.bg=mycolors, col=change.lightness(mycolors,0.5), title=title, bty="n")
  
  graphics.off()
  #reset the plotting margins
  par(op)
  maxlength.monoexp = max(sapply(output.monoexp,length))
  output.monoexp = lapply(output.monoexp, function(x){length(x)=maxlength.monoexp})
  output.monoexp = as.data.frame(output.monoexp,check.names=FALSE)
  write.table(output.monoexp, file=paste(filename,"-TPVs-monoexp.csv",sep=""), row.names=FALSE, na="", sep=",")
}, error=function(e){cat("error in monoexp plotting - ", e$message, "\n")})

tryCatch({
  print("robust monoexp")
  i <- 0
  if(output_pdf){
    pdf(paste(filename, "-TPVs-robustmonoexp.pdf",sep=""), width=image_smallpdf_width, height=image_smallpdf_height, pointsize=7);
  }else{
    png(paste(filename, "-TPVs-robustmonoexp.png",sep=""), width=image_width, height=image_height);
  }
  op <- par(mar = c(5,7,1,1) + 0.1) ## default is c(5,4,4,2) + 0.1
  plot(NULL, ylim=ylim, xlim=xlim,cex.axis=1.4, log="y", xlab="", ylab="", yaxt="n")
  title(ylab = "Small perturbation lifetime (s)", cex.lab = 1.7, line = 4)
  title(xlab = "Light bias (V)", cex.lab = 1.7, line = 3)
  
  eaxis(side=2,at=c(1e-10,1e-9,1e-8,1e-7,1e-6,1e-5,1e-4,1e-3,1e-2,0.1,1,10,100,1e3), cex.axis=1.4)
  minor.tick(nx=10)
  
  lapply(dirs, function(x) {print(x);
    subdirs <- list.dirs(path=x, recursive=F)
    subdirs.tpv <- subdirs[grep("tpv", subdirs, ignore.case=T)]
    fulloutput <- read.table(file.path(subdirs.tpv,"output-robustmonoexp.txt"), header=TRUE);
    n<-tail(grep("file",fulloutput[,1]),n=1)
    output <- read.table(file.path(subdirs.tpv,"output-robustmonoexp.txt"), header=TRUE, skip=ifelse(length(n),n,0)); 
    output.robustmonoexp[[paste("Voc",sub("_.*","",sub("^0","",x)),sep="")]] <<- signif(output$Voc,5)
    output.robustmonoexp[[sub("_.*","",sub("^0","",x))]] <<- signif(output$T,5)
    points(output$Voc, output$T, pch=21+(i%%5), col=change.lightness(mycolors[i+1],0.5), bg=mycolors[i+1], cex=1.5);
    i <<- i+1
  })
  legend(x="topright",inset=0.05,legend,pt.cex=2, pt.lwd=1.5,cex=1.5, pch=seq(21,25), pt.bg=mycolors, col=change.lightness(mycolors,0.5), title=title, bty="n")
  
  graphics.off()
  #reset the plotting margins
  par(op)
  maxlength.robustmonoexp = max(sapply(output.robustmonoexp,length))
  output.robustmonoexp = lapply(output.robustmonoexp, function(x){length(x)=maxlength.robustmonoexp})
  output.robustmonoexp = as.data.frame(output.robustmonoexp,check.names=FALSE)
  write.table(output.robustmonoexp, file=paste(filename,"-TPVs-robustmonoexp.csv",sep=""), row.names=FALSE, na="", sep=",")
}, error=function(e){cat("error in robust monoexp plotting - ", e$message, "\n")})




