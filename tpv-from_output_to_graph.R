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

tpvFromOutputToGraph <- function(tpvdir="tpv")
{
  print("TPV: PLOTTING")
  library(Hmisc)
  library(sfsmisc)
  require(minpack.lm)
  require(robustbase)
  require(RColorBrewer)

  directory <- tail(strsplit(getwd(), "/")[[1]], n=2)
  name<- directory[1]
  mycolors=brewer.pal(8,"Dark2")

change.lightness <- function(col, lightness=1){
  if(missing(col))
    stop("Please provide a vector of colours.")
  apply(sapply(col, col2rgb, alpha=TRUE)/255, 2,
        function(x)
          rgb(x[1]*lightness, x[2]*lightness, x[3]*lightness, alpha=x[4]))
}

  print("biexp")
  fulloutput <- read.table(file.path(tpvdir,"output-biexp.txt"), header=TRUE);#,stringsAsFactors=F);
  n<-tail(grep("file",fulloutput[,1]),n=1)
  tpv <- read.table(file.path(tpvdir,"output-biexp.txt"), header=TRUE, skip=ifelse(length(n),n,0)); 
  if(output_pdf){
    pdf(file.path(tpvdir, paste(name, "-tpv-biexp.pdf",sep="")), width=image_smallpdf_width, height=image_smallpdf_height, pointsize=7);
  }else{
    png(file.path(tpvdir, paste(name, "-tpv-biexp.png",sep="")), width=image_width, height=image_height);
  }
  plot(1, ylim=c(min(tpv$T1,tpv$T2),max(tpv$T2,tpv$T1)), xlim=c(min(tpv$Voc),max(tpv$Voc)), log="y", xlab="Voc (V)", ylab="Life-time (s)", main=paste(name, "TPV biexp"), cex.axis=1.5, cex.lab=1.5);
  #errbar(tpv$Voc, tpv$T1, tpv$T1+tpv$T1.error, tpv$T1-tpv$T1.error,log="y",add=TRUE, pch="")
  #errbar(tpv$Voc, tpv$T2, tpv$T2+tpv$T2.error, tpv$T2-tpv$T2.error,log="y",add=TRUE, pch="")
  points(tpv$Voc, tpv$T1, pch=20, col="red")#, cex=4*tpv$A1/(tpv$A1+tpv$A2));
  points(tpv$Voc, tpv$T2, pch=20, col="black")#, cex=4*tpv$A2/(tpv$A1+tpv$A2));
  graphics.off()
  
  
  print("monoexp")
  
  fulloutput <- read.table(file.path(tpvdir,"output-monoexp.txt"), header=TRUE);
  n<-tail(grep("file",fulloutput[,1]),n=1)
  tpv <- read.table(file.path(tpvdir,"output-monoexp.txt"), header=TRUE, skip=ifelse(length(n),n,0));
  tpv <- tpv[with(tpv, order(tpv$Voc)),]
  lo<-loess(tpv$T~tpv$Voc,span=0.5)
  
  tpvEnding <- tpv[round(length(tpv$V)*0.5):length(tpv$V),]
  fit <- nlsLM(T~B*exp(C*Voc), start=list(B=4,C=-15), data=tpvEnding)
  tryCatch({
    fit <- nlrob(T~B*exp(C*Voc), start=list(B=coef(fit)["B"],C=coef(fit)["C"]), data=tpvEnding)
  }, error=function(e) print("Failed robust fit"))
  if(output_pdf){
    pdf(file.path(tpvdir, paste(name, "-tpv-monoexp.pdf",sep="")), width=image_smallpdf_width, height=image_smallpdf_height, pointsize=7);
  }else{
    png(file.path(tpvdir, paste(name, "-tpv-monoexp.png",sep="")), width=image_width, height=image_height);
  }
  par(mar=c(5.1,5,4.1,2.1))
  plot(tpv$Voc, tpv$T, col="black", log="y", xlab=bquote("V"["oc"]~"(V)"), ylab="Life-time (s)", cex.axis=1.4, cex.lab=1.4)#, main=paste(name, "TPV monoexp"));
  #errbar(tpv$Voc, tpv$T, tpv$T+tpv$T.error, tpv$T-tpv$T.error,log="y",add=TRUE)
  lines(tpv$Voc, predict(lo), lwd=1,col="red")
  lines(tpvEnding$V,predict(fit),col="green",lwd=2)
  graphics.off()
  
  write.table(t(c("Tau0","beta")), file=file.path(tpvdir,"output-monoexp-fit.txt"), append=FALSE, col.names=F, row.names=F);
  output <- t(c(coef(fit)[[1]], coef(fit)[[2]]))
  write.table(output, file=file.path(tpvdir,"output-monoexp-fit.txt"), append=TRUE, col.names=F, row.names=F)
  
  
  
  if(file.exists(file.path(tpvdir, "output-mixedbimono.txt")))
  {
    print("mixed biexp monoexp")
    fulloutput <- read.table(file.path(tpvdir,"output-mixedbimono.txt"), header=TRUE, fill = TRUE);
    n<-tail(grep("file",fulloutput[,1]),n=1)
    tpv <- read.table(file.path(tpvdir,"output-mixedbimono.txt"), header=TRUE, skip=ifelse(length(n),n,0), fill=TRUE); 
  if(output_pdf){
    pdf(file.path(tpvdir, paste(name, "-tpv-mixedbimono.pdf",sep="")), width=image_smallpdf_width, height=image_smallpdf_height, pointsize=7);
  }else{
    png(file.path(tpvdir, paste(name, "-tpv-mixedbimono.png",sep="")), width=image_width, height=image_height);
  }
  op <- par(mar = c(5,7,1,1) + 0.1) ## default is c(5,4,4,2) + 0.1
    plot(NULL, ylim=c(min(tpv$T1,tpv$T2[!is.na(tpv$T2)]),max(tpv$T2[!is.na(tpv$T2)],tpv$T1)), xlim=c(min(tpv$Voc),max(tpv$Voc)), log="y", xlab="Light bias (V)", ylab="", cex.axis=1.4, cex.lab=1.7, yaxt="n");#, main=paste(name, "TPV biexp and monoexp")
  title(ylab = "Small perturbation life-time (s)", cex.lab = 1.7, line = 4)
  eaxis(side=2,at=c(1e-10,1e-9,1e-8,1e-7,1e-6,1e-5,1e-4,1e-3,1e-2,0.1,1,10,100,1e3), cex.axis=1.4)
  minor.tick(nx=10)
    points(tpv$Voc, tpv$T1, pch=21, col=change.lightness(mycolors[1],0.5), bg=mycolors[1], cex=1.5)
    points(tpv$Voc, tpv$T2, pch=21, col=change.lightness(mycolors[1],0.5), cex=1.5)    
    #par(new = TRUE)
    #plot(NULL, ylim=c(min(tpv$A1,tpv$A2[!is.na(tpv$A2)]),max(tpv$A2[!is.na(tpv$A2)],tpv$A1)), xlim=c(min(tpv$Voc),max(tpv$Voc)), axes = FALSE, bty = "n", xlab = "", ylab = "")
    #points(tpv$Voc, tpv$A1, pch=20, col="red")
    #points(tpv$Voc, tpv$A2, pch=1, col="red")
    #axis(side=4, at = pretty(range(c(tpv$A1, tpv$A2))))
    #mtext("y0", side=4, line=3)
    graphics.off()
  }
}
