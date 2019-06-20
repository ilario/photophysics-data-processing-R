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

#title=gsub("-","\n\n",gsub("_"," ",name))
filename=gsub(",","",gsub(":","",name))

ylim=lim.TPVCE.lifetime
xlim=lim.TPVCE.charge
ylim_nogeom=lim.TPVCE.nogeom.lifetime
xlim_nogeom=lim.TPVCE.nogeom.charge

output=list()
output.nogeom=list()
output.nogeom.total=list()

library(robustbase)
library(RColorBrewer)
library(minpack.lm)
library(sfsmisc)
library(Hmisc)

q=1.60217662e-19

## Add an alpha value to a colour https://gist.github.com/mages/5339689
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


dirs <- list.dirs(recursive=FALSE)
dirs <- sub("./","",dirs)
#remove everything after the last underscore, then convert the remaining underscores to spaces
legendlist=sub("^0","",gsub("_"," ",sub("_((?!_).)+$","",dirs, perl=TRUE)))

# try to obtain the color from the file name
mycolors=gsub(".*-col_","",dirs[grepl("-col_", dirs)])
# if the color is not set, use the default one
if(!length(mycolors)){mycolors=brewer.pal(8,"Dark2")}

i <- 0
if(output_pdf){
  pdf(paste(filename,"-TPVCEs.pdf",sep=""), width=image_smallpdf_width, height=image_smallpdf_height, pointsize=7)
}else{
  png(paste(filename,"-TPVCEs.png",sep=""), width=image_width, height=image_height)
}
op <- par(mar = c(4,6,1,1) + 0.1) ## default is c(5,4,4,2) + 0.1
plot(NULL,xlim=xlim,ylim=ylim,xlab="", ylab="",log="xy", yaxt="n", xaxt="n")
#line is for introducing more space between label and axis
title(ylab = "Small perturbation life-time (s)", cex.lab = 1.7, line = 4)
title(xlab = bquote("Charge per area (C/cm"^"2"*")"), cex.lab = 1.7, line = 3)

eaxis(side=2,at=c(1e-11, 1e-10,1e-9,1e-8,1e-7,1e-6,1e-5,1e-4,1e-3,1e-2,0.1,1,10,100,1e3), cex.axis=1.4)
#xtick = 10^(floor(log10(xlim[2])))
#eaxis(side=1,at=c(1e-11, 1e-10,1e-9,1e-8,1e-7,1e-6,1e-5,1e-4,1e-3,1e-2,0.1,1,10,100,1e3), cex.axis=1.4)
eaxis(side=1, cex.axis=1.4)

#for x linear
#eaxis(side=1,at=seq(0, xlim[2], xtick), cex.axis=1.4)
#minor.tick(nx=10)

lapply(dirs, function(x) {print(x);
  subdirs <- list.dirs(path=x, recursive=F)
  subdirs.ce <- subdirs[grep("ce", subdirs, ignore.case=T)]
  subdirs.tpv <- subdirs[grep("tpv", subdirs, ignore.case=T)]
  a <- read.table(file.path(subdirs.ce,"outputChargeDensityCE.txt"),header=T,stringsAsFactors=F)
  lo <- loess(a$ChargeDensityCE~a$Voc,span=0.9)
  tryCatch({
    lin <- lm(ChargeDensityCE ~ 0 + Voc, data=a)
    expfit <- lin;
  }, error=function(e) {cat("FAILED LINEAR FIT ", e$message, "\n")});
  tryCatch({
    expfit <- nlsLM(ChargeDensityCE~ exp(B)*Voc+exp(C)*(exp(exp(D)*Voc)-1), start=list(B=log(max(a$ChargeDensityCE)/max(a$Voc)),C=log(1e-10),D=2), data=a)
  }, error=function(e) {cat("FAILED non-robust FIT ", e$message, "\n")});
  tryCatch({
    expfit <- nlsLM(ChargeDensityCE~ exp(B)*Voc+exp(C)*(exp(exp(D)*Voc)-1), start=list(B=log(coef(lin)[[1]]),C=log(1e-10),D=2), data=a)
  }, error=function(e) {cat("FAILED second non-robust FIT ", e$message, "\n")});
  tryCatch({
    expfit <- nlrob(ChargeDensityCE~ exp(B)*Voc+exp(C)*(exp(exp(D)*Voc)-1), start=list(B=coef(expfit)["B"],C=coef(expfit)["C"],D=coef(expfit)["D"]), data=a)
  }, error=function(e) {cat("FAILED robust FIT ", e$message, "\n")});
  
  
  temp=dev.cur()
  dev.set(temp+1)
  png(paste("debug_tpvce_", x, ".png",sep=""))
  plot(a$Voc,a$ChargeDensityCE)
  lines(a$Voc, predict(expfit))
  dev.off(temp+1)
  dev.set(temp)
  
  
  fulloutput <- read.table(file.path(subdirs.tpv,"output-robustmonoexp.txt"), header=TRUE);
  n<-tail(grep("file",fulloutput[,1]),n=1)
  tpv <- read.table(file.path(subdirs.tpv,"output-robustmonoexp.txt"), header=TRUE, skip=ifelse(length(n),n,0));
  #importante che la variabile in new abbia lo stesso nome di quella fittata
  new <- data.frame(Voc = tpv$Voc)
  charge <- (predict(lo, tpv$Voc) + predict(expfit, new))/2
  new2 <- data.frame(Voc = tpv$Voc[is.na(charge)])
  charge[is.na(charge)] <- predict(expfit,new2)
  output[[paste("Charge",sub("nm","",sub("_.*","",sub("^0","",x))),sep="")]] <<- signif(charge,5)
  output[[sub("_.*","",sub("^0","",x))]] <<- signif(tpv$T,5)
  points(charge, tpv$T, bg=add.alpha(mycolors[i+1],0.5), col=change.lightness(mycolors[i+1],0.5), cex=1.5, pch=21+(i%%5));
  index_shown_charge = which(charge >= xlim[1] & charge <= xlim[2])
  shown_charge = charge[index_shown_charge]
  shown_T = tpv$T[index_shown_charge]
  
  weights = (min(shown_T)/shown_T)^2
  
  if(length(shown_T) < 11 || length(shown_charge) < 11){stop("TPVCE: you need wider plot limits!")}
  # set to zero all but last 10 points' weight
  weights[1:(length(weights)-10)] = 0
  
  #just in case...
  rm(powerlaw)
  

  
  j=1
  while(!exists("powerlaw") && j < 1000){
    j <- j + 0.1
    #start <- list(y0=log(5e-10*runif(1,1/j,j)), A=log(1e-105*runif(1,1/j,j)), alpha=-13*runif(1,1/j,j))
    start <- list(A=log(q^-11*1e-105*runif(1,1/j,j)), alpha=-11*runif(1,1/j,j))
    tryCatch({
      #powerlaw <- nlsLM(shown_T~exp(y0)+exp(A)*shown_charge^alpha, start=start, weights=weights)
      powerlaw <- nlsLM(shown_T~exp(A)*(shown_charge/q)^alpha, start=start, weights=weights)
      #check convergence and sum the p-values
    }, error=function(e) {cat("FAILED POWERLAW FIT ", e$message, "\n");
    })
    #summary fails if the fit was done on no data, with some chol2inv error
    if(exists("powerlaw")){
      print("Checking powerlaw result")
      if(!summary(powerlaw)$convInfo$isConv || sum(coef(summary(powerlaw))[,"Pr(>|t|)"]) > 1){
        rm(powerlaw)
      }
    }
  }
  if(exists("powerlaw")){
    lines(shown_charge, predict(powerlaw, shown_charge), lwd=2, col=change.lightness(mycolors[i+1],0.5))
    capture.output(summary(powerlaw), file=paste(x, "-tpvce-fit.txt", sep=""),  append=TRUE);
  }
  
  i <<- i+1
})
legend(x="topright",inset=0.05,legendlist,pch=seq(21,25), pt.bg=mycolors, lwd=2, pt.lwd=1.5, pt.cex=2, col=change.lightness(mycolors,0.5),cex=1.5, title=title,bg="gray90", bty="n")
graphics.off()

#reset the plotting margins
par(op)

maxlength = max(sapply(output,length))
output = lapply(output, function(x){length(x)=maxlength; print(x)})
output = as.data.frame(output,check.names=FALSE)
write.table(output, file=paste(filename,"-TPVCEs.csv",sep=""), row.names=FALSE, na="", sep=",")






if(output_pdf){
  pdf(paste(filename,"-TPVCEs-nogeom.pdf",sep=""), width=image_smallpdf_width, height=image_smallpdf_height, pointsize=7)
}else{
  png(paste(filename,"-TPVCEs-nogeom.png",sep=""), width=image_width, height=image_height)
}

op <- par(mar = c(4,6,1,1) + 0.1) ## default is c(5,4,4,2) + 0.1
plot(NULL,xlim=xlim_nogeom,ylim=ylim_nogeom,xlab="", ylab="",log="xy", yaxt="n", xaxt="n")

#line is for introducing more space between label and axis
title(ylab = "Small perturbation life-time (s)", cex.lab = 1.7, line = 4)
title(xlab = bquote("Charge per area (C/cm"^"2"*")"), cex.lab = 1.7, line = 3)
eaxis(side=2,at=c(1e-16, 1e-15, 1e-14, 1e-13, 1e-12, 1e-11, 1e-10,1e-9,1e-8,1e-7,1e-6,1e-5,1e-4,1e-3,1e-2,0.1,1,10,100,1e3), cex.axis=1.4)
#xtick = 10^(floor(log10(xlim_nogeom[2])))
eaxis(side=1,at=c(1e-16, 1e-15, 1e-14, 1e-13, 1e-12, 1e-11, 1e-10,1e-9,1e-8,1e-7,1e-6,1e-5,1e-4,1e-3,1e-2,0.1,1,10,100,1e3), cex.axis=1.4)

# preallocate recombination orders array
recombination_orders_nogeom = numeric(length(dirs))
i <- 0
lapply(dirs, function(x) {print(x);
  subdirs <- list.dirs(path=x, recursive=F)
  subdirs.ce <- subdirs[grep("ce", subdirs, ignore.case=T)]
  subdirs.tpv <- subdirs[grep("tpv", subdirs, ignore.case=T)]
  a <- read.table(file.path(subdirs.ce,"outputChargeDensityCE.txt"),header=T,stringsAsFactors=F)
  lo <- loess(a$ChargeDensityCE~a$Voc,span=0.9)
  tryCatch({
    lin <- lm(ChargeDensityCE ~ 0 + Voc, data=a)
    expfit <- lin;
  }, error=function(e) {cat("FAILED LINEAR FIT ", e$message, "\n")});
  tryCatch({
    expfit <- nlsLM(ChargeDensityCE~ exp(B)*Voc+exp(C)*(exp(exp(D)*Voc)-1), start=list(B=log(max(a$ChargeDensityCE)/max(a$Voc)),C=log(1e-10),D=2), data=a)
  }, error=function(e) {cat("FAILED non-robust FIT ", e$message, "\n")});
  tryCatch({
    expfit <- nlsLM(ChargeDensityCE~ exp(B)*Voc+exp(C)*(exp(exp(D)*Voc)-1), start=list(B=log(coef(lin)[[1]]),C=log(1e-10),D=2), data=a)
  }, error=function(e) {cat("FAILED second non-robust FIT ", e$message, "\n")});
  tryCatch({
    expfit <- nlrob(ChargeDensityCE~ exp(B)*Voc+exp(C)*(exp(exp(D)*Voc)-1), start=list(B=coef(expfit)["B"],C=coef(expfit)["C"],D=coef(expfit)["D"]), data=a)
  }, error=function(e) {cat("FAILED robust FIT ", e$message, "\n")});
  
  fulloutput <- read.table(file.path(subdirs.tpv,"output-robustmonoexp.txt"), header=TRUE);
  n<-tail(grep("file",fulloutput[,1]),n=1)
  tpv <- read.table(file.path(subdirs.tpv,"output-robustmonoexp.txt"), header=TRUE, skip=ifelse(length(n),n,0));
  
  charge_nogeom <- exp(coef(expfit)[2])*(exp(exp(coef(expfit)[3])*tpv$Voc)-1)
  output.nogeom[[paste("Charge",sub("nm","",sub("_.*","",sub("^0","",x))),sep="")]] <<- signif(charge_nogeom,5)
  output.nogeom[[sub("_.*","",sub("^0","",x))]] <<- signif(tpv$T,5)

  index_shown_charge_nogeom = which(charge_nogeom >= xlim_nogeom[1] & charge_nogeom <= xlim_nogeom[2])
  shown_charge_nogeom = charge_nogeom[index_shown_charge_nogeom]
  shown_T_nogeom = tpv$T[index_shown_charge_nogeom]
  
  weights_nogeom = (min(shown_T_nogeom)/shown_T_nogeom)^2
  
  if(length(shown_T_nogeom) < 11 || length(shown_charge_nogeom) < 11){graphics.off(); stop("TPVCE_nogeom: you need wider plot limits!")}
  # set to zero all but last 10 points' weight
  weights_nogeom[1:(length(weights_nogeom)-10)] = 0
  
  #just in case...
  rm(powerlaw_nogeom)

  j=1
  while(!exists("powerlaw_nogeom") && j < 1000){
    j <- j + 0.1
    #start_nogeom <- list(y0=log(5e-7*runif(1,1/j,j)), A=log(1e-28*runif(1,1/j,j)), alpha=-3.2*runif(1,1/j,j))
    start_nogeom <- list(A=log(1e19*runif(1,1/j,j)), alpha=-3.2*runif(1,1/j,j))
    
    tryCatch({
      #powerlaw_nogeom <- nlsLM(shown_T_nogeom~exp(y0)+exp(A)*shown_charge_nogeom^alpha, start=start_nogeom, weights=weights_nogeom)
      powerlaw_nogeom <- nlsLM(shown_T_nogeom~exp(A)*(shown_charge_nogeom/q)^alpha, start=start_nogeom, weights=weights_nogeom)
      #powerlaw_nogeom_rob <- nlrob(shown_T_nogeom~exp(y0)+exp(A)*shown_charge_nogeom^alpha, start=list(y0=coef(powerlaw_nogeom)["y0"], A=coef(powerlaw_nogeom)["A"], alpha=coef(powerlaw_nogeom)["alpha"]), data=data.frame(shown_T_nogeom=shown_T_nogeom, shown_charge_nogeom=shown_charge_nogeom), weights=weights_nogeom)

    }, error=function(e) {cat("FAILED POWERLAW FIT ", e$message, "\n");
    })
    #summary fails if the fit was done on no data, with some chol2inv error
    if(exists("powerlaw_nogeom")){
      print("Checking powerlaw result")
      #check convergence and sum the p-values
      if(!summary(powerlaw_nogeom)$convInfo$isConv || sum(coef(summary(powerlaw_nogeom))[,"Pr(>|t|)"]) > 1){
        rm(powerlaw_nogeom)
      }
    }
  }
  if(exists("powerlaw_nogeom")){
    lines(shown_charge_nogeom, predict(powerlaw_nogeom, shown_charge_nogeom), lwd=2, col=add.alpha(change.lightness(mycolors[i+1],0.8),0.8))
    capture.output(summary(powerlaw_nogeom), file=paste(x, "-tpvce-nogeom-fit.txt", sep=""),  append=TRUE);
  }
  recombination_orders_nogeom[i+1] <<- 1-coef(powerlaw_nogeom)["alpha"]

  points(charge_nogeom, tpv$T, bg=add.alpha(mycolors[i+1],0.5), col=change.lightness(mycolors[i+1],0.5), cex=1.5, pch=21+(i%%5));
  
  i <<- i+1
})
legend(x="bottomleft",inset=0.02,legendlist,pch=seq(21,25), pt.bg=mycolors, lwd=2, pt.lwd=1.5, pt.cex=2, col=change.lightness(mycolors,0.5),cex=1.5, title=title,bg="gray90", bty="n")
graphics.off()

#reset the plotting margins
par(op)

output.nogeom = lapply(output.nogeom, function(x){length(x)=maxlength; print(x)})
output.nogeom = as.data.frame(output.nogeom,check.names=FALSE)
write.table(output.nogeom, file=paste(filename,"-TPVCEs-nogeom.csv",sep=""), row.names=FALSE, na="", sep=",")







if(output_pdf){
  pdf(paste(filename,"-TPVCEs-nogeom_total.pdf",sep=""), width=image_smallpdf_width, height=image_smallpdf_height, pointsize=7)
}else{
  png(paste(filename,"-TPVCEs-nogeom_total.png",sep=""), width=image_width, height=image_height)
}

op <- par(mar = c(4,6,1,1) + 0.1) ## default is c(5,4,4,2) + 0.1
plot(NULL,xlim=xlim_nogeom,ylim=ylim_nogeom,xlab="", ylab="",log="xy", yaxt="n", xaxt="n")

title(ylab = "Total carrier life-time (s)", cex.lab = 1.7, line = 4)
title(xlab = bquote("Charge per area (C/cm"^"2"*")"), cex.lab = 1.7, line = 3)
eaxis(side=2,at=c(1e-16, 1e-15, 1e-14, 1e-13, 1e-12, 1e-11, 1e-10,1e-9,1e-8,1e-7,1e-6,1e-5,1e-4,1e-3,1e-2,0.1,1,10,100,1e3), cex.axis=1.4)
#xtick = 10^(floor(log10(xlim_nogeom[2])))
eaxis(side=1,at=c(1e-16, 1e-15, 1e-14, 1e-13, 1e-12, 1e-11, 1e-10,1e-9,1e-8,1e-7,1e-6,1e-5,1e-4,1e-3,1e-2,0.1,1,10,100,1e3), cex.axis=1.4)

i <- 0
lapply(dirs, function(x) {print(x);
  subdirs <- list.dirs(path=x, recursive=F)
  subdirs.ce <- subdirs[grep("ce", subdirs, ignore.case=T)]
  subdirs.tpv <- subdirs[grep("tpv", subdirs, ignore.case=T)]
  a <- read.table(file.path(subdirs.ce,"outputChargeDensityCE.txt"),header=T,stringsAsFactors=F)
  lo <- loess(a$ChargeDensityCE~a$Voc,span=0.9)
  tryCatch({
    lin <- lm(ChargeDensityCE ~ 0 + Voc, data=a)
    expfit <- lin;
  }, error=function(e) {cat("FAILED LINEAR FIT ", e$message, "\n")});
  tryCatch({
    expfit <- nlsLM(ChargeDensityCE~ exp(B)*Voc+exp(C)*(exp(exp(D)*Voc)-1), start=list(B=log(max(a$ChargeDensityCE)/max(a$Voc)),C=log(1e-10),D=2), data=a)
  }, error=function(e) {cat("FAILED non-robust FIT ", e$message, "\n")});
  tryCatch({
    expfit <- nlsLM(ChargeDensityCE~ exp(B)*Voc+exp(C)*(exp(exp(D)*Voc)-1), start=list(B=log(coef(lin)[[1]]),C=log(1e-10),D=2), data=a)
  }, error=function(e) {cat("FAILED second non-robust FIT ", e$message, "\n")});
  tryCatch({
    expfit <- nlrob(ChargeDensityCE~ exp(B)*Voc+exp(C)*(exp(exp(D)*Voc)-1), start=list(B=coef(expfit)["B"],C=coef(expfit)["C"],D=coef(expfit)["D"]), data=a)
  }, error=function(e) {cat("FAILED robust FIT ", e$message, "\n")});
  
  fulloutput <- read.table(file.path(subdirs.tpv,"output-robustmonoexp.txt"), header=TRUE);
  n<-tail(grep("file",fulloutput[,1]),n=1)
  tpv <- read.table(file.path(subdirs.tpv,"output-robustmonoexp.txt"), header=TRUE, skip=ifelse(length(n),n,0));
  
  tpv$Ttotal = tpv$T * recombination_orders_nogeom[i+1]
  
  charge_nogeom <- exp(coef(expfit)[2])*(exp(exp(coef(expfit)[3])*tpv$Voc)-1)
  output.nogeom.total[[paste("Charge",sub("nm","",sub("_.*","",sub("^0","",x))),sep="")]] <<- signif(charge_nogeom,5)
  output.nogeom.total[[sub("_.*","",sub("^0","",x))]] <<- signif(tpv$Ttotal,5)

  index_shown_charge_nogeom = which(charge_nogeom >= xlim_nogeom[1] & charge_nogeom <= xlim_nogeom[2])
  shown_charge_nogeom = charge_nogeom[index_shown_charge_nogeom]
  shown_T_nogeom = tpv$Ttotal[index_shown_charge_nogeom]

  weights_nogeom = (min(shown_T_nogeom)/shown_T_nogeom)^2
  
  if(length(shown_T_nogeom) < 11 || length(shown_charge_nogeom) < 11){graphics.off(); stop("TPVCE_nogeom_total: you need wider plot limits!")}
  # set to zero all but last 10 points' weight
  weights_nogeom[1:(length(weights_nogeom)-10)] = 0
  
  #just in case...
  rm(powerlaw_nogeom)
  

  j=1
  while(!exists("powerlaw_nogeom") && j < 1000){
    j <- j + 0.1
    #start_nogeom <- list(y0=log(5e-7*runif(1,1/j,j)), A=log(1e-28*runif(1,1/j,j)), alpha=-3.2*runif(1,1/j,j))
    start_nogeom <- list(A=log(1e19*runif(1,1/j,j)), alpha=-3.2*runif(1,1/j,j))
    
    tryCatch({
      #powerlaw_nogeom <- nlsLM(shown_T_nogeom~exp(y0)+exp(A)*shown_charge_nogeom^alpha, start=start_nogeom, weights=weights_nogeom)
      powerlaw_nogeom <- nlsLM(shown_T_nogeom~exp(A)*(shown_charge_nogeom/q)^alpha, start=start_nogeom, weights=weights_nogeom)
      #powerlaw_nogeom_rob <- nlrob(shown_T_nogeom~exp(y0)+exp(A)*shown_charge_nogeom^alpha, start=list(y0=coef(powerlaw_nogeom)["y0"], A=coef(powerlaw_nogeom)["A"], alpha=coef(powerlaw_nogeom)["alpha"]), data=data.frame(shown_T_nogeom=shown_T_nogeom, shown_charge_nogeom=shown_charge_nogeom), weights=weights_nogeom)
    }, error=function(e) {cat("FAILED POWERLAW FIT ", e$message, "\n");
    })
    #summary fails if the fit was done on no data, with some chol2inv error
    if(exists("powerlaw_nogeom")){
      print("Checking powerlaw result")
#check convergence and sum the p-values
      if(!summary(powerlaw_nogeom)$convInfo$isConv || sum(coef(summary(powerlaw_nogeom))[,"Pr(>|t|)"]) > 1){
        rm(powerlaw_nogeom)
      }
    }
  }
  if(exists("powerlaw_nogeom")){
    lines(shown_charge_nogeom, predict(powerlaw_nogeom, shown_charge_nogeom), lwd=2, col=add.alpha(change.lightness(mycolors[i+1],0.8),0.8))
    capture.output(summary(powerlaw_nogeom), file=paste(x, "-tpvce-nogeom-total-fit.txt", sep=""),  append=TRUE);
  }
  
  points(charge_nogeom, tpv$Ttotal, bg=add.alpha(mycolors[i+1],0.5), col=change.lightness(mycolors[i+1],0.5), cex=1.5, pch=21+(i%%5));
  
  i <<- i+1
})
legend(x="bottomleft",inset=0.02,legendlist,pch=seq(21,25), pt.bg=mycolors, lwd=2, pt.lwd=1.5, pt.cex=2, col=change.lightness(mycolors,0.5),cex=1.5, title=title,bg="gray90", bty="n")
graphics.off()

#reset the plotting margins
par(op)


output.nogeom.total = lapply(output.nogeom.total, function(x){length(x)=maxlength; print(x)})
output.nogeom.total = as.data.frame(output.nogeom.total,check.names=FALSE)
write.table(output.nogeom.total, file=paste(filename,"-TPVCEs-nogeom-total.csv",sep=""), row.names=FALSE, na="", sep=",")


