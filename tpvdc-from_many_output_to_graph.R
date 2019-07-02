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
ylim_nogeom=lim.TPVDC.nogeom.lifetime
xlim_nogeom=lim.TPVDC.nogeom.charge

output=list()
output.nogeom=list()
output.nogeom.total=list()

library(robustbase)
library(RColorBrewer)
library(minpack.lm)
library(sfsmisc)
library(Hmisc)

q=1.60217662e-19

write.table(t(c("file","RecConstant","RecOrder","CEexpprefactor","EquilibriumLifetime")), file="outputTPVDC-nogeom-fit.txt", append=FALSE, col.names=F, row.names=F);
write.table(t(c("file","RecConstant","RecOrder","CEexpprefactor","EquilibriumLifetime")), file="outputTPVDC-nogeom-total-fit.txt", append=FALSE, col.names=F, row.names=F);

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
  pdf(paste(filename,"-TPVDCs.pdf",sep=""), width=image_smallpdf_width, height=image_smallpdf_height, pointsize=7)
}else{
  png(paste(filename,"-TPVDCs.png",sep=""), width=image_width, height=image_height)
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
  
  filex <- file.path(subdirs.tpv, "output-robustmonoexp.txt")
  
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
  points(charge, tpv$T, bg=add.alpha(mycolors[i+1],0.5), col=change.lightness(mycolors[i+1],0.5), cex=1.5, pch=21+(i%%5));
  i <<- i+1
})
legend(x="topright",inset=0.05,legendlist,pch=seq(21,25), pt.bg=mycolors, lwd=2, pt.lwd=1.5, pt.cex=2, col=change.lightness(mycolors,0.5), cex=1.5, title=title, bg="gray90", bty="n")
graphics.off()
#reset the plotting margins
par(op)

maxlength = max(sapply(output,length))
output = lapply(output, function(x){length(x)=maxlength; print(x)})
output = as.data.frame(output,check.names=FALSE)
write.table(output, file=paste(filename,"-TPVDCs.csv",sep=""), row.names=FALSE, na="", sep=",")

# preallocate recombination orders array
recombination_orders_nogeom = numeric(length(dirs))
i <- 0
if(output_pdf){
  pdf(paste(filename,"-TPVDCs-nogeom.pdf",sep=""), width=image_smallpdf_width, height=image_smallpdf_height, pointsize=7)
}else{
  png(paste(filename,"-TPVDCs-nogeom.png",sep=""), width=image_width, height=image_height)
}
op <- par(mar = c(4,6,1,1) + 0.1) ## default is c(5,4,4,2) + 0.1
plot(NULL,xlim=xlim_nogeom,ylim=ylim_nogeom,xlab="", ylab="",log="xy", yaxt="n", xaxt="n")

#line is for introducing more space between label and axis
title(ylab = "Small perturbation life-time (s)", cex.lab = 1.7, line = 4)
title(xlab = bquote("Charge per area (C/cm"^"2"*")"), cex.lab = 1.7, line = 3)
eaxis(side=2,at=c(1e-16, 1e-15, 1e-14, 1e-13, 1e-12, 1e-11, 1e-10,1e-9,1e-8,1e-7,1e-6,1e-5,1e-4,1e-3,1e-2,0.1,1,10,100,1e3), cex.axis=1.4)
#xtick = 10^(floor(log10(xlim_nogeom[2])))
eaxis(side=1,at=c(1e-16, 1e-15, 1e-14, 1e-13, 1e-12, 1e-11, 1e-10,1e-9,1e-8,1e-7,1e-6,1e-5,1e-4,1e-3,1e-2,0.1,1,10,100,1e3), cex.axis=1.4)

lapply(dirs, function(x) {print(x);
  subdirs <- list.dirs(path=x, recursive=F)
  subdirs.tpv <- subdirs[grep("tpv", subdirs, ignore.case=T)]
  a <- read.table(file.path(x,"outputDCcharge-nogeom.txt"),header=T,stringsAsFactors=F)

  rm(expfit)
  startlist=list(C=1e-12,D=9)
  tryCatch({
    expfit <- nlsLM(ChargeDensityDC~ C*(exp(D*Voc)-1), start=startlist, data=a)
    tryCatch({
      expfit <- nlrob(ChargeDensityDC~ C*(exp(D*Voc)-1), start=list(C=coef(expfit)[1],D=coef(expfit)[2]), data=a)
    }, error=function(e) {print("FAILED nogeom FIT ROBUST")});
  }, error=function(e) {print("FAILED nogeom FIT non-robust")});
  
  temp=dev.cur()
  dev.set(temp+1)
  png(paste("debug_tpvdc_", x, ".png",sep=""))
  plot(a$Voc,a$ChargeDensityDC)
  lines(a$Voc, predict(expfit))
  dev.off(temp+1)
  dev.set(temp)
  
  filex <- file.path(subdirs.tpv, "output-robustmonoexp.txt")
  fulloutput <- read.table(filex, header=TRUE);
  n<-tail(grep("file",fulloutput[,1]),n=1)
  tpv <- read.table(filex, header=TRUE, skip=ifelse(length(n),n,0));
  #importante che la variabile in new abbia lo stesso nome di quella fittata
  new <- data.frame(Voc = tpv$Voc)
  
  charge_nogeom <- predict(expfit,new)
#  lo <- loess(a$ChargeDensityDC~a$Voc,span=0.9)
#  charge_nogeom <- (predict(lo,tpv$Voc)+predict(expfit,new))/2
#  new2 <- data.frame(Voc = tpv$Voc[is.na(charge_nogeom)])
#  charge_nogeom[is.na(charge_nogeom)] <- predict(expfit,new2)
  
  output.nogeom[[paste("Charge",sub("nm","",sub("_.*","",sub("^0","",x))),sep="")]] <<- signif(charge_nogeom,5)
  output.nogeom[[sub("_.*","",sub("^0","",x))]] <<- signif(tpv$T,5)

  index_shown_charge_nogeom = which(charge_nogeom >= xlim_nogeom[1] & charge_nogeom <= xlim_nogeom[2])
  shown_charge_nogeom = charge_nogeom[index_shown_charge_nogeom]
  shown_T_nogeom = tpv$T[index_shown_charge_nogeom]
  
  weights_nogeom = (min(shown_T_nogeom)/shown_T_nogeom)^2
  # set to zero all but last 10 points' weight
  weights_nogeom = c(integer(length(weights_nogeom)-10), weights_nogeom[length(weights_nogeom)-(9:0)])

  
  #just in case...
  rm(powerlaw_nogeom)
  
  if(length(shown_T_nogeom) < 4 || length(shown_charge_nogeom) < 4){graphics.off(); stop("TPVDC_nogeom: you need wider plot limits!")}
  j=1
  while(!exists("powerlaw_nogeom") && j < 1000){
    j <- j + 0.1
    start_nogeom <- list(A=log(1e19*runif(1,1/j,j)), alpha=-3.2*runif(1,1/j,j))
    
    tryCatch({
      powerlaw_nogeom <- nlsLM(shown_T_nogeom~exp(A)*(shown_charge_nogeom/q)^alpha, start=start_nogeom, weights=weights_nogeom)
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
    capture.output(summary(powerlaw_nogeom), file=paste(x, "-tpvdc-nogeom-fit.txt", sep=""),  append=TRUE);
    
    # DCexpprefactor
    n_eq = coef(expfit)["C"]/q
    # RecOrder
    Phi = 1-coef(powerlaw_nogeom)["alpha"]
    # RecConstant
    k = 1/exp(coef(powerlaw_nogeom)["A"])
    # EquilibriumLifetime
    tau0 = (n_eq^(1-Phi))/k
    outputTPVDC_nogeom_fit <- data.frame(x, k, Phi, n_eq, tau0);
    write.table(outputTPVDC_nogeom_fit, file="outputTPVDC-nogeom-fit.txt", append=TRUE, col.names=F, row.names=F, quote=F);
    
    recombination_orders_nogeom[i+1] <<- 1-coef(powerlaw_nogeom)["alpha"]
  }else{
    recombination_orders_nogeom[i+1] <<- 0
  }
  
  points(charge_nogeom, tpv$T, bg=add.alpha(mycolors[i+1],0.5), col=change.lightness(mycolors[i+1],0.5), cex=1.5, pch=21+(i%%5));
  
  i <<- i+1
})
legend(x="topright",inset=0.05,legendlist,pch=seq(21,25), pt.bg=mycolors, lwd=2, pt.lwd=1.5, pt.cex=2, col=change.lightness(mycolors,0.5),cex=1.5, title=title,bg="gray90", bty="n")
graphics.off()
#reset the plotting margins
par(op)

output.nogeom = lapply(output.nogeom, function(x){length(x)=maxlength; print(x)})
output.nogeom = as.data.frame(output.nogeom,check.names=FALSE)
write.table(output.nogeom, file=paste(filename,"-TPVDCs-nogeom.csv",sep=""), row.names=FALSE, na="", sep=",")



i <- 0
if(output_pdf){
  pdf(paste(filename,"-TPVDCs-nogeom_total.pdf",sep=""), width=image_smallpdf_width, height=image_smallpdf_height, pointsize=7)
}else{
  png(paste(filename,"-TPVDCs-nogeom_total.png",sep=""), width=image_width, height=image_height)
}
op <- par(mar = c(4,6,1,1) + 0.1) ## default is c(5,4,4,2) + 0.1
plot(NULL,xlim=xlim_nogeom,ylim=ylim_nogeom,xlab="", ylab="",log="xy", yaxt="n", xaxt="n")

#line is for introducing more space between label and axis
title(ylab = "Total carrier life-time (s)", cex.lab = 1.7, line = 4)
title(xlab = bquote("Charge per area (C/cm"^"2"*")"), cex.lab = 1.7, line = 3)
eaxis(side=2,at=c(1e-16, 1e-15, 1e-14, 1e-13, 1e-12, 1e-11, 1e-10,1e-9,1e-8,1e-7,1e-6,1e-5,1e-4,1e-3,1e-2,0.1,1,10,100,1e3), cex.axis=1.4)
#xtick = 10^(floor(log10(xlim_nogeom[2])))
eaxis(side=1,at=c(1e-16, 1e-15, 1e-14, 1e-13, 1e-12, 1e-11, 1e-10,1e-9,1e-8,1e-7,1e-6,1e-5,1e-4,1e-3,1e-2,0.1,1,10,100,1e3), cex.axis=1.4)

lapply(dirs, function(x) {print(x);
  subdirs <- list.dirs(path=x, recursive=F)
  subdirs.tpv <- subdirs[grep("tpv", subdirs, ignore.case=T)]
  a <- read.table(file.path(x,"outputDCcharge-nogeom.txt"),header=T,stringsAsFactors=F)
  startlist=list(C=1e-10,D=8)
  tryCatch({
    expfit <- nlsLM(ChargeDensityDC~ C*(exp(D*Voc)-1), start=startlist, data=a)
    tryCatch({
      expfit <- nlrob(ChargeDensityDC~ C*(exp(D*Voc)-1), start=list(C=coef(expfit)[1],D=coef(expfit)[2]), data=a)
    }, error=function(e) {print("FAILED nogeom FIT ROBUST")});
  }, error=function(e) {print("FAILED nogeom FIT non-robust")});
  
  filex <- file.path(subdirs.tpv, "output-robustmonoexp.txt")
  fulloutput <- read.table(filex, header=TRUE);
  n<-tail(grep("file",fulloutput[,1]),n=1)
  tpv <- read.table(filex, header=TRUE, skip=ifelse(length(n),n,0));
  
  tpv$Ttotal = tpv$T * recombination_orders_nogeom[i+1]
  
  #importante che la variabile in new abbia lo stesso nome di quella fittata
  new <- data.frame(Voc = tpv$Voc)
  
  charge_nogeom <- predict(expfit,new)
  
#    lo <- loess(a$ChargeDensityDC~a$Voc,span=0.9)
#  charge_nogeom <- (predict(lo,tpv$Voc)+predict(expfit,new))/2
#  new2 <- data.frame(Voc = tpv$Voc[is.na(charge_nogeom)])
#  charge_nogeom[is.na(charge_nogeom)] <- predict(expfit,new2)
  
  output.nogeom.total[[paste("Charge",sub("nm","",sub("_.*","",sub("^0","",x))),sep="")]] <<- signif(charge_nogeom,5)
  output.nogeom.total[[sub("_.*","",sub("^0","",x))]] <<- signif(tpv$Ttotal,5)

  index_shown_charge_nogeom = which(charge_nogeom >= xlim_nogeom[1] & charge_nogeom <= xlim_nogeom[2])
  shown_charge_nogeom = charge_nogeom[index_shown_charge_nogeom]
  shown_T_nogeom = tpv$Ttotal[index_shown_charge_nogeom]
  
  weights_nogeom = (min(shown_T_nogeom)/shown_T_nogeom)^2
  # set to zero all but last 10 points' weight
  weights_nogeom = c(integer(length(weights_nogeom)-10), weights_nogeom[length(weights_nogeom)-(9:0)])
  
  #just in case...
  rm(powerlaw_nogeom)
  
  if(length(shown_T_nogeom) < 4 || length(shown_charge_nogeom) < 4){graphics.off(); stop("TPVDC_nogeom: you need wider plot limits!")}
  j=1
  while(!exists("powerlaw_nogeom") && j < 1000){
    j <- j + 0.1
    start_nogeom <- list(A=log(1e19*runif(1,1/j,j)), alpha=-3.2*runif(1,1/j,j))
    
    tryCatch({
      powerlaw_nogeom <- nlsLM(shown_T_nogeom~exp(A)*(shown_charge_nogeom/q)^alpha, start=start_nogeom, weights=weights_nogeom)
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
    capture.output(summary(powerlaw_nogeom), file=paste(x, "-tpvdc-nogeom-total-fit.txt", sep=""),  append=TRUE);
    
    # DCexpprefactor
    n_eq = coef(expfit)["C"]/q
    # RecOrder
    Phi = 1-coef(powerlaw_nogeom)["alpha"]
    # RecConstant
    k = 1/exp(coef(powerlaw_nogeom)["A"])
    # EquilibriumLifetime
    tau0 = (n_eq^(1-Phi))/k
    outputTPVDC_nogeom_total_fit <- data.frame(x, k, Phi, n_eq, tau0);
    write.table(outputTPVDC_nogeom_total_fit, file="outputTPVDC-nogeom-total-fit.txt", append=TRUE, col.names=F, row.names=F, quote=F);
  }

  points(charge_nogeom, tpv$Ttotal, bg=add.alpha(mycolors[i+1],0.5), col=change.lightness(mycolors[i+1],0.5), cex=1.5, pch=21+(i%%5));
  
  i <<- i+1
})
legend(x="bottomleft",inset=0,legendlist,pch=seq(21,25), pt.bg=mycolors, lwd=2, pt.lwd=1.5, pt.cex=2, col=change.lightness(mycolors,0.5),cex=1.5, title=title,bg="gray90", bty="n")
graphics.off()
#reset the plotting margins
par(op)

output.nogeom.total = lapply(output.nogeom.total, function(x){length(x)=maxlength; print(x)})
output.nogeom.total = as.data.frame(output.nogeom.total,check.names=FALSE)
write.table(output.nogeom.total, file=paste(filename,"-TPVDCs-nogeom-total.csv",sep=""), row.names=FALSE, na="", sep=",")