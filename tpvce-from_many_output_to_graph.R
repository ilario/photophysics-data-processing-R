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

library(robustbase)
library(RColorBrewer)
library(minpack.lm)
library(sfsmisc)
library(Hmisc)

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
legend=sub("_.*","",sub("^0","",dirs))

# try to obtain the color from the file name
mycolors=gsub(".*-col_","",dirs[grepl("-col_", dirs)])
# if the color is not set, use the default one
if(!length(mycolors)){mycolors=brewer.pal(8,"Dark2")}

i <- 0
if(output_pdf){
	pdf(paste(filename,"-TPVCEs.pdf",sep=""), width=image_bigpdf_width, height=image_bigpdf_height, pointsize=7)
}else{
	jpeg(quality=98, paste(filename,"-TPVCEs.jpg",sep=""), width=image_width, height=image_height)
}
op <- par(mar = c(5,6,1,1) + 0.1) ## default is c(5,4,4,2) + 0.1
plot(1,xlim=xlim,ylim=ylim,xlab="", ylab="",log="y", yaxt="n", xaxt="n")
#line is for introducing more space between label and axis
title(ylab = "Charge carrier lifetime (s)", cex.lab = 1.7, line = 4)
title(xlab = bquote("Charge (C/cm"^"2"*")"), cex.lab = 1.7, line = 4)

eaxis(side=2,at=c(1e-11, 1e-10,1e-9,1e-8,1e-7,1e-6,1e-5,1e-4,1e-3,1e-2,0.1,1,10,100,1e3), cex.axis=1.4)
xtick = 10^(floor(log10(xlim[2])))
#eaxis(side=1,at=c(1e-11, 1e-10,1e-9,1e-8,1e-7,1e-6,1e-5,1e-4,1e-3,1e-2,0.1,1,10,100,1e3), cex.axis=1.4)
#for x linear
eaxis(side=1,at=seq(0, xlim[2], xtick), cex.axis=1.4)
minor.tick(nx=10)

lapply(dirs, function(x) {print(x);
 subdirs <- list.dirs(path=x, recursive=F)
 subdirs.ce <- subdirs[grep("ce", subdirs, ignore.case=T)]
 subdirs.tpv <- subdirs[grep("tpv", subdirs, ignore.case=T)]
 a <- read.table(paste(subdirs.ce,"/outputChargeDensityCE.txt",sep=""),header=T,stringsAsFactors=F)
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
 #tryCatch({
	  expfit <- nlrob(ChargeDensityCE~ exp(B)*Voc+exp(C)*(exp(exp(D)*Voc)-1), start=list(B=coef(expfit)["B"],C=coef(expfit)["C"],D=coef(expfit)["D"]), data=a)
 #}, error=function(e) {cat("FAILED robust FIT ", e$message, "\n")});

fulloutput <- read.table(paste(subdirs.tpv,"/output-robustmonoexp.txt",sep=""), header=TRUE);
n<-tail(grep("file",fulloutput[,1]),n=1)
tpv <- read.table(paste(subdirs.tpv,"/output-robustmonoexp.txt",sep=""), header=TRUE, skip=ifelse(length(n),n,0));
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
weights= (min(shown_T)/shown_T)^3
#just in case...
rm(powerlaw)

if(length(shown_T) < 4 || length(shown_charge) < 4){stop("you need wider plot limits!")}

j=1
while(!exists("powerlaw") && j < 1000){
	j <- j + 0.1
	start <- list(y0=log(5e-7*runif(1,1/j,j)), A=log(1e-28*runif(1,1/j,j)), alpha=-3.2*runif(1,1/j,j))
	tryCatch({
		powerlaw <- nlsLM(shown_T~exp(y0)+exp(A)*shown_charge^alpha, start=start, weights=weights)
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
	capture.output(summary(powerlaw), file=paste(x, "-tpvce-fit", sep=""),  append=TRUE);
}

i <<- i+1
})
legend(x="topright",inset=0.05,legend,pch=seq(21,25), pt.bg=mycolors, lwd=2, pt.lwd=1.5, pt.cex=2, col=change.lightness(mycolors,0.5),cex=1.5, title=title,bg="gray90", bty="n")
graphics.off()

#reset the plotting margins
par(op)

maxlength = max(sapply(output,length))
output = lapply(output, function(x){length(x)=maxlength; print(x)})
output = as.data.frame(output,check.names=FALSE)
write.table(output, file=paste(filename,"-TPVCEs.csv",sep=""), row.names=FALSE, na="", sep=",")

if(output_pdf){
	pdf(paste(filename,"-TPVCEs_nogeom.pdf",sep=""), width=image_smallpdf_width, height=image_smallpdf_height, pointsize=7)
}else{
	jpeg(quality=98, paste(filename,"-TPVCEs_nogeom.jpg",sep=""), width=image_width, height=image_height)
}

op <- par(mar = c(5,6,1,1) + 0.1) ## default is c(5,4,4,2) + 0.1
plot(1,xlim=xlim_nogeom,ylim=ylim_nogeom,xlab="", ylab="",log="xy", yaxt="n", xaxt="n")
#line is for introducing more space between label and axis
title(ylab = "Charge carrier lifetime (s)", cex.lab = 1.7, line = 4)
title(xlab = bquote("Charge per area (C/cm"^"2"*")"), cex.lab = 1.7, line = 4)

eaxis(side=2,at=c(1e-11, 1e-10,1e-9,1e-8,1e-7,1e-6,1e-5,1e-4,1e-3,1e-2,0.1,1,10,100,1e3), cex.axis=1.4)
xtick = 10^(floor(log10(xlim_nogeom[2])))
eaxis(side=1,at=c(1e-11, 1e-10,1e-9,1e-8,1e-7,1e-6,1e-5,1e-4,1e-3,1e-2,0.1,1,10,100,1e3), cex.axis=1.4)
i <- 0
lapply(dirs, function(x) {print(x);
 subdirs <- list.dirs(path=x, recursive=F)
 subdirs.ce <- subdirs[grep("ce", subdirs, ignore.case=T)]
 subdirs.tpv <- subdirs[grep("tpv", subdirs, ignore.case=T)]
 a <- read.table(paste(subdirs.ce,"/outputChargeDensityCE.txt",sep=""),header=T,stringsAsFactors=F)
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
 #tryCatch({
	  expfit <- nlrob(ChargeDensityCE~ exp(B)*Voc+exp(C)*(exp(exp(D)*Voc)-1), start=list(B=coef(expfit)["B"],C=coef(expfit)["C"],D=coef(expfit)["D"]), data=a)
 #}, error=function(e) {cat("FAILED robust FIT ", e$message, "\n")});

fulloutput <- read.table(paste(subdirs.tpv,"/output-robustmonoexp.txt",sep=""), header=TRUE);
n<-tail(grep("file",fulloutput[,1]),n=1)
tpv <- read.table(paste(subdirs.tpv,"/output-robustmonoexp.txt",sep=""), header=TRUE, skip=ifelse(length(n),n,0));

charge_nogeom <- exp(coef(expfit)[2])*(exp(exp(coef(expfit)[3])*tpv$Voc)-1)
points(charge_nogeom, tpv$T, bg=add.alpha(mycolors[i+1],0.5), col=change.lightness(mycolors[i+1],0.5), cex=1.5, pch=21+(i%%5));
index_shown_charge_nogeom = which(charge_nogeom >= xlim_nogeom[1] & charge_nogeom <= xlim_nogeom[2])
shown_charge_nogeom = charge_nogeom[index_shown_charge_nogeom]
shown_T_nogeom = tpv$T[index_shown_charge_nogeom]
weights_nogeom = (1/(shown_T_nogeom/min(shown_T_nogeom)))^3

#just in case...
rm(powerlaw_nogeom)

if(length(shown_T_nogeom) < 4 || length(shown_charge_nogeom) < 4){graphics.off(); stop("you need wider plot limits!")}
j=1
while(!exists("powerlaw_nogeom") && j < 1000){
	j <- j + 0.1
	start_nogeom <- list(y0=log(5e-7*runif(1,1/j,j)), A=log(1e-28*runif(1,1/j,j)), alpha=-3.2*runif(1,1/j,j))
	tryCatch({
		powerlaw_nogeom <- nlsLM(shown_T_nogeom~exp(y0)+exp(A)*shown_charge_nogeom^alpha, start=start_nogeom, weights=weights_nogeom)
		#check convergence and sum the p-values
	}, error=function(e) {cat("FAILED POWERLAW FIT ", e$message, "\n");
	})
		#summary fails if the fit was done on no data, with some chol2inv error
		if(exists("powerlaw_nogeom")){
			print("Checking powerlaw result")
			if(!summary(powerlaw_nogeom)$convInfo$isConv || sum(coef(summary(powerlaw_nogeom))[,"Pr(>|t|)"]) > 1){
				rm(powerlaw_nogeom)
			}
		}
}
if(exists("powerlaw_nogeom")){
	lines(shown_charge_nogeom, predict(powerlaw_nogeom, shown_charge_nogeom), lwd=2, col=change.lightness(mycolors[i+1],0.5))
	capture.output(summary(powerlaw_nogeom), file=paste(x, "-tpvce-nogeom-fit", sep=""),  append=TRUE);
}

i <<- i+1
})
legend(x="topright",inset=-0.01,legend,pch=seq(21,25), pt.bg=mycolors, lwd=2, pt.lwd=1.5, pt.cex=2, col=change.lightness(mycolors,0.5),cex=1.5, title=title,bg="gray90", bty="n")
graphics.off()

#reset the plotting margins
par(op)

