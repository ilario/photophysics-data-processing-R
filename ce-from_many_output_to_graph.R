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
library(robustbase)
library(sfsmisc)
library(Hmisc)
require(minpack.lm)

ylim=lim.CE.charge
xlim=lim.CE.voltage

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

output=list()

i <- 0
dirs <- list.dirs(recursive=FALSE)
dirs <- sub("./","",dirs)
legend=sub("_.*","",sub("^0","",dirs))

# try to obtain the color from the file name
mycolors=gsub(".*-col_","",dirs[grepl("-col_", dirs)])
# if the color is not set, use the default one
if(!length(mycolors)){mycolors=brewer.pal(8,"Dark2")}

data <- lapply(dirs, function(x) {print(x);
 subdirs <- list.dirs(path=x, recursive=F)
 subdirs.ce <- subdirs[grep("ce", subdirs, ignore.case=T)]
 a <- read.table(paste(subdirs.ce,"/outputChargeDensityCE.txt",sep=""),header=T,stringsAsFactors=F)
 output[[paste("Voc",sub("nm","",sub("_.*","",sub("^0","",x))),sep="")]] <<- a$Voc
 a <- a[with(a, order(a$Voc)),]
 lin <- lm(ChargeDensityCE ~ 0 + Voc, data=a)
 tryCatch({
 expfit <- nlsLM(ChargeDensityCE~ exp(B)*Voc+exp(C)*(exp(exp(D)*Voc)-1), start=list(B=log(max(a$ChargeDensityCE)/max(a$Voc)),C=log(1e-10),D=2), data=a)
 }, error=function(e) cat("Failed fit ", e$message, "\n"))
 tryCatch({
 expfit <- nlsLM(ChargeDensityCE~ exp(B)*Voc+exp(C)*(exp(exp(D)*Voc)-1), start=list(B=log(coef(lin)[[1]]),C=log(1e-10),D=2), data=a)
 }, error=function(e) cat("Failed fit ", e$message, "\n"))
 #tryCatch({
 expfit <- nlrob(ChargeDensityCE~ exp(B)*Voc+exp(C)*(exp(exp(D)*Voc)-1), start=list(B=coef(expfit)[[1]],C=coef(expfit)[[2]],D=coef(expfit)[[3]]), data=a)
 #}, error=function(e) cat("Failed robust fit ", e$message, "\n"))

if(exists("expfit")){ 
 print(expfit)
 g <- predict(expfit,a$Voc)
 onlyexp <- exp(coef(expfit)[2])*(exp(exp(coef(expfit)[3])*a$Voc)-1)
} else {
 print("using linear fitting")
 g <- predict(lin,data.frame(Voc=a$Voc))
 onlyexp <- rep(0, length(a$Voc))
}
 a$g <- g
 output[[sub("_.*","",sub("^0","",x))]] <<- signif(g,5)
 a$onlyexp <- onlyexp
 a})
names(data) <- dirs

maxlength = max(sapply(output,length))
output = lapply(output, function(x){length(x)=maxlength; print(x)})
output = as.data.frame(output,check.names=FALSE)
write.table(output, file=paste(filename,"-CEs.csv",sep=""), row.names=FALSE, na="", sep=",")

i<-0
if(output_pdf){
	pdf(paste(filename,"-CEs.pdf",sep=""), width=image_smallpdf_width, height=image_smallpdf_height, pointsize=7)
}else{
	png(paste(filename,"-CEs.png",sep=""), width=image_width, height=image_height)
}
op <- par(mar = c(5,7.5,1,1) + 0.1) ## default is c(5,4,4,2) + 0.1 
plot(NULL,xlim=xlim,ylim=ylim,xlab="", ylab="", cex.axis=1.4, yaxt="n");
title(ylab = bquote("Charge per area (C/cm"^"2"*")"), cex.lab = 1.7, line = 5.5)
title(xlab = "Light bias (V)", cex.lab = 1.7, line = 3)

eaxis(side=2, cex.axis=1.4)
minor.tick(nx=10, ny=10)
lapply(dirs, function(x) {print(x);
 points(data[[x]]$Voc, data[[x]]$ChargeDensityCE, col=change.lightness(mycolors[i+1],0.5), bg=add.alpha(mycolors[i+1],0.5), pch=21+(i%%5), cex=1.5)
 lines(data[[x]]$Voc, data[[x]]$g, col=change.lightness(mycolors[i+1],0.5),lwd=2)
 lines(data[[x]]$Voc, data[[x]]$onlyexp, col=mycolors[i+1],lwd=2)
 i <<- i+1
})
legend(x="topleft",inset=0.05,legend, pch=seq(21,25), pt.bg=mycolors, col=change.lightness(mycolors,0.5), pt.cex=2, cex=1.5, pt.lwd=1.5, lwd=3, title=title, bty="n")
graphics.off()
#reset the plotting margins
par(op)
