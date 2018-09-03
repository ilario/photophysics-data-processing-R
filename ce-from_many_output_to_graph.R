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
colors=gsub(".*-col_","",dirs)
# if the color is not set, use the default one
if(!length(colors[1])){colors=colorRampPalette(c("red","orange","springgreen","royalblue"))(max(length(dirs),3))}

data <- lapply(dirs, function(x) {print(x);
 a <- read.table(paste(x,"/ce/outputChargeDensityCE.txt",sep=""),header=T,stringsAsFactors=F)
 output[[paste("Voc",sub("nm","",sub("_.*","",sub("^0","",x))),sep="")]] <<- a$Voc
 a <- a[with(a, order(a$Voc)),]
 lin <- lm(ChargeDensityCE ~ 0 + Voc, data=a)
 tryCatch({
 exp <- nlsLM(ChargeDensityCE~ exp(B)*Voc+exp(C)*(exp(exp(D)*Voc)-1), start=list(B=log(max(a$ChargeDensityCE)/max(a$Voc)),C=log(1e-10),D=2), data=a)
 }, error=function(e) print("Failed fit"))
 tryCatch({
 exp <- nlsLM(ChargeDensityCE~ exp(B)*Voc+exp(C)*(exp(exp(D)*Voc)-1), start=list(B=log(coef(lin)[[1]]),C=log(1e-10),D=2), data=a)
 }, error=function(e) print("Failed fit"))
 tryCatch({
 exp <- nlrob(ChargeDensityCE~ exp(B)*Voc+exp(C)*(exp(exp(D)*Voc)-1), start=list(B=coef(exp)[[1]],C=coef(exp)[[2]],D=coef(exp)[[3]]), data=a)
 }, error=function(e) print("Failed robust fit"))
 
 print(exp)
 g <- predict(exp,a$Voc)
 a$g <- g
 output[[sub("_.*","",sub("^0","",x))]] <<- signif(g,5)
 onlyexp <- exp(coef(exp)[2])*(exp(exp(coef(exp)[3])*a$Voc)-1)
 a$onlyexp <- onlyexp
 a})
names(data) <- dirs

maxlength = max(sapply(output,length))
output = lapply(output, function(x){length(x)=maxlength; print(x)})
output = as.data.frame(output,check.names=FALSE)
write.table(output, file=paste(filename,"-CEs.csv",sep=""), row.names=FALSE, na="", sep=",")

#jpeg(quality=98, paste(filename,"-CEs-onlyexp.jpg",sep=""), width=image_width, height=image_height)
#op <- par(mar = c(5,8,4,2) + 0.1) ## default is c(5,4,4,2) + 0.1 
#plot(NULL,xlim=xlim,ylim=c(1e-10,2.5e-8),cex.main=1.5,xlab="",ylab="",  cex.lab=2, cex.axis=1.5);
#title(ylab = bquote("Charge density (C/cm"^"2"*")"), cex.lab = 2, line = 5)
#title(xlab = "Voltage (V)", cex.lab = 2, line = 3)
#
##eaxis(side=2,at=c(1e-10,1e-9,1e-8,1e-7,1e-6,1e-5,1e-4,1e-3,1e-2,0.1,1,10,100,1e3), cex.axis=1.5)
##minor.tick(nx=10)
#lapply(dirs, function(x) {print(x);
## points(data[[x]]$Voc, data[[x]]$ChargeDensityCE, lwd=0.2, bg=add.alpha(colors[i+1],0.5), pch=21+(i%%5), cex=2)
##just the exponential part
#lines(data[[x]]$Voc, data[[x]]$onlyexp, col=colors[i+1],lwd=3)
#lines(data[[x]]$Voc, data[[x]]$onlylin, col=colors[i+1],lwd=3)
#lines(data[[x]]$Voc, data[[x]]$explinsum, col=colors[i+1],lwd=3)
## mtext(bquote(.(gsub("-outputChargeDensityCE.txt","",x))~": n" == .(signif(exp$coefficients["A"],3)) + 
##	      .(signif(exp$coefficients["C"],3)) ~ "e" ^ {.(signif(exp$coefficients["D"],3))~V}),side=3,line=-(i*2+4),cex=1.5,col=colors[i+1])
## mtext(bquote(.(x)~": n" == .(signif(exp$coefficients["A"],3)) + 
##	      .(signif(exp$coefficients["C"],3)) ~ "e" ^ {.(signif(exp$coefficients["D"],3))~V}),side=3,line=-(i*2+4),cex=1.5,col=colors[i+1])
# i <<- i+1
#})
#legend(x="topleft",inset=0.05,legend, pch=seq(21,25), pt.bg=colors, col=colors, pt.cex=2, cex=2, pt.lwd=2, lwd=4, title=title, bg="gray90", bty="n")
#graphics.off()
##reset the plotting margins
#par(op)

i<-0
jpeg(quality=98, paste(filename,"-CEs.jpg",sep=""), width=image_width, height=image_height)
op <- par(mar = c(5,8,4,2) + 0.1) ## default is c(5,4,4,2) + 0.1 
plot(NULL,xlim=xlim,ylim=ylim,cex.main=1.5,xlab="", ylab="", cex.lab=2, cex.axis=1.5, yaxt="n");
title(ylab = bquote("Charge density (C/cm"^"2"*")"), cex.lab = 2, line = 5)
title(xlab = "Voltage (V)", cex.lab = 2, line = 3)

eaxis(side=2, cex.axis=1.5)
minor.tick(nx=10, ny=10)
lapply(dirs, function(x) {print(x);
 points(data[[x]]$Voc, data[[x]]$ChargeDensityCE, lwd=0.2, bg=add.alpha(colors[i+1],0.5), pch=21+(i%%5), cex=2)
 lines(data[[x]]$Voc, data[[x]]$g, col=change.lightness(colors[i+1],0.4),lwd=3)
 lines(data[[x]]$Voc, data[[x]]$onlyexp, col=colors[i+1],lwd=3)
# mtext(bquote(.(gsub("-outputChargeDensityCE.txt","",x))~": n" == .(signif(exp$coefficients["A"],3)) + 
#	      .(signif(exp$coefficients["C"],3)) ~ "e" ^ {.(signif(exp$coefficients["D"],3))~V}),side=3,line=-(i*2+4),cex=1.5,col=colors[i+1])
# mtext(bquote(.(x)~": n" == .(signif(exp$coefficients["A"],3)) + 
#	      .(signif(exp$coefficients["C"],3)) ~ "e" ^ {.(signif(exp$coefficients["D"],3))~V}),side=3,line=-(i*2+4),cex=1.5,col=colors[i+1])
 i <<- i+1
})
legend(x="topleft",inset=0.05,legend, pch=seq(21,25), pt.bg=colors, col=colors, pt.cex=2, cex=2, pt.lwd=2, lwd=4, title=title, bg="gray90",bty="n")
graphics.off()
#reset the plotting margins
par(op)
