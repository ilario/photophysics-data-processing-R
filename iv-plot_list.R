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

library(RColorBrewer)
library(sfsmisc)
library(Hmisc)

#name=""
#title=""
filename=gsub(",","",gsub(":","",name))

change.lightness <- function(col, lightness=1){
	if(missing(col))
		stop("Please provide a vector of colours.")
	apply(sapply(col, col2rgb, alpha=TRUE)/255, 2,
		function(x)
			rgb(x[1]*lightness, x[2]*lightness, x[3]*lightness, alpha=x[4]))
}

#fileslist=list.files(pattern="*-forward.*.txt|*-reverse.*.txt")
fileslist=names(mydata)
reverselist=fileslist[!grepl("forward", fileslist) & !grepl("dark", fileslist, ignore.case=T)]
cellslist=sub("-reverse","",sub(".txt","",reverselist))
cellslist=cellslist[!duplicated(cellslist)]
legendlist=sub("_.*","",sub("^0","",cellslist))
mycolors=gsub(".*-col_","",cellslist[grepl("-col_", cellslist)])

if(!length(mycolors)){
mycolors<-brewer.pal(8,"Dark2")
}

#jpeg(quality=98, paste(filename,"-revIVs.jpg",sep=""), width=image_width, height=image_height);
pdf(paste(filename,"-revIVs.pdf",sep=""), width=image_smallpdf_width, height=image_smallpdf_height, pointsize=7);
par(mar=c(5.1,6,1,1))
plot(NULL,xlim=lim.IV.voltage,ylim=lim.IV.current,xlab="Voltage (V)",ylab=bquote("Current Density (mA/cm"^"2"*")"), cex.lab=1.7, yaxt="n", xaxt="n");
eaxis(side=2, cex.axis=1.4)
eaxis(side=1, cex.axis=1.4)
minor.tick(nx=10, ny=10)

abline(h=0, col="gray50");abline(v=0, col="gray50")

i = 1
lapply(reverselist, function(x){ 
V=mydata[[x]]$Voltage_V
J=mydata[[x]]$Current_mA/0.09
lines(V, J, lwd=2, col=mycolors[i])
Jmarkers=J[V*10 == floor(V*10)]
Vmarkers=V[V*10 == floor(V*10)]
points(Vmarkers, Jmarkers, col=change.lightness(mycolors[i],0.5), bg=mycolors[i], cex=1.5, pch=20+(i%%5))
i <<- i+1
})

legend(x="topleft",inset=c(0.15,0.2),legendlist, lty=c(rep(1,length(fileslist))), lwd=2, pt.cex=2, pt.lwd=1.5, pt.bg=mycolors, cex=1.5, col=change.lightness(mycolors,0.5), bty="n", pch=seq(21,25))
graphics.off()

