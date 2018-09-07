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
forwardlist=fileslist[grepl("forward", fileslist) & !grepl("dark", fileslist)]
reverselist=fileslist[!grepl("forward", fileslist) & !grepl("dark", fileslist)]
forwarddarklist=fileslist[grepl("forward", fileslist) & grepl("dark", fileslist)]
reversedarklist=fileslist[!grepl("forward", fileslist) & grepl("dark", fileslist)]
cellslist=sub("-forward|-reverse","",sub(".txt","",fileslist))
cellslist=cellslist[!duplicated(cellslist)]
legendlist=sub("_.*","",sub("^0","",cellslist))
colors=gsub(".*-col_","",cellslist)
#colors=colorRampPalette(c("red","orange","springgreen","royalblue"))(max(length(fileslist),3))

jpeg(quality=98, paste(filename,"-IVs.jpg",sep=""), width=image_width, height=image_height);
par(mar=c(5.1,7,2,2.1))
plot(NULL,xlim=lim.IV.voltage,ylim=lim.IV.current,cex.main=1.5,xlab="Voltage (V)",ylab=bquote("Current Density (mA/cm"^"2"*")"), cex.lab=2,cex.axis=1.5, yaxt="n", xaxt="n");
eaxis(side=2, cex.axis=1.5)
eaxis(side=1, cex.axis=1.5)
minor.tick(nx=10, ny=10)

abline(h=0, col="gray50");abline(v=0, col="gray50")
i = 1
lapply(forwarddarklist, function(x){ 
V=mydata[[x]]$Voltage_V
J=mydata[[x]]$Current_mA/0.09
lines(V, J, lwd=4, col=colors[i])
Jmarkers=J[V*10 == floor(V*10)]
Vmarkers=V[V*10 == floor(V*10)]
points(Vmarkers, Jmarkers, col=change.lightness(colors[i],0.5), cex=2, pch=">")#21+(i%%5))
i <<- i+1
})
i = 1
lapply(reversedarklist, function(x){ 
V=mydata[[x]]$Voltage_V
J=mydata[[x]]$Current_mA/0.09
lines(V, J, lwd=4, col=colors[i])
Jmarkers=J[V*10 == floor(V*10)]
Vmarkers=V[V*10 == floor(V*10)]
points(Vmarkers, Jmarkers, col=change.lightness(colors[i],0.5), cex=2, pch="<")#21+(i%%5))
i <<- i+1
})
i = 1
lapply(forwardlist, function(x){ 
V=mydata[[x]]$Voltage_V
J=mydata[[x]]$Current_mA/0.09
lines(V, J, lwd=4, col=colors[i])
Jmarkers=J[V*10 == floor(V*10)]
Vmarkers=V[V*10 == floor(V*10)]
points(Vmarkers, Jmarkers, col=change.lightness(colors[i],0.5), cex=2, pch=">")#21+(i%%5))
i <<- i+1
})
i = 1
lapply(reverselist, function(x){ 
V=mydata[[x]]$Voltage_V
J=mydata[[x]]$Current_mA/0.09
lines(V, J, lwd=4, col=colors[i])
Jmarkers=J[V*10 == floor(V*10)]
Vmarkers=V[V*10 == floor(V*10)]
points(Vmarkers, Jmarkers, col=change.lightness(colors[i],0.5), cex=2, pch="<")#21+(i%%5))
i <<- i+1
})

legend(x="topleft",inset=c(0.15,0.2),legendlist, lty=c(rep(1,length(fileslist)),1,2), lwd=6, pt.cex=2, pt.lwd=2, pt.bg=colors, cex=2, col=colors, bty="n")
graphics.off()

