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

if(!exists("name")){
  name=basename(getwd())
}
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
forwardlist=fileslist[grepl("forward", fileslist) & !grepl("dark", fileslist, ignore.case=T)]
reverselist=fileslist[!grepl("forward", fileslist) & !grepl("dark", fileslist, ignore.case=T)]
forwarddarklist=fileslist[grepl("forward", fileslist) & grepl("dark", fileslist, ignore.case=T)]
reversedarklist=fileslist[!grepl("forward", fileslist) & grepl("dark", fileslist, ignore.case=T)]
cellslist=sub("-forward|-reverse","",sub(".txt","",fileslist))
cellslist=cellslist[!duplicated(cellslist)]
#remove everything after the last underscore, then convert the remaining underscores to spaces
legendlist=gsub("_"," ",sub("_((?!_).)+$","",cellslist, perl=TRUE))
mycolors=gsub(".*-col_","",cellslist[grepl("-col_", cellslist)])

if(!length(mycolors)){
  mycolors<-brewer.pal(8,"Dark2")
}

if(output_pdf){
  pdf(paste(filename,"-IVs.pdf",sep=""), width=image_midpdf_width, height=image_midpdf_height, pointsize=7)
}else{
  png(paste(filename,"-IVs.png",sep=""), width=image_width, height=image_height);
}
par(mar=c(5.1,6,1,1))
plot(NULL,xlim=lim.IV.voltage,ylim=lim.IV.current,xlab="Voltage (V)",ylab=bquote("Current Density (mA/cm"^"2"*")"), cex.lab=1.7, yaxt="n", xaxt="n");
eaxis(side=2, cex.axis=1.4)
eaxis(side=1, cex.axis=1.4)
minor.tick(nx=10, ny=10)

abline(h=0, col="gray50");abline(v=0, col="gray50")
#i = 1
#lapply(forwarddarklist, function(x){ 
#  V=mydata[[x]]$Voltage_V
#  J=mydata[[x]]$Current_mA/0.09
#  lines(V, J, lwd=2, col=mycolors[i], lty=2)
#  indexes=seq(i*5,300+i*5,20)
#  Jmarkers=J[indexes]
#  Vmarkers=V[indexes]
#  #Jmarkers=J[V*10 == floor(V*10)]
#  #Vmarkers=V[V*10 == floor(V*10)]
#  points(Vmarkers, Jmarkers, col=mycolors[i], cex=1.5, pch=20+(i%%5), lwd=2, bg="white")
#  i <<- i+1
#})
#i = 1
#lapply(reversedarklist, function(x){ 
#  V=mydata[[x]]$Voltage_V
#  J=mydata[[x]]$Current_mA/0.09
#  lines(V, J, lwd=2, col=mycolors[i])
#  indexes=seq(i*5,300+i*5,20)
#  Jmarkers=J[indexes]
#  Vmarkers=V[indexes]
#  #Jmarkers=J[V*10 == floor(V*10)]
#  #Vmarkers=V[V*10 == floor(V*10)]
#  points(Vmarkers, Jmarkers, col=change.lightness(mycolors[i],0.5), bg=mycolors[i], cex=1.5, pch=20+(i%%5))
#  i <<- i+1
#})
i = 1
lapply(forwardlist, function(x){ 
  V=mydata[[x]]$Voltage_V
  J=mydata[[x]]$Current_mA/0.09
  lines(V, J, lwd=2, col=mycolors[i], lty=2)
  i <<- i+1
})
i = 1
lapply(reverselist, function(x){ 
  V=mydata[[x]]$Voltage_V
  J=mydata[[x]]$Current_mA/0.09
  lines(V, J, lwd=2, col=mycolors[i])
  i <<- i+1
})
i = 1
lapply(forwardlist, function(x){ 
  V=mydata[[x]]$Voltage_V
  J=mydata[[x]]$Current_mA/0.09
  indexes=seq(i*2+10,300+i*2+10,20)
  Jmarkers=J[indexes]
  Vmarkers=V[indexes]
  points(Vmarkers, Jmarkers, col=change.lightness(mycolors[i],0.7), cex=1.5, pch=20+(i%%5), lwd=2, bg="white")
  i <<- i+1
})
i = 1
lapply(reverselist, function(x){ 
  V=rev(mydata[[x]]$Voltage_V)#reversed for having a meaningful indexes numbering, the same direction as in forward
  J=rev(mydata[[x]]$Current_mA/0.09)
  indexes=seq(i*2,300+i*2,20)
  Jmarkers=J[indexes]
  Vmarkers=V[indexes]
  points(Vmarkers, Jmarkers, col=change.lightness(mycolors[i],0.5), bg=mycolors[i], cex=1.5, pch=20+(i%%5))
  i <<- i+1
})

legend(x="topleft",inset=c(0.15,0.2),legendlist, lty=c(rep(1,length(fileslist))), lwd=3, pt.cex=2, pt.lwd=1.5, pt.bg=mycolors, cex=1.5, col=change.lightness(mycolors,0.5),bg="white",box.col="grey90", pch=seq(21,25),title=title)
graphics.off()

