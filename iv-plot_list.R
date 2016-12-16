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

name=""
title=gsub("_"," ",tail(unlist(strsplit(name,"-")),1))
filename=gsub(",","",gsub(":","",name))

fileslist = c("", "", "")
legendlist = c("", "", "")
colors = c("red","green", "blue")

i = 0
jpeg(quality=98, paste(filename,"-IVs.jpg",sep=""), width=640, height=480);
par(mar=c(5.1,7,2,2.1))
plot(NULL,xlim=c(-0.1,1),ylim=c(-16,5),cex.main=1.5,xlab="Voltage (V)",ylab=bquote("Current Density (mA/cm"^"2"*")"), cex.lab=1.5,cex.axis=1.2, yaxt="n", xaxt="n")#, main=name);
eaxis(side=2, cex.axis=1.2)
eaxis(side=1, cex.axis=1.2)
minor.tick(nx=10, ny=10)

abline(h=0, col="gray50");abline(v=0, col="gray50")
lapply(fileslist, function(x){print(x); fwd=paste(x, "-forward.txt", sep=""); rev=paste( x, "-reverse.txt", sep="");
#fwdV=mydata[[fwd]]$Voltage_V
#fwdJ=mydata[[fwd]]$Current_mA/0.09
revV=mydata[[rev]]$Voltage_V
revJ=mydata[[rev]]$Current_mA/0.09
##lines(fwdV, fwdJ, lwd=3, col=colors[i+1])
lines(revV, revJ, lwd=3, col=colors[i+1])
#fwdV=fwdV[c(T,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F)]
#fwdJ=fwdJ[c(T,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F)]
revJ=revJ[revV*10 == floor(revV*10)]
revV=revV[revV*10 == floor(revV*10)]
##points(fwdV, fwdJ, lwd=1, bg=colors[i+1], cex=2, pch=21+i)
points(revV, revJ, bg=colors[i+1], cex=2, pch=21+i)
i <<- i+1
})
legend(x="topleft",inset=c(0.15,0.05),legendlist, lty=c(1,1,1,1), pch=seq(21,25), lwd=4, pt.cex=2, pt.lwd=2, pt.bg=colors, cex=1.5, col=colors, title=title,bg="gray90")
graphics.off()

