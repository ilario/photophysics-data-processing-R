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

name=""
title=gsub("-","\n\n",gsub("_"," ",name))
filename=gsub(",","",gsub(":","",name))

fileslist = c("", "", "", "")
legendlist = c("", "", "", "")
#library(RColorBrewer)
#colors=brewer.pal(max(length(fileslist),3),"Set1")
colors = c("red","green", "blue")

i = 1
png(paste(filename,"-IVs.png",sep=""), width=640, height=640);
par(mar=c(5.1,7,2,2.1))
plot(NULL,xlim=c(-0.1,1.1),ylim=c(-23,5),cex.main=1.5,xlab="Voltage (V)",ylab="Current (mA)", cex.lab=1.5,cex.axis=1.2, yaxt="n", xaxt="n")#, main=name);
eaxis(side=2, cex.axis=1.2)
eaxis(side=1, cex.axis=1.2)
minor.tick(nx=10, ny=10)

lapply(fileslist, function(x){print(x); fwd=paste(x, "-forward.txt", sep=""); rev=paste( x, "-reverse.txt", sep="");
fwdV=mydata[[fwd]]$Voltage_V
fwdJ=mydata[[fwd]]$Current_mA
revV=mydata[[rev]]$Voltage_V
revJ=mydata[[rev]]$Current_mA
lines(fwdV, fwdJ, lwd=3, col=colors[i])#brewer.pal(9,"Reds")[i])
lines(revV, revJ, lwd=3, lty=2, col=colors[i])#brewer.pal(9,"Reds")[i], )
fwdV=fwdV[c(T,F,F,F,F)]
fwdJ=fwdJ[c(T,F,F,F,F)]
revV=revV[c(T,F,F,F,F)]
revJ=revJ[c(T,F,F,F,F)]
points(fwdV, fwdJ, lwd=3, col=colors[i])
points(revV, revJ, lwd=3, col=colors[i])
i <<- i+1
})
abline(h=0);abline(v=0)
legend(x="topleft",inset=0.1,legendlist, lty=c(1,1,1,1), lwd=4,pch=seq(21,25), lwd=4, pt.cex=2, pt.lwd=2, pt.bg=colors, cex=1.5, col=colors, title=title,bg="gray90")
graphics.off()

