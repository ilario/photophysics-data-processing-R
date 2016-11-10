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
fileslist = c("", "", "", "")
legendlist = c("", "", "", "")
#library(RColorBrewer)
#colors=brewer.pal(max(length(fileslist),3),"Set1")
colors = c("red","green", "blue")

i = 1
png(paste(gsub(" ", "_", name),"-IVs.png",sep=""), width=640, height=640);
plot(NULL,xlim=c(-0.1,1.1),ylim=c(-2.1,0.5),cex.main=1.5,xlab="Voltage (V)",ylab="Current (mA)", cex.lab=1.5, main=name);
lapply(fileslist, function(x){print(x); fwd=paste(x, "-forward.txt", sep=""); rev=paste( x, "-reverse.txt", sep="");
lines(mydata[[fwd]]$Voltage_V, mydata[[fwd]]$Current_mA, lwd=3, col=colors[i])#brewer.pal(9,"Reds")[i])
lines(mydata[[rev]]$Voltage_V, mydata[[rev]]$Current_mA, lwd=3, lty=2, col=colors[i])#brewer.pal(9,"Reds")[i], )
i <<- i+1
})
abline(h=0);abline(v=0)
legend(x="topleft",inset=0.05,legendlist,lty=c(1,1,1,1), lwd=4,col=colors, cex=2)
graphics.off()

