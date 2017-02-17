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

trashfornullmessages <- lapply(names(mydata), function(x){print(x);cell = substr(x,1,nchar(x)-12);
      if(!file.exists(paste(cell,".png",sep=""))){
       fwd=gsub("reverse","forward",x);rev=gsub("forward","reverse",x);dark=gsub("-bis|-tris|-quater|-penta|-exa|-epta|-octa|-nona|-deca|-undeca|-dodeca","",gsub("forward","dark-forward",x)); title=paste(cell, "(",results[[x]][5],")"); 
       png(paste(cell,".png",sep=""), width=640, height=640);
#      pdf(paste(cell,".pdf",sep=""),title=title, width=8, height=8);
       plot(NULL,xlim=c(min(min(mydata[[fwd]]$Voltage_V),min(mydata[[rev]]$Voltage_V)),max(max(mydata[[fwd]]$Voltage_V),max(mydata[[rev]]$Voltage_V))),ylim=c(-2.1,2),main=title,cex.main=1.5,xlab="Voltage (V)",ylab="Current (mA)");
       lines(mydata[[fwd]]$Voltage_V, mydata[[fwd]]$Current_mA, col="red", lwd=2); 
       tryCatch({arrows(mydata[[fwd]]$Voltage_V[1]-0.01, mydata[[fwd]]$Current_mA[1], mydata[[fwd]]$Voltage_V[1], mydata[[fwd]]$Current_mA[1], col="red", length=0.1); arrows(tail(mydata[[fwd]]$Voltage_V, n=1), tail(mydata[[fwd]]$Current_mA, n=1),tail(mydata[[fwd]]$Voltage_V, n=1)+0.01, tail(mydata[[fwd]]$Current_mA, n=1), col="red", length=0.1)}, error=function(e) {print("missing fwd");}); 
       tryCatch({arrows(mydata[[rev]]$Voltage_V[1]+0.01, mydata[[rev]]$Current_mA[1], mydata[[rev]]$Voltage_V[1], mydata[[rev]]$Current_mA[1], col="red", length=0.1); arrows(tail(mydata[[rev]]$Voltage_V, n=1), tail(mydata[[rev]]$Current_mA, n=1),tail(mydata[[rev]]$Voltage_V, n=1)-0.01, tail(mydata[[rev]]$Current_mA, n=1), col="red", length=0.1)}, error=function(e) {print("missing rev");}); 
       lines(mydata[[rev]]$Voltage_V, mydata[[rev]]$Current_mA, col="red", lwd=2, lty=2); 
       lines(mydata[[dark]]$Voltage_V, mydata[[dark]]$Current_mA, col="blue", lwd=2); 
       abline(h=0);abline(v=0);mtext(paste(" Forward:", "Jsc",results[[fwd]][1],"Voc",results[[fwd]][2],"FF",results[[fwd]][3],"eff",results[[fwd]][4]), side=3, line=-4, adj=0, col="red", cex=1.5);mtext(paste(" Reverse:", "Jsc",results[[rev]][1],"Voc",results[[rev]][2],"FF",results[[rev]][3],"eff",results[[rev]][4]), side=3, line=-6, adj=0, col="red", cex=1.5);graphics.off()}})

