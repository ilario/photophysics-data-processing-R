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

directory=tail(strsplit(getwd(), "/")[[1]], n=1)
name=directory

files <- list.files(path=".", pattern="^ig.*\\.txt$");

illumination <- function(names){
  values <- lapply(names, function(name){
    splittedname <- strsplit(name, "-")
    value="1 sun"
    if(grepl("dark",name)){value="dark"}
    if(grepl("sun",name)){value=paste(as.numeric(gsub("sun","",unlist(splittedname)[grepl("sun",unlist(splittedname))])), " sun")}
    return(value)})
  return(unlist(values))
}



suns <- illumination(files)
files <- files[order(suns)]
uniquesortsuns <- unique(sort(suns))

filesNoDark = files[!grepl("dark", files)];
filesDark = files[grepl("dark", files)];
colors=colorRampPalette(c("blue","red"))(max(length(uniquesortsuns[uniquesortsuns != "dark"]),3))
print(uniquesortsuns[uniquesortsuns != "dark"])

png(paste(name,"-suns.png",sep=""), width=1000, height=800);
plot(NULL,xlim=c(-0.2,1.1),ylim=c(-2.1,1.2),cex.main=1.5,xlab="Voltage (V)",ylab="Current (mA)", cex.lab=1.5)#, main=paste(name, "at various Light Intensities"));
lines(mydata[[filesDark[1]]]$Voltage_V, mydata[[filesDark[1]]]$Current_mA, lwd=1, col="black")
i<-1
lapply(filesNoDark[grepl("reverse",files)], function(x){print(x); print(colors[i]); 
  lines(mydata[[x]]$Voltage_V, mydata[[x]]$Current_mA, lwd=1, col=colors[i]);
  i<<-i+1
})
i<-1
lapply(filesNoDark[grepl("forward",files)], function(x){print(x); print(colors[i]); 
  lines(mydata[[x]]$Voltage_V, mydata[[x]]$Current_mA, lwd=1, lty=2, col=colors[i])
  i<<-i+1
})

abline(h=0);abline(v=0)
legend(x="topleft",inset=c(0.2,0.05),legend=c("dark", uniquesortsuns[uniquesortsuns != "dark"]), col=c("#000000", colors), cex=1.5, title="Illumination Intensity:", lwd=4, bty="n")
graphics.off()

