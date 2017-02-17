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
#files <- grep(name,files,value=T)
#files <- grep("sun",files,value=T)
#splittedfilename <- strsplit(files, "-")
#suns <- gsub("sun","",unique(unlist(splittedfilename)[grepl("sun",unlist(splittedfilename))]))

illumination <- function(names){
	values <- lapply(names, function(name){
				 splittedname <- strsplit(name, "-")
				 value=100
				 if(grepl("dark",name)){value=0}
				 if(grepl("sun",name)){value=100*as.numeric(gsub("sun","",unlist(splittedname)[grepl("sun",unlist(splittedname))]))}
				 return(value)})
	return(unlist(values))
}

colorsfwd=colorRampPalette(c("black","blue","green"))(max(length(files[grepl("forward",files)]),3))
colorsrev=colorRampPalette(c("black","blue","green"))(max(length(files[grepl("reverse",files)]),3))

suns <- illumination(files)
files <- files[order(suns)]
uniquesortsuns <- unique(sort(suns))

i<-1
png(paste(name,"-suns.png",sep=""), width=640, height=640);
plot(NULL,xlim=c(-0.2,1.2),ylim=c(-1.4,1),cex.main=1.5,xlab="Voltage (V)",ylab="Current (mA)", cex.lab=1.5)#, main=paste(name, "at various Light Intensities"));
lapply(files[grepl("forward",files)], function(x){print(x); 
lines(mydata[[x]]$Voltage_V, mydata[[x]]$Current_mA, lwd=1, col=colorsfwd[i])
i<<-i+1
})
i<-1
lapply(files[grepl("reverse",files)], function(x){print(x); 
lines(mydata[[x]]$Voltage_V, mydata[[x]]$Current_mA, lwd=1, lty=2, col=colorsrev[i]);
i<<-i+1
})
abline(h=0);abline(v=0)
legend(x="topleft",inset=c(0.2,0.05),legend=uniquesortsuns, col=colorRampPalette(c("black","blue","green"))(length(uniquesortsuns)), cex=1.5, title="Illumination Intensity:", lwd=4, bty="n")
graphics.off()

