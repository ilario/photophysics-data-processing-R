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

source("~/software/r/PvCurvesR/extractdata-curves-vi-separated_files.R")
source("~/software/r/PvCurvesR/iv-generate_mydata.R")

files <- list.files(path=".", pattern="\\.txt$");
files <- grep(name,files,value=T)
#files <- grep("sun",files,value=T)

png(paste(name,"-suns.png",sep=""), width=640, height=640);
plot(NULL,xlim=c(-0.1,1),ylim=c(-1.7,0.5),cex.main=1.5,xlab="Voltage (V)",ylab="Current (mA)", cex.lab=1.5, main=paste(name, "at various Light Intensities"));
lapply(files[grepl("forward",files)], function(x){print(x); 
lines(mydata[[x]]$Voltage_V, mydata[[x]]$Current_mA, lwd=2)
})
lapply(files[grepl("reverse",files)], function(x){print(x); 
lines(mydata[[x]]$Voltage_V, mydata[[x]]$Current_mA, lwd=2, lty=2);
})
abline(h=0);abline(v=0)
graphics.off()

