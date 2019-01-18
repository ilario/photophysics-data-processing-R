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

tpvOnlyDeltaV <- function(tpvdir="tpv")
{
  print("TPV ONLY DELTA V: CALCULATING")
  
  noiseTime=5e-8
  
  files <- list.files(path=tpvdir, pattern="^TPV.*\\.txt.table$");
  mydata <- lapply(file.path(tpvdir,files), read.table, header=FALSE, col.names=c("time","voltage"));
  files <- sub(".txt.table","",files);
  names(mydata) <- files;
  ## output for importing
  write.table(t(c("file","Voc","deltaV")), file=file.path(tpvdir,"outputDeltaV.txt"), append=FALSE, col.names=F, row.names=F);
  write.table(t(c("file","Voc","deltaV")), file=file.path(tpvdir,"outputDeltaVloess.txt"), append=FALSE, col.names=F, row.names=F);
  write.table(t(c("file","Voc","deltaV")), file=file.path(tpvdir,"outputDeltaVfirstPoints.txt"), append=FALSE, col.names=F, row.names=F);
  
  trashfornullmessages <- lapply(files, function(x) {
    message(x);	
    peaktime <- mydata[[x]]$time[which.max(mydata[[x]]$voltage)];
    voltage2 <- subset(mydata[[x]], time >= peaktime, select=voltage); 
    time2 <- subset(mydata[[x]], time >= peaktime, select=time);
    temp <- data.frame(time2, voltage2);
    endtime = tail(temp$time, n=1)
    startingvoltage <- mean(mydata[[x]]$voltage[0:200]);
    maxVoltageIndex <- which.max(mydata[[x]]$voltage)
    deltavoltage <- mydata[[x]]$voltage[maxVoltageIndex] - startingvoltage;
    outputDeltaV <- t(c(x, startingvoltage, deltavoltage));
    write.table(outputDeltaV, file=file.path(tpvdir,"outputDeltaV.txt"), append=TRUE, col.names=F, row.names=F, quote=F);
    
    tempsubset2 <- subset(mydata[[x]], time > endtime/120, select=c(time,voltage));
    lo2 <- loess(tempsubset2$voltage~tempsubset2$time, span=0.02);
    deltaVloess <- max(predict(lo2))-startingvoltage
    outputDeltaVloess <- t(c(x, startingvoltage, deltaVloess));
    write.table(outputDeltaVloess, file=file.path(tpvdir,"outputDeltaVloess.txt"), append=TRUE, col.names=F, row.names=F, quote=F);
    
    afterPeakNoiseTimePoints = time2 < peaktime + noiseTime
    afterPeakNoiseVoltagePoints = voltage2[afterPeakNoiseTimePoints]
    print(paste("DeltaVfirstPoints averaging on", length(afterPeakNoiseVoltagePoints), "points"))
    deltaVfirstPoints <- mean(afterPeakNoiseVoltagePoints) - startingvoltage
    outputDeltaVfirstPoints <- t(c(x, startingvoltage, deltaVfirstPoints));
    write.table(outputDeltaVfirstPoints, file=file.path(tpvdir,"outputDeltaVfirstPoints.txt"), append=TRUE, col.names=F, row.names=F, quote=F);
  })}


