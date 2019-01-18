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

import.iv.separated <- function(dir = ".", pattern="", pattern.excl="\\.png$", list.excl="output.txt", header=TRUE, skip=22){
  files.all <- list.files(path=dir, pattern=pattern)
  #	files.all <- grep(pattern=grep_pattern, files.all)
  files.pattern.excl <- list.files(path=dir, pattern=pattern.excl)
  files.excl <- append(files.pattern.excl, list.excl)
  files <- files.all[!files.all %in% files.excl]
  mydata <- lapply(files, function(x){
    tryCatch({read.table(x, header=header, skip=skip)}, error=function(e) {write(paste("Import error with file", x), stderr())})})
  names(mydata) <- files
  mydata <- Filter(function(x) length(x) > 0, mydata)
  return(mydata)
}

#for importing all files in the current directory into a dataframe
#mydata <- import.iv.separated(pattern.excl="\\.pdf$")

#for analyzing a measurement
#extract.iv(mytable$Volts, mytable$Photocurrent.mA.)

#for analyzing a dataframe with data
#results <- lapply(names(mydata), function(x){extract.iv(mydata[[x]]$Voltage_V, mydata[[x]]$Current_mA, comment="")})
#names(results) <- names(mydata)

#for analyzing a dataframe with data and add comment based on filename
#results <- lapply(names(mydata), function(x){y <- factor(strsplit(x,"-")[[1]][3], c(30,31,32,33,34,11,20,22,38)); levels(y) <- list(`with QD spin2000rpm` = c(22,31), `with QD spin3000rpm` = c(30, 32), `with OMeTAD` = c(34,38), "with QD and MPA" = c(11,20)); extract.iv(mydata[[x]]$Voltage_V, mydata[[x]]$Current_mA, comment=toString(y))})
#names(results) <- names(mydata)

#for writing data to a file
#the header
#write(paste(format("directory",width=10), format("sample",width=20), format("reverse",width=8), format("Jsc",width=9), format("Voc",width=9), format("FF",width=6), "efficiency"), file="output.txt", append=TRUE);
#
#saving to a file the result from a SINGLE measurement
#write(extract.iv(ig1.5.1.r.bis$Volts, ig1.5.1.r.bis$Photocurrent.mA., formatted.output=TRUE, directory=20101020, sample="ig1-5-1-r-bis", reverse=1), file="output.txt", append=TRUE)

#for determining the date:
#date <- strsplit(tail(strsplit(getwd(), "/")[[1]], n=1), "-")[[1]][[1]]
#for determining if reverse:
#reverse <- as.integer(grepl("-r(-|$)", x));
#reverse <- as.integer(grepl("reverse", x));
#
#saving to a file the results from a dataframe of measurements
#lapply(names(mydata), function(x){write(extract.iv(mydata[[x]]$Voltage_V, mydata[[x]]$Current_mA, formatted.output=TRUE, directory=20101020, sample=x, reverse=as.integer(grepl("reverse", x))), file="output.txt", append=TRUE)})

extract.iv <- function(voltage, current, current.positive=FALSE, cell.surface=0.09, irradiance=100, formatted.output=FALSE, directory, sample, reverse, comment=""){
  efficiency <- extract.iv.power(voltage, current, current.positive, cell.surface, irradiance);
  jsc <- extract.iv.jsc(voltage, current, current.positive, cell.surface);
  voc <- extract.iv.voc(voltage, current);
  ff <- efficiency / (jsc * voc);
  if(formatted.output){
    result.formatted <- paste(as.character(directory), as.character(sample), reverse, as.character(format(jsc,digits=4)), as.character(format(voc,digits=4)), as.character(format(ff,digits=2)), as.character(format(efficiency,digits=4)), gsub(" ", "-", comment), sep="\t");
    return(result.formatted);
  }else{
    names <- c("Jsc", "Voc", "FF", "efficiency", "comment");
    result <- c(format(jsc,digits=4), format(voc,digits=4), format(ff,digits=2), format(efficiency,digits=4), comment);
    names(result) <- names;
    return(result);
  }
}


extract.iv.power <- function(voltage, current, current.positive=FALSE, cell.surface=0.09, irradiance=100){
  power <- voltage * current * (as.numeric(current.positive)*2-1);
  efficiency <- 100 * (max(power) / cell.surface) / irradiance;
  return(efficiency);
}

extract.iv.jsc <- function(voltage, current, current.positive=FALSE, cell.surface=0.09){
  pos.jsc <- which.min(abs(voltage));
  length <- length(voltage);
  range.jsc <- max(pos.jsc-3,1):min(pos.jsc+3,length);
  data.jsc <- data.frame(voltage[range.jsc], current[range.jsc]);
  fit2 <- lm(data.jsc$current ~ data.jsc$voltage + I(data.jsc$voltage^2));
  jsc <- coef(fit2)[1] / cell.surface * (as.numeric(current.positive)*2-1);
  return(jsc);
}

extract.iv.voc <- function(voltage, current){
  pos.voc <- which.min(abs(current));
  length <- length(voltage);
  range.voc <- max((pos.voc-3),1):min(pos.voc+3,length);
  data.voc <- data.frame(voltage[range.voc], current[range.voc]);
  fit2inv <- lm(data.voc$voltage ~ data.voc$current + I(data.voc$current^2));
  voc <- coef(fit2inv)[1];
  return(voc);
}


