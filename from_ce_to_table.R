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

fromCeToTable <- function(cedir="ce")
{
  print("CE: ADDING TIME COLUMN IN .table FILES")
  files <- list.files(path=cedir, pattern="^CE.*\\.txt$");
  mydata <- lapply(file.path(cedir, files), read.table, header=FALSE, skip=9, col.names=c("voltage"));
  names(mydata) <- files;
  
  trashfornullmessages <- lapply(files, function(x) {
    #	if(!file.exists(paste(x, ".table", sep=""))){
    message(x);
    header = read.table(file.path(cedir, x), skip=4, header=FALSE, nrows=3)$V2
    #		header = as.numeric(system(paste("head -7 '", file.path(cedir, x), "' | tail -3|sed 's/\r$//' | cut -f2 -d' '", sep=""), intern = TRUE))
    mydata[[x]]$time = seq(0,header[3]-1)*header[1]+header[2];
    write.table(mydata[[x]][,c("time","voltage")], file.path(cedir, paste(x,".table", sep="")), col.names=F, row.names=F, quote=F)
    #	}
  })
}
