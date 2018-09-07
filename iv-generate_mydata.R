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

source("~/software/photophysics-data-processing-R/extractdata-curves-vi-separated_files.R")

mydata <- import.iv.separated(pattern.excl="\\.png$", pattern=".*forward.*.txt$|.*reverse.*.txt", list.excl="output.txt")
#mydata <- import.iv.separated(pattern.excl="\\.pdf$")
file.create("output.txt")
results <- lapply(names(mydata), function(x){
	write(extract.iv(mydata[[x]]$Voltage_V, mydata[[x]]$Current_mA, cell.surface=0.09, formatted.output=TRUE, directory=getwd(), sample=x, reverse=as.integer(grepl("reverse", x))), file="output.txt", append=TRUE);
	extract.iv(mydata[[x]]$Voltage_V, mydata[[x]]$Current_mA, comment="", cell.surface=0.09)})
names(results) <- names(mydata)
