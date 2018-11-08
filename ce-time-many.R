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

#title=gsub("_"," ",tail(unlist(strsplit(name,"-")),1))
filename=gsub(",","",gsub(":","",name))

library(RColorBrewer)
library(robustbase)
library(sfsmisc)
library(Hmisc)
require(minpack.lm)

output=list()
ce.max = 0

i <- 0
dirs <- list.dirs(recursive=FALSE)
dirs <- sub("./","",dirs)
legend=sub("_.*","",sub("^0","",dirs))

# try to obtain the color from the file name
mycolors=gsub(".*-col_","",dirs[grepl("-col_", dirs)])
# if the color is not set, use the default one
if(!length(mycolors)){mycolors=brewer.pal(8,"Dark2")}

data <- lapply(dirs, function(x) {print(x);
 subdirs <- list.dirs(path=x, recursive=F)
 subdirs.ce <- subdirs[grep("ce", subdirs, ignore.case=T)]
ce.files <- list.files(path=subdirs.ce, pattern="^CE.*\\.txt.table$");
ce.file <- tail(ce.files, n=1)
print(ce.file)
print(x)
ce <- read.table(file.path(subdirs.ce, ce.file), header=F)
ce.max <<- max(ce.max, max(ce$V2))
 ce})
names(data) <- dirs

ylim=c(-0.05*ce.max, ce.max)
xlim=c(-5e-8,5e-6)

png(paste(filename, "-CE_times.png", sep=""), width=image_width, height=image_height)
par(mar=c(5.1,5,2,2.1))
plot(NULL, xlim=xlim, ylim=ylim, cex.main=1.5, xlab="Time (s)", ylab="Current (a.u.)", cex.lab=1.5, cex.axis=1.2, yaxt="n")
lapply(dirs, function(x) {print(x);
 lines(data[[x]]$V1, data[[x]]$V2, col=mycolors[i+1], lwd=2)
 i <<- i+1
})
legend(x="topright",inset=0.05, legend, col=mycolors, cex=1.5, lwd=4, title=title, bg="gray90"
)
graphics.off()



