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

#for i in `ls`; do cp ../$i/tpv/output* $i; done

#name=""

#title=gsub("_"," ",tail(unlist(strsplit(name,"-")),1))
filename=gsub(",","",gsub(":","",name))

ylim=lim.TPV.lifetime
xlim=lim.TPV.voltage

output.biexp=list()
output.monoexp=list()

#library(Hmisc)
#library(RColorBrewer)
library(sfsmisc)
library(Hmisc)


dirs <- list.dirs(recursive=FALSE)
dirs <- sub("./","",dirs)
legend=sub("_.*","",sub("^0","",dirs))

# try to obtain the color from the file name
colors=gsub(".*-col_","",dirs)
# if the color is not set, use the default one
if(!length(colors[1])){colors=colorRampPalette(c("red","orange","springgreen","royalblue"))(max(length(dirs),3))}

print("errori derivano da avere molti header in output, bisogna pulirlo")

print("biexp")

i <- 0
jpeg(quality=98, paste(filename, "-TPVs-biexp.jpg",sep=""), width=image_width, height=image_height);
op <- par(mar = c(5,7,4,2) + 0.1) ## default is c(5,4,4,2) + 0.1
plot(NULL, ylim=ylim, xlim=xlim, cex.main=1.5,cex.axis=1.5,cex.lab=2, log="y", xlab="", ylab="", yaxt="n")#, main=paste(name, "TPV biexp"));
title(ylab = "Charge carrier lifetime (s)", cex.lab = 2, line = 4)
title(xlab = "Voltage (V)", cex.lab = 2, line = 4)

eaxis(side=2,at=c(1e-10,1e-9,1e-8,1e-7,1e-6,1e-5,1e-4,1e-3,1e-2,0.1,1,10,100,1e3), cex.axis=1.5)
minor.tick(nx=10)

lapply(dirs, function(x) {print(x);
fulloutput <- read.table(paste(x,"/tpv/output-biexp.txt",sep=""), header=TRUE);#,stringsAsFactors=F);
n<-tail(grep("file",fulloutput[,1]),n=1)
output <- read.table(paste(x,"/tpv/output-biexp.txt",sep=""), header=TRUE, skip=ifelse(length(n),n,0)); 
output.biexp[[paste("Voc",sub("nm","",sub("_.*","",sub("^0","",x))),sep="")]] <<- signif(output$Voc,5)
output.biexp[[paste(sub("_.*","",sub("^0","",x)),"T1",sep="")]] <<- signif(output$T1,5)
output.biexp[[paste(sub("_.*","",sub("^0","",x)),"T2",sep="")]] <<- signif(output$T2,5)
points(output$Voc, output$T1, lwd=1, pch=21+(i%%5), bg=colors[i+1], cex=2);
points(output$Voc, output$T2, lwd=1, pch=21+(i%%5), col=colors[i+1], cex=2);
i <<- i+1
})
legend(x="topright",inset=0.05,legend,pt.cex=2, pt.lwd=2, cex=2, pch=seq(21,25), pt.bg=colors,title=title, bty="n", col=colors)
graphics.off()
#reset the plotting margins
par(op)

maxlength.biexp = max(sapply(output.biexp,length))
output.biexp = lapply(output.biexp, function(x){length(x)=maxlength.biexp; print(x)})
output.biexp = as.data.frame(output.biexp,check.names=FALSE)
write.table(output.biexp, file=paste(filename,"-TPVs-biexp.csv",sep=""), row.names=FALSE, na="", sep=",")


print("monoexp")
i <- 0
jpeg(quality=98, paste(filename, "-TPVs-monoexp.jpg",sep=""), width=image_width, height=image_height);
op <- par(mar = c(5,7,4,2) + 0.1) ## default is c(5,4,4,2) + 0.1
plot(NULL, ylim=ylim, xlim=xlim, cex.main=1.5,cex.axis=1.5,cex.lab=2, log="y", xlab="", ylab="", yaxt="n")#, main=paste(name, "TPV monoexp"));
title(ylab = "Charge carrier lifetime (s)", cex.lab = 2, line = 4)
title(xlab = "Voltage (V)", cex.lab = 2, line = 4)

eaxis(side=2,at=c(1e-10,1e-9,1e-8,1e-7,1e-6,1e-5,1e-4,1e-3,1e-2,0.1,1,10,100,1e3), cex.axis=1.5)
minor.tick(nx=10)

lapply(dirs, function(x) {print(x);
fulloutput <- read.table(paste(x,"/tpv/output-monoexp.txt",sep=""), header=TRUE);
n<-tail(grep("file",fulloutput[,1]),n=1)
output <- read.table(paste(x,"/tpv/output-monoexp.txt",sep=""), header=TRUE, skip=ifelse(length(n),n,0)); 
output.monoexp[[paste("Voc",sub("nm","",sub("_.*","",sub("^0","",x))),sep="")]] <<- signif(output$Voc,5)
output.monoexp[[sub("_.*","",sub("^0","",x))]] <<- signif(output$T,5)
points(output$Voc, output$T, lwd=1, pch=21+(i%%5), bg=colors[i+1], cex=2);
i <<- i+1
})
legend(x="topright",inset=0.05,legend,pt.cex=2, pt.lwd=2,cex=2, pch=seq(21,25), pt.bg=colors, col=colors, title=title, bty="n")

graphics.off()
#reset the plotting margins
par(op)

maxlength.monoexp = max(sapply(output.monoexp,length))
output.monoexp = lapply(output.monoexp, function(x){length(x)=maxlength.monoexp})
output.monoexp = as.data.frame(output.monoexp,check.names=FALSE)
write.table(output.monoexp, file=paste(filename,"-TPVs-monoexp.csv",sep=""), row.names=FALSE, na="", sep=",")

