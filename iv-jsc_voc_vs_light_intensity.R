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
directory=tail(strsplit(getwd(), "/")[[1]], n=1)
name=directory


library(robustbase)
a <- read.table("output.txt",header=F,stringsAsFactors=F)
#file <- strsplit(a$V1,"/")

illumination <- function(names){
values <- lapply(names, function(name){
splittedname <- strsplit(name, "-")
value=100
if(grepl("dark",name)){value=0}
if(grepl("sun",name)){value=100*as.numeric(gsub("sun","",unlist(splittedname)[grepl("sun",unlist(splittedname))]))}
return(value)})
return(unlist(values))
}

a$suns <- illumination(a$V2)
a<-a[with(a, order(a$suns)), ]

#suns <- as.numeric(gsub("sun","",unlist(splittedfilename)[grepl("sun",unlist(splittedfilename))]))

#if(length(a$V1) != length(suns))
#{ suns[length(suns)+1] <- 1; suns[length(suns)+1] <- 1}

#suns <- suns * 100

#a$suns <- suns

fitJsc<-nlrob(V4 ~ B*suns^C, start=list(B=20, C=1), data=a)
alfa=signif(fitJsc$coefficients["C"],3)

png(paste(name,"-Jsc_vs_LI.png",sep=""), width=640, height=640)
par(mar=c(7,7,4.1,2.1))
plot(a$suns,a$V4,cex.main=2,ylab=bquote("J"["sc"]~"(mA/cm"^"2"*")"),xlab="Light Intensity (%sun)", cex.lab=2,cex.axis=1.5)#, main=paste(name,"Jsc vs LI"))
lines(a$suns, predict(fitJsc))
mtext(bquote("Jsc" == .(signif(fitJsc$coefficients["B"],3)) * " LI" ^ .(alfa)),side=3,line=-3,cex=2)
mtext(bquote(alpha == .(alfa)),side=3,line=-5,cex=2)
graphics.off()

b=a[a$suns != 0,]
fitVoc<-lmrob(b$V5 ~ log(b$suns)) #log = ln; log10 = Log

nid <- signif(fitVoc$coefficients[2],3)

png(paste(name,"-Voc_vs_LI.png",sep=""), width=640, height=640)
par(mar=c(5.5,5,4.1,2.1))
plot(b$suns,b$V5,cex.main=2,ylab=bquote("V"["oc"]~"(V)"),xlab="Light Intensity (%sun)", cex.lab=2,cex.axis=1.5, log="x")#, main=paste(name,"Voc vs LI"))
lines(b$suns, predict(fitVoc))
mtext(bquote("Voc" == .(signif(fitVoc$coefficients[1],3)) + .(signif(nid,3)) ~ "ln(LI)"),side=3,line=-3,cex=2)
mtext(bquote("n"["id"] == .(signif(nid/0.02585,3))),side=3,line=-5,cex=2)
graphics.off()

