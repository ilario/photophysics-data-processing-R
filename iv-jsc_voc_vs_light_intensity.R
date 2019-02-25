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

library(robustbase)
library(minpack.lm)
a <- read.table("output.txt",header=F,stringsAsFactors=F)
names(a) = c("dir", "file", "reverse", "jsc", "voc", "ff", "pce")

illumination <- function(names){
  values <- lapply(names, function(name){
    splittedname <- strsplit(name, "-")
    value=100
    if(grepl("dark",name, ignore.case=T)){value=0}
    if(grepl("sun",name, ignore.case=T)){value=100*as.numeric(gsub("sun","",unlist(splittedname)[grepl("sun",unlist(splittedname), ignore.case=T)]))}
    return(value)})
  return(unlist(values))
}

a$suns <- illumination(a$file)
a<-a[with(a, order(a$suns)), ]

fitJsc<-nlrob(jsc ~ B*suns^C, start=list(B=20, C=1), data=a)
alfa=signif(fitJsc$coefficients["C"],3)

png(paste(name,"-Jsc_vs_LI.png",sep=""), width=640, height=640)
par(mar=c(7,7,4.1,2.1))
plot(a$suns,a$jsc,cex.main=2,ylab=bquote("J"["SC"]~"(mA/cm"^"2"*")"),xlab="Light Intensity (%sun)", cex.lab=2,cex.axis=1.5)
lines(a$suns, predict(fitJsc))
mtext(bquote("J"["SC"] == .(signif(fitJsc$coefficients["B"],3)) * " LI%" ^ .(alfa)),side=3,line=-3,cex=2)
mtext(bquote(alpha == .(alfa)),side=3,line=-5,cex=2)
graphics.off()

b=a[a$suns != 0,]
print(b)
fitVoc<-lmrob(b$voc ~ log(b$suns))

nid <- fitVoc$coefficients[2]

png(paste(name,"-Voc_vs_LI.png",sep=""), width=640, height=640)
par(mar=c(5.5,5,4.1,2.1))
plot(b$suns,b$voc,cex.main=2,ylab=bquote("V"["OC"]~"(V)"),xlab="Light Intensity (%sun)", cex.lab=2,cex.axis=1.5, log="x")
lines(b$suns, predict(fitVoc))
mtext(bquote("V"["OC"] == .(signif(fitVoc$coefficients[1],3)) + .(signif(nid,3)) ~ "ln(LI%)"),side=3,line=-3,cex=2)
mtext(bquote("n"["id"] == .(signif(nid/0.02585,3))),side=3,line=-5,cex=2)
graphics.off()

bReverse = b[as.logical(b$reverse),]
print(bReverse)
fitVocReverse<-lmrob(bReverse$voc ~ log(bReverse$suns))
nidReverse <- fitVocReverse$coefficients[2]

png(paste(name,"-VocReverse_vs_LI.png",sep=""), width=640, height=640)
par(mar=c(5.5,5,4.1,2.1))
plot(bReverse$suns,bReverse$voc,cex.main=2,ylab=bquote("V"["OC"]~"(V)"),xlab="Light Intensity (%sun)", cex.lab=2,cex.axis=1.5, log="x")
lines(bReverse$suns, predict(fitVocReverse))
mtext(bquote("V"["OC"] == .(signif(fitVocReverse$coefficients[1],3)) + .(signif(nidReverse,3)) ~ "ln(LI%)"),side=3,line=-3,cex=2)
mtext(bquote("n"["id"] == .(signif(nidReverse/0.02585,3))),side=3,line=-5,cex=2)
graphics.off()

fitVoc2<-nlsLM(voc ~ A*log(suns/D + 1), start=list(A=3E-2, D=4E-10), data=b, lower=c(0,1E-30), algorithm="port")
tryCatch({
	fitVoc2<-nlrob(voc ~ A*log(suns/D + 1), start=coef(fitVoc2), data=b, lower=c(0,1E-30), algorithm="port")
}, error=function(e) cat("Robust fit: Error ", e$message, "\n"))
coef2 = coef(fitVoc2)

png(paste(name,"-Voc_vs_LI2.png",sep=""), width=640, height=640)
par(mar=c(5.5,5,4.1,2.1))
plot(b$suns,b$voc,cex.main=2,ylab=bquote("V"["OC"]~"(V)"),xlab="Light Intensity (%sun)", cex.lab=2,cex.axis=1.5, log="x")
lines(b$suns, predict(fitVoc2))
mtext(bquote("V"["OC"] == .(signif(coef2[["A"]],3)) ~ "ln(LI% / "~.(signif(coef2[["D"]],3))~" + 1)"),side=3,line=-3,cex=2)
mtext(bquote("n"["id"] == .(signif(coef2[["A"]]/0.02585,3))),side=3,line=-5,cex=2)
graphics.off()

fitVocReverse2 <- nlsLM(voc ~ A*log(suns/D + 1), start=list(A=4E-2, D=4E-10), data=bReverse, lower=c(0,1E-30), algorithm="port")
tryCatch({
	fitVocReverse2<-nlrob(voc ~ A*log(suns/D + 1), start=coef(fitVocReverse2), data=bReverse, lower=c(0,1E-30), algorithm="port")
}, error=function(e) cat("Robust fit: Error ", e$message, "\n"))
coefReverse2 <- coef(fitVocReverse2)
print(coefReverse2)

png(paste(name,"-VocReverse_vs_LI2.png",sep=""), width=640, height=640)
par(mar=c(5.5,5,4.1,2.1))
plot(bReverse$suns,bReverse$voc,cex.main=2,ylab=bquote("V"["OC"]~"(V)"),xlab="Light Intensity (%sun)", cex.lab=2,cex.axis=1.5, log="x")
lines(bReverse$suns, predict(fitVocReverse2))
mtext(bquote("V"["OC"] == .(signif(coefReverse2[["A"]],3)) ~ "ln(LI% / "~.(signif(coefReverse2[["D"]],3))~" + 1)"),side=3,line=-3,cex=2)
mtext(bquote("n"["id"] == .(signif(coefReverse2[["A"]]/0.02585,3))),side=3,line=-5,cex=2)
graphics.off()

fitVoc3 <- nlsLM(voc ~ A*log(jsc/E + 1), start=list(A=4E-2, E=4E-10), data=b, lower=c(0,1E-30), algorithm="port")
tryCatch({
	fitVoc3<-nlrob(voc ~ A*log(jsc/E + 1), start=coef(fitVoc3), data=b, lower=c(0,1E-30), algorithm="port")
}, error=function(e) cat("Robust fit: Error ", e$message, "\n"))
coef3 <- coef(fitVoc3)
print(coef3)

png(paste(name,"-Voc_vs_LI3.png",sep=""), width=640, height=640)
par(mar=c(5.5,5,4.1,2.1))
plot(b$jsc,b$voc,cex.main=2,ylab=bquote("V"["OC"]~"(V)"),xlab=bquote("J"["SC"]~" (mA/cm"^"2"*")"), cex.lab=2,cex.axis=1.5, log="x")
lines(b$jsc, predict(fitVoc3))
mtext(bquote("V"["OC"] == .(signif(coef3[["A"]],3)) ~ "ln(J"["SC"]~" / "~.(signif(coef3[["E"]]))~" + 1)"),side=3,line=-3,cex=2)
mtext(bquote("n"["id"] == .(signif(coef3[["A"]]/0.02585,3))),side=3,line=-5,cex=2)
graphics.off()

fitVocReverse3 <- nlsLM(voc ~ A*log(jsc/E + 1), start=list(A=4E-2, E=4E-10), data=bReverse, lower=c(0,1E-30), algorithm="port")
tryCatch({
	fitVocReverse3<-nlrob(voc ~ A*log(jsc/E + 1), start=coef(fitVocReverse3), data=bReverse, lower=c(0,1E-30), algorithm="port")
}, error=function(e) cat("Robust fit: Error ", e$message, "\n"))
coefReverse3 <- coef(fitVocReverse3)

png(paste(name,"-VocReverse_vs_LI3.png",sep=""), width=640, height=640)
par(mar=c(5.5,5,4.1,2.1))
plot(bReverse$jsc,bReverse$voc,cex.main=2,ylab=bquote("V"["OC"]~"(V)"),xlab=bquote("J"["SC"]~" (mA/cm"^"2"*")"), cex.lab=2,cex.axis=1.5, log="x")
lines(bReverse$jsc, predict(fitVocReverse3))
mtext(bquote("V"["OC"] == .(signif(coefReverse3[["A"]],3)) ~ "ln(J"["SC"]~" / "~.(signif(coefReverse3[["E"]]))~" + 1)"),side=3,line=-3,cex=2)
mtext(bquote("n"["id"] == .(signif(coefReverse3[["A"]]/0.02585,3))),side=3,line=-5,cex=2)
graphics.off()
