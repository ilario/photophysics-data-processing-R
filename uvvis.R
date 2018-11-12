library(RColorBrewer)
library(sfsmisc)
library(Hmisc)
library(robustbase)

filename = "************FILL ME HERE***********"
colorIndex = ***********FILL ME HERE***********

wavelength_lim=c(250, 800)

mylegend= c(paste(filename, ", experimental", sep=""), paste(filename, ", simulated", sep=""))

mycolors=brewer.pal(8,"Dark2")

experimental_full <- read.csv(paste(filename, "-abs.csv",sep=""))

simulated_full <- read.csv(list.files(path=".", pattern=paste(filename, ".*simulatedabs.*\\.csv$", sep="")), header=T)

experimental = subset(experimental_full, Wavelength > wavelength_lim[1] & Wavelength < wavelength_lim[2])
simulated = subset(simulated_full, Wavelength > wavelength_lim[1] & Wavelength < wavelength_lim[2])

# I assume that the spectra reaches zero at some point
experimental$Abs = experimental$Abs - quantile(experimental$Abs, 0.05)
experimental$Epsilon = experimental$Epsilon - quantile(experimental$Epsilon, 0.05)

experimental$Norm <- experimental$Abs/max(experimental$Abs)
simulated$Norm <- simulated$Epsilon/max(simulated$Epsilon)

experimental$eV = 1240/experimental$Wavelength
simulated$eV = 1240/simulated$Wavelength
eV_lim=c(1240/wavelength_lim[2], 1240/wavelength_lim[1])

experimental$taucDirect = (experimental$eV*experimental$Norm)^2
simulated$taucDirect = (simulated$eV*simulated$Norm)^2

# could have some NaN
experimental$taucIndirect = (experimental$eV*experimental$Norm)^0.5
simulated$taucIndirect = (simulated$eV*simulated$Norm)^0.5

if(output_pdf){
	pdf(paste(filename,"-UVvis.pdf",sep=""), width=image_bigpdf_width, height=image_bigpdf_height, pointsize=7)
}else{
	png(paste(filename,"-UVvis.png",sep=""), width=image_width, height=image_height)
}
op <- par(mar = c(5,7.5,1,1) + 0.1) ## default is c(5,4,4,2) + 0.1 
plot(NULL,xlim=wavelength_lim,ylim=c(0,1), xlab="", ylab="", cex.axis=1.4, yaxt="n");
title(ylab = "Normalized Absorption", cex.lab = 1.7, line = 5.5)
title(xlab = "Wavelength (nm)", cex.lab = 1.7, line = 3)

eaxis(side=2, cex.axis=1.4)
minor.tick(nx=10, ny=10)

lines(experimental$Wavelength, experimental$Norm, col=mycolors[colorIndex],lwd=2)
lines(simulated$Wavelength, simulated$Norm, col=mycolors[colorIndex],lwd=2, lty=2)

legend(x="topright", inset=0.05, mylegend, col=mycolors[colorIndex], cex=1.5, lwd=2, lty=c(1,2), bty="n")
graphics.off()
#reset the plotting margins
par(op)



if(output_pdf){
	pdf(paste(filename,"-UVvis_epsilon.pdf",sep=""), width=image_bigpdf_width, height=image_bigpdf_height, pointsize=7)
}else{
	png(paste(filename,"-UVvis_epsilon.png",sep=""), width=image_width, height=image_height)
}
op <- par(mar = c(5,7.5,1,1) + 0.1) ## default is c(5,4,4,2) + 0.1 
plot(NULL,xlim=wavelength_lim,ylim=c(0,max(simulated$Epsilon, experimental$Epsilon)), xlab="", ylab="", cex.axis=1.4, yaxt="n");
title(ylab = bquote("Molar attenuation coefficient (L mol"^"-1"~" cm"^"-1"~")"), cex.lab = 1.7, line = 5.5)
title(xlab = "Wavelength (nm)", cex.lab = 1.7, line = 3)

eaxis(side=2, cex.axis=1.4)
minor.tick(nx=10, ny=10)

lines(experimental$Wavelength, experimental$Epsilon, col=mycolors[colorIndex],lwd=2)
lines(simulated$Wavelength, simulated$Epsilon, col=mycolors[colorIndex],lwd=2, lty=2)

legend(x="topright", inset=0.05, mylegend, col=mycolors[colorIndex], cex=1.5, lwd=2, lty=c(1,2), bty="n")
graphics.off()
#reset the plotting margins
par(op)




fitDirectLim1 = (max(experimental$eV)-min(experimental$eV))*0.5 + min(experimental$eV)
fitDirectData1 = subset(experimental, eV<fitDirectLim1 & taucDirect<max(taucDirect)/2.5)
fitDirectData2 = subset(fitDirectData1, taucDirect>max(taucDirect)/4)
fitDirectLim2 = fitDirectData2$eV[which.max(fitDirectData2$taucDirect)]
fitDirectData = subset(fitDirectData2, eV<fitDirectLim2)

fitDirect=lmrob(taucDirect ~ eV, fitDirectData, control=lmrob.control(tuning.psi=2));
gapDirect <- -coef(fitDirect)[1]/coef(fitDirect)[2]

fitSimDirectLim1 = (max(simulated$eV)-min(simulated$eV))*0.5 + min(simulated$eV)
fitSimDirectData1 = subset(simulated, eV<fitSimDirectLim1 & taucDirect<max(taucDirect)/2.5)
fitSimDirectData2 = subset(fitSimDirectData1, taucDirect>max(taucDirect)/4)
fitSimDirectLim2 = fitSimDirectData2$eV[which.max(fitSimDirectData2$taucDirect)]
fitSimDirectData = subset(fitSimDirectData2, eV<fitSimDirectLim2)

fitSimDirect=lmrob(taucDirect ~ eV, fitSimDirectData, control=lmrob.control(tuning.psi=2));
gapSimDirect <- -coef(fitSimDirect)[1]/coef(fitSimDirect)[2]

if(output_pdf){
	pdf(paste(filename,"-UVvis_tauc_direct.pdf",sep=""), width=image_bigpdf_width, height=image_bigpdf_height, pointsize=7)
}else{
	png(paste(filename,"-UVvis_tauc_direct.png",sep=""), width=image_width, height=image_height)
}
op <- par(mar = c(5,7.5,1,1) + 0.1) ## default is c(5,4,4,2) + 0.1
plot(NULL,xlim=c(eV_lim[1], max(fitDirectLim2, fitSimDirectLim2)), ylim=c(0,max(fitDirectData$taucDirect, fitSimDirectData$taucDirect)), xlab="", ylab="", cex.axis=1.4, yaxt="n");
title(ylab = "", cex.lab = 1.7, line = 5.5)
title(xlab = "Photon energy (eV)", cex.lab = 1.7, line = 3)

eaxis(side=2, cex.axis=1.4)
minor.tick(nx=10, ny=10)

lines(experimental$eV, experimental$taucDirect, col=mycolors[colorIndex],lwd=2)
lines(simulated$eV, simulated$taucDirect, col=mycolors[colorIndex],lwd=2, lty=2)
lines(fitDirectData$eV, fitDirectData$taucDirect)
lines(fitSimDirectData$eV, fitSimDirectData$taucDirect)

print(fitDirect)
abline(fitDirect)

print(fitSimDirect)
abline(fitSimDirect)

abline(h=0);
mtext(paste("Direct BG", signif(gapDirect,4), "eV"), side=3, line=-4, cex=1.5)
mtext(paste("Sim Direct BG", signif(gapSimDirect,4), "eV"), side=3, line=-5.5, cex=1.5)

legend(x="topleft", inset=0.05, mylegend, col=mycolors[colorIndex], cex=1.5, lwd=2, lty=c(1,2), bty="n")
graphics.off()
#reset the plotting margins
par(op)



# remove NaN due to neg values
fitIndirectData1 = experimental[!is.nan(experimental$taucIndirect),]

fitIndirectLim <- (max(fitIndirectData1$eV)-min(fitIndirectData1$eV))*0.5 + min(fitIndirectData1$eV)
fitIndirectData2 = subset(fitIndirectData1, eV<fitIndirectLim & taucIndirect<max(taucIndirect)/2.5)
fitIndirectData3 = subset(fitIndirectData2, taucIndirect>max(taucIndirect)/4)
fitIndirectLim2 = fitIndirectData3$eV[which.max(fitIndirectData3$taucIndirect)]
fitIndirectData = subset(fitIndirectData3, eV<fitIndirectLim2)

fitIndirect=lmrob(taucIndirect ~ eV, fitIndirectData, control=lmrob.control(tuning.psi=2));
gapIndirect <- -coef(fitIndirect)[1]/coef(fitIndirect)[2]


fitSimIndirectLim <- (max(simulated$eV)-min(simulated$eV))*0.5 + min(simulated$eV)
fitSimIndirectData1 = subset(simulated, eV<fitSimIndirectLim & taucIndirect<max(taucIndirect)/2.5)
fitSimIndirectData2 = subset(fitSimIndirectData1, taucIndirect>max(taucIndirect)/4)
fitSimIndirectLim2 = fitSimIndirectData2$eV[which.max(fitSimIndirectData2$taucIndirect)]
fitSimIndirectData = subset(fitSimIndirectData2, eV<fitSimIndirectLim2)

fitSimIndirect=lmrob(taucIndirect ~ eV, fitSimIndirectData, control=lmrob.control(tuning.psi=2));
gapSimIndirect <- -coef(fitSimIndirect)[1]/coef(fitSimIndirect)[2]

if(output_pdf){
	pdf(paste(filename,"-UVvis_tauc_indirect.pdf",sep=""), width=image_bigpdf_width, height=image_bigpdf_height, pointsize=7)
}else{
	png(paste(filename,"-UVvis_tauc_indirect.png",sep=""), width=image_width, height=image_height)
}
op <- par(mar = c(5,7.5,1,1) + 0.1) ## default is c(5,4,4,2) + 0.1
plot(NULL,xlim=c(eV_lim[1], max(fitIndirectLim2, fitSimIndirectLim2)), ylim=c(0,max(fitIndirectData$taucIndirect, fitSimIndirectData$taucIndirect)), xlab="", ylab="", cex.axis=1.4, yaxt="n");
title(ylab = "", cex.lab = 1.7, line = 5.5)
title(xlab = "Photon energy (eV)", cex.lab = 1.7, line = 3)
lines(fitIndirectData$eV, fitIndirectData$taucIndirect)
eaxis(side=2, cex.axis=1.4)
minor.tick(nx=10, ny=10)

lines(experimental$eV, experimental$taucIndirect, col=mycolors[colorIndex],lwd=2)
lines(simulated$eV, simulated$taucIndirect, col=mycolors[colorIndex],lwd=2, lty=2)
lines(fitIndirectData$eV, fitIndirectData$taucIndirect)
lines(fitSimIndirectData$eV, fitSimIndirectData$taucIndirect)

print(fitIndirect)
abline(fitIndirect)

print(fitSimIndirect)
abline(fitSimIndirect)

abline(h=0);
mtext(paste("Indirect BG", signif(gapIndirect,4), "eV"), side=3, line=-4, cex=1.5)
mtext(paste("Sim Indirect BG", signif(gapSimIndirect,4), "eV"), side=3, line=-5.5, cex=1.5)

legend(x="topleft", inset=0.05, mylegend, col=mycolors[colorIndex], cex=1.5, lwd=2, lty=c(1,2), bty="n")
graphics.off()
#reset the plotting margins
par(op)
