library(RColorBrewer)
library(sfsmisc)
library(Hmisc)
library(robustbase)

############################# Settings

filename = "************FILL ME HERE***********"
uvvis_dir = "data-UVvis"
pl_dir = "data-PL"
colorIndex = ************FILL ME HERE***********
  
  uvvis_wavelength_lim=c(250, 800)
pl_wavelength_lim=c(600, 850)

mylegend = c(paste(filename, ", experimental absorbance", sep=""), paste(filename, ", simulated absorbance", sep=""))
mylegend_Tauc = c(paste(filename, ", experimental Tauc plot", sep=""), paste(filename, ", simulated Tauc plot", sep=""))
mylegend_pl = c(paste(filename, ", experimental absorbance", sep=""), paste(filename, ", simulated absorbance", sep=""), paste(filename, ", photoluminescence", sep=""))
mylegend_Tauc_pl = c(paste(filename, "experimental Tauc plot"), paste(filename, "simulated Tauc plot"), paste(filename, "photoluminescence"))

mycolors=brewer.pal(8,"Dark2")

############################# General data UVvis

experimental_full <- read.csv(file.path(uvvis_dir, paste(filename, "-abs.csv",sep="")))

simulated_full <- read.csv(list.files(path=uvvis_dir, pattern=paste(filename, ".*simulatedabs.*\\.csv$", sep=""), full.names=T), header=T)

experimental = subset(experimental_full, Wavelength > uvvis_wavelength_lim[1] & Wavelength < uvvis_wavelength_lim[2])
simulated = subset(simulated_full, Wavelength > uvvis_wavelength_lim[1] & Wavelength < uvvis_wavelength_lim[2])

# I assume that the spectra reaches zero at some point
experimental$Abs = experimental$Abs - quantile(experimental$Abs, 0.05)
experimental$Epsilon = experimental$Epsilon - quantile(experimental$Epsilon, 0.05)

maxUVvisExpAbs = max(experimental$Abs)
experimental$Norm <- experimental$Abs/maxUVvisExpAbs
maxUVvisSim = max(simulated$Epsilon)
simulated$Norm <- simulated$Epsilon/maxUVvisSim
maxUVvisExpEpsilon = max(experimental$Epsilon)
maxUVvisEpsilon = max(maxUVvisExpEpsilon, maxUVvisExpEpsilon)

############################# General data photoluminescence

pl_full <- read.table(list.files(path=pl_dir, pattern=paste(filename, ".*\\.dat$", sep=""), full.names=T), header=T)[-1, ]
pl_full$Wavelength = as.numeric(as.character(pl_full$Wavelength))
pl_full$S1 = as.numeric(as.character(pl_full$S1))
#pl_full$S1c = as.numeric(as.character(pl_full$S1c))

print(pl_full)

pl = subset(pl_full, Wavelength > pl_wavelength_lim[1] & Wavelength < pl_wavelength_lim[2])

pl$Norm <- pl$S1/max(pl$S1)

pl$eV = 1240/pl$Wavelength

eV_pl_lim=c(1240/pl_wavelength_lim[2], 1240/pl_wavelength_lim[1])

pl_peak_eV = 1240/pl$Wavelength[which.max(pl$S1)]

############################# UVvis with photoluminescence data

pl$NormToEpsilon <- pl$Norm*maxUVvisEpsilon*0.7

############################# General data Tauc

experimental$eV = 1240/experimental$Wavelength
simulated$eV = 1240/simulated$Wavelength
eV_uvvis_lim=c(1240/uvvis_wavelength_lim[2], 1240/uvvis_wavelength_lim[1])

############################# Tauc direct data

experimental$taucDirect = (experimental$eV*experimental$Epsilon)^2
simulated$taucDirect = (simulated$eV*simulated$Epsilon)^2

fitDirectLim1 = (max(experimental$eV)-min(experimental$eV))*0.5 + min(experimental$eV)
fitDirectData1 = subset(experimental, eV<fitDirectLim1 & taucDirect<max(taucDirect)/2.5)
fitDirectData2 = subset(fitDirectData1, taucDirect>max(taucDirect)/4)
fitDirectLim2 = fitDirectData2$eV[which.max(fitDirectData2$taucDirect)]
fitDirectData = subset(fitDirectData2, eV<fitDirectLim2)

fitDirect=lmrob(taucDirect ~ eV, fitDirectData, control=lmrob.control(tuning.psi=1));
gapDirect <- -coef(fitDirect)[1]/coef(fitDirect)[2]
maxTaucDirect = max(fitDirectData$taucDirect)

fitSimDirectLim1 = (max(simulated$eV)-min(simulated$eV))*0.5 + min(simulated$eV)
fitSimDirectData1 = subset(simulated, eV<fitSimDirectLim1 & taucDirect<max(taucDirect)/2.5)
fitSimDirectData2 = subset(fitSimDirectData1, taucDirect>max(taucDirect)/4)
fitSimDirectLim2 = fitSimDirectData2$eV[which.max(fitSimDirectData2$taucDirect)]
fitSimDirectData = subset(fitSimDirectData2, eV<fitSimDirectLim2)

fitSimDirect=lmrob(taucDirect ~ eV, fitSimDirectData, control=lmrob.control(tuning.psi=1));
gapSimDirect <- -coef(fitSimDirect)[1]/coef(fitSimDirect)[2]
maxSimTaucDirect = max(fitSimDirectData$taucDirect)

############################# Tauc direct with photoluminescence data

pl$NormToTaucDirect <- pl$Norm*max(maxTaucDirect, maxSimTaucDirect)*0.7

############################# Tauc indirect data

# could have some NaN
experimental$taucIndirect = (experimental$eV*experimental$Epsilon)^0.5
simulated$taucIndirect = (simulated$eV*simulated$Epsilon)^0.5

# remove NaN due to neg values
fitIndirectData1 = experimental[!is.nan(experimental$taucIndirect),]

fitIndirectLim <- (max(fitIndirectData1$eV)-min(fitIndirectData1$eV))*0.5 + min(fitIndirectData1$eV)
fitIndirectData2 = subset(fitIndirectData1, eV<fitIndirectLim & taucIndirect<max(taucIndirect)/2.5)
fitIndirectData3 = subset(fitIndirectData2, taucIndirect>max(taucIndirect)/4)
fitIndirectLim2 = fitIndirectData3$eV[which.max(fitIndirectData3$taucIndirect)]
fitIndirectData = subset(fitIndirectData3, eV<fitIndirectLim2)

fitIndirect=lmrob(taucIndirect ~ eV, fitIndirectData, control=lmrob.control(tuning.psi=1));
gapIndirect <- -coef(fitIndirect)[1]/coef(fitIndirect)[2]
maxTaucIndirect = max(fitIndirectData$taucIndirect)

fitSimIndirectLim <- (max(simulated$eV)-min(simulated$eV))*0.5 + min(simulated$eV)
fitSimIndirectData1 = subset(simulated, eV<fitSimIndirectLim & taucIndirect<max(taucIndirect)/2.5)
fitSimIndirectData2 = subset(fitSimIndirectData1, taucIndirect>max(taucIndirect)/4)
fitSimIndirectLim2 = fitSimIndirectData2$eV[which.max(fitSimIndirectData2$taucIndirect)]
fitSimIndirectData = subset(fitSimIndirectData2, eV<fitSimIndirectLim2)

fitSimIndirect=lmrob(taucIndirect ~ eV, fitSimIndirectData, control=lmrob.control(tuning.psi=1));
gapSimIndirect <- -coef(fitSimIndirect)[1]/coef(fitSimIndirect)[2]
maxSimTaucIndirect = max(fitSimIndirectData$taucIndirect)

############################# UVvis

if(output_pdf){
  pdf(paste(filename,"-UVvis.pdf",sep=""), width=image_bigpdf_width, height=image_bigpdf_height, pointsize=7)
}else{
  png(paste(filename,"-UVvis.png",sep=""), width=image_width, height=image_height)
}
op <- par(mar = c(5,7.5,1,4) + 0.1) ## default is c(5,4,4,2) + 0.1 
plot(NULL,xlim=uvvis_wavelength_lim,ylim=c(0,1), xlab="", ylab="", cex.axis=1.4, yaxt="n");
title(ylab = "Normalized Absorption", cex.lab = 1.7, line = 5.5)
title(xlab = "Wavelength (nm)", cex.lab = 1.7, line = 3)
mtext("Photoluminescence (a.u.)", cex=1.7, side=4, line=2)

eaxis(side=2, cex.axis=1.4)
minor.tick(nx=10, ny=10)

lines(experimental$Wavelength, experimental$Norm, col=mycolors[colorIndex],lwd=2)
lines(simulated$Wavelength, simulated$Norm, col=mycolors[colorIndex],lwd=2, lty=2)
lines(pl$Wavelength, pl$Norm*0.7, col=mycolors[colorIndex],lwd=2, lty=3)

legend(x="topright", inset=0.05, mylegend_pl, col=mycolors[colorIndex], cex=1.5, lwd=2, lty=c(1,2,3), bty="n")
graphics.off()
#reset the plotting margins
par(op)

############################# UVvis epsilon

if(output_pdf){
  pdf(paste(filename,"-UVvis_epsilon.pdf",sep=""), width=image_bigpdf_width, height=image_bigpdf_height, pointsize=7)
}else{
  png(paste(filename,"-UVvis_epsilon.png",sep=""), width=image_width, height=image_height)
}
op <- par(mar = c(5,7.5,1,4) + 0.1) ## default is c(5,4,4,2) + 0.1 
plot(NULL,xlim=uvvis_wavelength_lim,ylim=c(0,max(simulated$Epsilon, experimental$Epsilon)), xlab="", ylab="", cex.axis=1.4, yaxt="n");
title(ylab = bquote("Molar attenuation coefficient (L mol"^"-1"~" cm"^"-1"~")"), cex.lab = 1.7, line = 5.5)
title(xlab = "Wavelength (nm)", cex.lab = 1.7, line = 3)
mtext("Photoluminescence (a.u.)", cex=1.7, side=4, line=2)

eaxis(side=2, cex.axis=1.4)
minor.tick(nx=10, ny=10)

lines(experimental$Wavelength, experimental$Epsilon, col=mycolors[colorIndex],lwd=2)
lines(simulated$Wavelength, simulated$Epsilon, col=mycolors[colorIndex],lwd=2, lty=2)
lines(pl$Wavelength, pl$NormToEpsilon, col=mycolors[colorIndex],lwd=2, lty=3)

legend(x="topright", inset=0.05, mylegend_pl, col=mycolors[colorIndex], cex=1.5, lwd=2, lty=c(1,2,3), bty="n")
graphics.off()
#reset the plotting margins
par(op)



############################# Tauc direct

if(output_pdf){
  pdf(paste(filename,"-UVvis_tauc_direct.pdf",sep=""), width=image_bigpdf_width, height=image_bigpdf_height, pointsize=7)
}else{
  png(paste(filename,"-UVvis_tauc_direct.png",sep=""), width=image_width, height=image_height)
}
op <- par(mar = c(5,7.5,1,1) + 0.1) ## default is c(5,4,4,2) + 0.1
xlimDirect=c(max(eV_uvvis_lim[1],min(fitSimDirectData$eV, fitDirectData$eV)-1), max(fitDirectLim2, fitSimDirectLim2))
ylimDirect=c(0,max(maxTaucDirect, maxSimTaucDirect))
plot(NULL,xlim=xlimDirect, ylim=ylimDirect, xlab="", ylab="", cex.axis=1.4, yaxt="n");
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
mtext(paste("Direct BG", signif(gapDirect,4), "eV"), side=3, line=-7, cex=1.5, adj=0.05)
mtext(paste("Sim Direct BG", signif(gapSimDirect,4), "eV"), side=3, line=-8.5, cex=1.5, adj=0.05)

legend(x="topleft", inset=0.05, mylegend_Tauc, col=mycolors[colorIndex], cex=1.5, lwd=2, lty=c(1,2), bty="n")
graphics.off()
#reset the plotting margins
par(op)



############################# Tauc indirect

if(output_pdf){
  pdf(paste(filename,"-UVvis_tauc_indirect.pdf",sep=""), width=image_bigpdf_width, height=image_bigpdf_height, pointsize=7)
}else{
  png(paste(filename,"-UVvis_tauc_indirect.png",sep=""), width=image_width, height=image_height)
}
op <- par(mar = c(5,7.5,1,1) + 0.1) ## default is c(5,4,4,2) + 0.1
xlimIndirect=c(max(eV_uvvis_lim[1], min(fitSimIndirectData$eV, fitIndirectData$eV)-1), max(fitIndirectLim2, fitSimIndirectLim2))
ylimIndirect=c(0,max(maxTaucIndirect, maxSimTaucIndirect))
plot(NULL,xlim=xlimIndirect, ylim=ylimIndirect, xlab="", ylab="", cex.axis=1.4, yaxt="n");
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
mtext(paste("Indirect BG", signif(gapIndirect,4), "eV"), side=3, line=-7.5, cex=1.5, adj=0.1)
mtext(paste("Sim Indirect BG", signif(gapSimIndirect,4), "eV"), side=3, line=-9, cex=1.5, adj=0.1)

legend(x="topleft", inset=0.05, mylegend_Tauc, col=mycolors[colorIndex], cex=1.5, lwd=2, lty=c(1,2), bty="n")
graphics.off()
#reset the plotting margins
par(op)


############################# Tauc direct with photoluminescence

if(output_pdf){
  pdf(paste(filename,"-UVvis_tauc_direct_with_PL.pdf",sep=""), width=image_bigpdf_width, height=image_bigpdf_height, pointsize=7)
}else{
  png(paste(filename,"-UVvis_tauc_direct_with_PL.png",sep=""), width=image_width, height=image_height)
}
op <- par(mar = c(5,3,1,7.5) + 0.1) ## default is c(5,4,4,2) + 0.1
xlimDirectPL=c(eV_pl_lim[1], max(fitDirectLim2, fitSimDirectLim2))
ylimDirectPL=c(0,max(maxTaucDirect, maxSimTaucDirect))
plot(NULL,xlim=xlimDirectPL, ylim=ylimDirectPL, xlab="", ylab="", cex.axis=1.4, yaxt="n");
title(ylab = "Photoluminescence (a.u.)", cex.lab = 1.7, line = 1)
title(xlab = "Photon energy (eV)", cex.lab = 1.7, line = 3)
mtext(bquote("(" ~ epsilon ~"h" ~nu ~")" ^ "2"), cex=1.7, side=4, line=6)

eaxis(side=4, cex.axis=1.4)
minor.tick(nx=10, ny=0)

lines(experimental$eV, experimental$taucDirect, col=mycolors[colorIndex],lwd=2)
lines(simulated$eV, simulated$taucDirect, col=mycolors[colorIndex],lwd=2, lty=2)

abline(fitDirect, col="gray")
abline(fitSimDirect, col="gray")

lines(pl$eV, pl$NormToTaucDirect, col=mycolors[colorIndex],lwd=2, lty=3)

abline(v=pl_peak_eV, col="gray")

abline(h=0);
mtext(paste("Exp. Direct Optical BG", signif(gapDirect,4), "eV"), side=3, line=-7.5, cex=1.5, adj=0.05)
mtext(paste("Sim. Direct Optical BG", signif(gapSimDirect,4), "eV"), side=3, line=-9, cex=1.5, adj=0.05)
mtext(paste("PL peak", signif(pl_peak_eV,4), "eV"), side=3, line=-10.5, cex=1.5, adj=0.05)

legend(x="topleft", inset=0, mylegend_Tauc_pl, col=mycolors[colorIndex], cex=1.5, lwd=2, lty=c(1,2,3), bty="n")
graphics.off()
#reset the plotting margins
par(op)

############################# Tauc indirect with photoluminescence

if(output_pdf){
  pdf(paste(filename,"-UVvis_tauc_indirect_with_PL.pdf",sep=""), width=image_bigpdf_width, height=image_bigpdf_height, pointsize=7)
}else{
  png(paste(filename,"-UVvis_tauc_indirect_with_PL.png",sep=""), width=image_width, height=image_height)
}
op <- par(mar = c(5,3,1,7.5) + 0.1) ## default is c(5,4,4,2) + 0.1
xlimIndirectPL=c(eV_pl_lim[1], fitIndirectLim2)
ylimIndirectPL=c(0,maxTaucIndirect)
plot(NULL,xlim=xlimIndirectPL, ylim=ylimIndirectPL, xlab="", ylab="", cex.axis=1.4, yaxt="n");
title(ylab = "Photoluminescence (a.u.)", cex.lab = 1.7, line=1)
title(xlab = "Photon energy (eV)", cex.lab = 1.7, line = 3)
mtext(bquote("(" ~ epsilon ~"h" ~nu ~")" ^ "0.5"), cex=1.7, side=4, line=5)

eaxis(side=4, cex.axis=1.4)
minor.tick(nx=10, ny=0)

lines(experimental$eV, experimental$taucIndirect, col=mycolors[colorIndex],lwd=2)
lines(simulated$eV, simulated$taucIndirect, col=mycolors[colorIndex],lwd=2, lty=2)

abline(fitIndirect, col="gray")
abline(fitSimIndirect, col="gray")

lines(pl$eV, pl$NormToTaucIndirect, col=mycolors[colorIndex],lwd=2, lty=3)

abline(v=pl_peak_eV, col="gray")

abline(h=0);
mtext(paste("Exp. Indirect Optical BG", signif(gapIndirect,4), "eV"), side=3, line=-7.5, cex=1.5, adj=0.05)
mtext(paste("Sim. Indirect Optical BG", signif(gapSimIndirect,4), "eV"), side=3, line=-9, cex=1.5, adj=0.05)
mtext(paste("PL peak", signif(pl_peak_eV,4), "eV"), side=3, line=-10.5, cex=1.5, adj=0.05)

legend(x="topleft", inset=0, mylegend_Tauc_pl, col=mycolors[colorIndex], cex=1.5, lwd=2, lty=c(1,2,3), bty="n")
graphics.off()
#reset the plotting margins
par(op)
