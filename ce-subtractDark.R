ceSubtractDark <- function(cedir="ce")
{
library(sfsmisc)
library(robustbase)
library(minpack.lm)

timeStart = 3e-8
noiseEndTime = 3e-7

rm("startList3")
options(error=function() { traceback(2); if(!interactive()) quit("no", status = 1, runLast = FALSE) })

print("CE: INTEGRATING")
files <- list.files(path=cedir, pattern="^CE.*\\.txt.table$");
mydata <- lapply(file.path(cedir,files), read.table, header=FALSE, col.names=c("time","voltage"));
files <- sub(".txt.table","",files);
names(mydata) <- files;
write.table(t(c("Voc","ChargeDensityCE")), file=file.path(cedir,"outputChargeDensityCE.txt"), append=FALSE, col.names=F, row.names=F);
write.table(t(c("Voc","CEmonoexpTime")), file=file.path(cedir,"outputMonoexpCE.txt"), append=FALSE, col.names=F, row.names=F);

print(files[1])
darkCEvoltageNoBaseline = mydata[[files[1]]]$voltage;
len<-length(darkCEvoltageNoBaseline)
darkStartVoltage <- mean(darkCEvoltageNoBaseline[1:600])
darkEndVoltage <- mean(tail(darkCEvoltageNoBaseline,-600))
darkBaseline <- seq(darkStartVoltage, darkEndVoltage, length.out=len)
darkCEvoltage <- darkCEvoltageNoBaseline - darkBaseline
darkCEtime = mydata[[files[1]]]$time;
deltaT = darkCEtime[2] - darkCEtime[1]

timeZeroIndex = match(0,darkCEtime);

timeStartIndex = which.min(abs(darkCEtime-timeStart))
darkCEtimeDecay = tail(darkCEtime, -timeStartIndex)
darkCEvoltageDecay = tail(darkCEvoltage, -timeStartIndex)

timeEndNoiseDecayIndex = which.min(abs(darkCEtimeDecay-noiseEndTime))
darkCEtimeDecayHead = head(darkCEtimeDecay, timeEndNoiseDecayIndex)

tryCatch({
	expfitDarkCE <- nlsLM(darkCEvoltageDecay ~ exp(C)*exp(-darkCEtimeDecay/D), start=list(C=log(max(darkCEvoltageDecay)),D=0.01*tail(darkCEtimeDecay, n=1)))
	tryCatch({
		previousC = coef(expfitDarkCE)["C"]
		previousD = coef(expfitDarkCE)["D"]
		expfitDarkCE <- nlrob(darkCEvoltageDecay ~ exp(C) * exp(-darkCEtimeDecay / D), start=list(C=previousC,D=previousD), data=data.frame(darkCEvoltageDecay =darkCEvoltageDecay, darkCEtimeDecay=darkCEtimeDecay))
	}, error=function(e) cat("Failed monoexponential robust fit", e$message, "\n"))

	plot(darkCEtimeDecay, darkCEvoltageDecay, log="x")
	darkCEvoltageDecay = darkCEvoltageDecay - predict(expfitDarkCE)
	lines(darkCEtimeDecay,predict(expfitDarkCE),col="red")
Sys.sleep(1)
}, error=function(e) cat("Failed monoexponential fit", e$message, "\n"))

darkLOESS = loess(darkCEvoltage~darkCEtime, span=0.002);
darkCEvoltageDecayLOESS = predict(darkLOESS)
#darkCEvoltageDecayLOESShead = head(darkCEvoltageDecayLOESS, timeEndNoiseDecayIndex)

darkCEvoltageDecayLOESSfun = approxfun(darkCEtime, darkCEvoltageDecayLOESS, method="linear", 0, 0)
#lines(darkCEtimeDecay, darkCEvoltageDecayLOESSfun(darkCEtimeDecay), col="green")

trashfornullmessages <- lapply(files, function(x) {
	message(x);
	if(mydata[[x]]$time[2] - mydata[[x]]$time[1] != deltaT){
		even_not_considering_different_time_windows_it_is_too_complex..._please_remove_this_file
	}

	startVoltage <- mean(mydata[[x]]$voltage[1:600])
	endVoltage <- mean(tail(mydata[[x]]$voltage,600))
	baseline <- seq(startVoltage, endVoltage, length.out=len)
	voltage2 <- mydata[[x]]$voltage - baseline
	
	current <- voltage2/50
	charge <- cumsum(current)*deltaT

	charge=charge-charge[match(0,mydata[[x]]$time)]
	totalcharge=mean(charge[round(length(charge)*0.9):round(length(charge)*0.95)])
	totalchargedensity=totalcharge/0.09

	voltageNoBaseline <- mydata[[x]]$voltage
#	voltagezeroDecay = tail(voltagezero, -timeZeroIndex)
#	timeDecay = tail(mydata[[x]]$time, -timeZeroIndex)
	currentNoBaseline <- voltageNoBaseline/50
	chargeNoBaseline <- cumsum(currentNoBaseline)*deltaT
	chargeNoBaseline=chargeNoBaseline-chargeNoBaseline[timeZeroIndex]

	voltageDecay = tail(voltage2, -timeStartIndex)

	expfitCE <- nlsLM(voltageDecay~ exp(C)*exp(-darkCEtimeDecay/D), start=list(C=log(max(voltageDecay)),D=0.01*tail(darkCEtimeDecay, n=1)))
tryCatch({
	expfitCE <- nlrob(voltageDecay~ exp(C)*exp(-darkCEtimeDecay/D), start=list(C=coef(expfitCE)["C"],D=coef(expfitCE)["D"]), data=data.frame(voltageDecay =voltageDecay, darkCEtimeDecay=darkCEtimeDecay))
}, error=function(e) cat("Failed monoexponential fit", e$message, "\n"))
	voltageDecayNoise = voltageDecay - predict(expfitCE)
	#for plotting debug graphs
	voltageDecayNoiseHead = head(voltageDecayNoise, timeEndNoiseDecayIndex)
	noiseLOESS = loess(voltageDecayNoise~darkCEtimeDecay, span=0.002);
	voltageDecayNoiseLOESS = predict(noiseLOESS)
	voltageDecayNoiseLOESShead = head(voltageDecayNoiseLOESS, timeEndNoiseDecayIndex)

	fitfunfun = function(cLinear, cBias, cDelay, cSlow, cSlow2, cSlow3) {
		time = cDelay + cSlow*darkCEtimeDecayHead + cSlow2*darkCEtimeDecayHead^2 + cSlow3*darkCEtimeDecayHead^3
		noise = cLinear * darkCEvoltageDecayLOESSfun(time) + cBias/time
		return(list(noise=noise, time=time))
	}
	fitfun = function(params){
		output = fitfunfun(params[1], params[2], params[3], params[4], params[5], params[6])
		differences = abs(output$noise - voltageDecayNoiseLOESShead)
		result = sum(differences)

		plot(darkCEtimeDecayHead, voltageDecayNoiseHead, col="green", pty="+")
		lines(darkCEtimeDecayHead, output$noise)
		return(result)
	}
	if(!exists("startList")){
		startList = c(1, 0, 0, 1, 0, 0);
	}
	fitNoise = optim(startList, fitfun)
	message("*****************fitNoise*****************")
	print(fitNoise)
	# convergence is good when convergence = 0
	if(fitNoise$convergence){

	fitfun0 = function(params){
		output = fitfunfun(params[1], params[2], params[3], fitNoise$par[4], fitNoise$par[5], fitNoise$par[6])
		differences = abs(output$noise - voltageDecayNoiseLOESShead)
		result = sum(differences)

		plot(darkCEtimeDecayHead, voltageDecayNoiseHead, col="green", pty="+")
		lines(darkCEtimeDecayHead, output$noise)
		return(result)
	}

		fitNoise0 = optim(head(startList,3), fitfun0)
		message("*****************fitNoise0*****************")
		print(fitNoise0)
		startList <- c(fitNoise0$par, tail(startList,3))

	fitfun1 = function(params){
		output = fitfunfun(fitNoise$par[1], fitNoise$par[2], fitNoise$par[3], params[1], params[2], params[3])
		differences = abs(output$noise - voltageDecayNoiseLOESShead)
		result = sum(differences)

		plot(darkCEtimeDecayHead, voltageDecayNoiseHead, col="green", pty="+")
		lines(darkCEtimeDecayHead, output$noise)
		return(result)
	}

		fitNoise1 = optim(tail(startList,3), fitfun1)
		message("*****************fitNoise1*****************")
		print(fitNoise1)
		startList <- c(fitNoise0$par, fitNoise1$par)
		fitNoise = optim(startList, fitfun)
		message("*****************fitNoise*****************")
		print(fitNoise)
	}
	startList <<- fitNoise$par
	message("HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH")

	fitOutput = fitfunfun(fitNoise$par[1], fitNoise$par[2], fitNoise$par[3], fitNoise$par[4], fitNoise$par[5], fitNoise$par[6]);
	noiseProfileFun = approxfun(darkCEtimeDecayHead, fitOutput$noise, method="linear", 0, 0)
	#noiseProfileExtended = fitOutput$noise
	#length(noiseProfileExtended) <- length(voltageDecay)
	#noiseProfileExtended[is.na(noiseProfileExtended)] <- 0
	noiseProfile = noiseProfileFun(darkCEtimeDecay)
	voltageMinusNoise = voltageDecay - noiseProfile

	#plot(darkCEtimeDecay, voltageDecay, xlim=c(timeStart,noiseEndTime), log="x")
	#lines(darkCEtimeDecay,predict(expfitCE),col="red")
	#lines(darkCEtimeDecay,voltageDecayNoiseLOESS,col="green")
	#lines(darkCEtime,darkCEvoltageDecayLOESS,col="blue")
	#lines(darkCEtimeDecay,noiseProfile,col="red")
	#Sys.sleep(1)

	currentMinusNoise <- voltageMinusNoise/50
	chargeMinusNoise <- cumsum(currentMinusNoise)*deltaT

	b<-strsplit(x, "_")
	c<-unlist(b)
	c2 <- c[grepl("mV",c)]
	d<-as.numeric(sub("mV.*", "", c2))
        outputChargeDensityCE <- t(c(d, totalchargedensity));
	write.table(outputChargeDensityCE, file=file.path(cedir,"outputChargeDensityCE.txt"), append=TRUE, col.names=F, row.names=F, quote=F);

	png(file.path(cedir,paste(x, ".png", sep="")), width=image_width, height=image_height)
	op <- par(mar = c(5,7,4,8.5) + 0.1) ## default is c(5,4,4,2) + 0.1
	
	xlim=c(1e-9,tail(mydata[[x]]$time,1))
	plot(mydata[[x]],type="l", ylab="", xlab="", xaxt="n", yaxt="n", xlim=xlim, log="x")
	title(ylab="Voltage (V)", cex.lab=2, line=4)
	title(xlab="Time (s)", cex.lab=2, line=3.5)
	mtext(bquote("Collected Charge Density (C/cm"^"2"*")"), cex=2, side=4,line=7,col="red")
	eaxis(side=1, cex.axis=1.5)
	eaxis(side=2, cex.axis=1.5)
	
	lines(mydata[[x]]$time, baseline, col="blue")
	lines(darkCEtimeDecay, voltageMinusNoise, col="green")

	ylim_charge=c(min(chargeNoBaseline, charge)/0.09, max(chargeNoBaseline, charge)/0.09)
	par(new=TRUE)
	plot(mydata[[x]]$time,charge/0.09, type="l", col="red", xaxt="n",yaxt="n",xlab="",ylab="", xlim=xlim, ylim=ylim_charge, log="x")
	abline(h=0,col="red")
	eaxis(4,col.ticks="red",col.axis="red", col="red", cex.axis=1.5)
	text(xlim[2]*0.75,totalchargedensity*0.9,labels=bquote(.(signif(totalchargedensity,3))~"C/cm"^"2"),cex=2,col="red")

	lines(mydata[[x]]$time,chargeNoBaseline/0.09, type="l", col="orange")
	lines(darkCEtimeDecay, chargeMinusNoise/0.09, type="l", col="green")

	graphics.off()
#reset the plotting margins
par(op)


write.table(outputChargeDensityCE, file=file.path(cedir,"outputChargeDensityCE.txt"), append=TRUE, col.names=F, row.names=F, quote=F);
})
}
