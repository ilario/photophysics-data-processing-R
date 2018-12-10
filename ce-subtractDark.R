ceSubtractDark <- function(cedir="ce")
{
library(sfsmisc)
library(robustbase)
library(minpack.lm)

timeStartMinimum = 3e-8
noiseEndTime = 4e-7

rm("startList")
options(error=function() { traceback(2); if(!interactive()) quit("no", status = 1, runLast = FALSE) })

print("CE: INTEGRATING")
files <- list.files(path=cedir, pattern="^CE.*\\.txt.table$");
mydata <- lapply(file.path(cedir,files), read.table, header=FALSE, col.names=c("time","voltage"));
files <- sub(".txt.table","",files);
names(mydata) <- files;
write.table(t(c("Voc","ChargeDensityCE")), file=file.path(cedir,"outputChargeDensityCE.txt"), append=FALSE, col.names=F, row.names=F);

print(files[1])
darkCEvoltageNoBaseline = mydata[[files[1]]]$voltage;
len<-length(darkCEvoltageNoBaseline)
darkStartVoltage <- mean(head(darkCEvoltageNoBaseline,600))
darkEndVoltage <- mean(tail(darkCEvoltageNoBaseline,600))
darkBaseline <- seq(darkStartVoltage, darkEndVoltage, length.out=len)
darkCEvoltage <- darkCEvoltageNoBaseline - darkBaseline
darkCEtime = mydata[[files[1]]]$time;

timeStartMinimumIndex = which.min(abs(darkCEtime - timeStartMinimum))
maxDarkVoltageAfterStartMinimumIndex = which.max(tail(darkCEvoltage, -timeStartMinimumIndex))
timeDarkPeak = darkCEtime[timeStartMinimumIndex + maxDarkVoltageAfterStartMinimumIndex]
timeStart = darkCEtime[timeStartMinimumIndex + maxDarkVoltageAfterStartMinimumIndex] * 0.9

timeDarkStartIndex = which.min(abs(darkCEtime-timeStart))
darkCEtimeDecay = tail(darkCEtime, -timeDarkStartIndex)
darkCEvoltageDecay = tail(darkCEvoltage, -timeDarkStartIndex)

timeDarkEndNoiseDecayIndex = which.min(abs(darkCEtimeDecay-noiseEndTime))
darkCEtimeDecayHead = head(darkCEtimeDecay, timeDarkEndNoiseDecayIndex)

tryCatch({
	#plot(darkCEtimeDecay, darkCEvoltageDecay)
	#lines(darkCEtimeDecay, exp(log(quantile(darkCEvoltageDecay,0.999)))*exp(-darkCEtimeDecay/2e-7), col="green")
	expfitDarkCE <- nlsLM(darkCEvoltageDecay ~ exp(C)*exp(-darkCEtimeDecay/D), start=list(C=log(quantile(darkCEvoltageDecay,0.999)),D=2e-7))
	#lines(darkCEtimeDecay, predict(expfitDarkCE), col="red")
	#Sys.sleep(2)
	tryCatch({
		expfitDarkCE <- nlrob(darkCEvoltageDecay ~ exp(C) * exp(-darkCEtimeDecay / D), start=list(C=log(quantile(darkCEvoltageDecay,0.999)),D=2e-7), data=data.frame(darkCEvoltageDecay =darkCEvoltageDecay, darkCEtimeDecay=darkCEtimeDecay))
	#lines(darkCEtimeDecay, predict(expfitDarkCE), col="orange")
	#Sys.sleep(2)
}, error=function(e) cat("Failed monoexponential robust fit", e$message, "\n"))

	plot(darkCEtimeDecay, darkCEvoltageDecay, log="x")
	darkCEvoltageDecay = darkCEvoltageDecay - predict(expfitDarkCE)
	lines(darkCEtimeDecay,predict(expfitDarkCE),col="red")
Sys.sleep(1)
}, error=function(e) cat("Failed monoexponential fit", e$message, "\n"))

darkLOESS = loess(darkCEvoltageDecay~darkCEtimeDecay, span=0.001);
darkCEvoltageDecayLOESS = predict(darkLOESS)

darkCEvoltageDecayLOESSfun = approxfun(darkCEtimeDecay, darkCEvoltageDecayLOESS, method="linear", 0, 0)
lines(darkCEtimeDecay, darkCEvoltageDecayLOESSfun(darkCEtimeDecay), col="green")

trashfornullmessages <- lapply(files, function(x) {
	message(x);
	deltaT = mydata[[x]]$time[2] - mydata[[x]]$time[1]
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
	currentNoBaseline <- voltageNoBaseline/50
	chargeNoBaseline <- cumsum(currentNoBaseline)*deltaT

	chargeNoBaseline=chargeNoBaseline-chargeNoBaseline[mydata[[x]]$time == 0]
	totalchargeNoBaseline=mean(chargeNoBaseline[round(length(chargeNoBaseline)*0.9):round(length(chargeNoBaseline)*0.95)])
	totalchargedensityNoBaseline=totalchargeNoBaseline/0.09

	timeStartIndex = which.min(abs(mydata[[x]]$time-timeStart))
	timeDecay = tail(mydata[[x]]$time, -timeStartIndex)
	voltageDecay = tail(voltage2, -timeStartIndex)

	expfitCE <- nlsLM(voltageDecay~ exp(C)*exp(-timeDecay/D), start=list(C=log(max(voltageDecay)),D=0.01*tail(timeDecay, n=1)))
tryCatch({
	expfitCE <- nlrob(voltageDecay~ exp(C)*exp(-darkCEtimeDecay/D), start=list(C=coef(expfitCE)["C"],D=coef(expfitCE)["D"]), data=data.frame(voltageDecay =voltageDecay, darkCEtimeDecay=darkCEtimeDecay))
}, error=function(e) cat("Failed monoexponential fit", e$message, "\n"))

	timeEndNoiseDecayIndex = which.min(abs(timeDecay-noiseEndTime))
	timeDecayHead = head(timeDecay, timeEndNoiseDecayIndex)
	voltageDecayHead = head(voltageDecay, timeEndNoiseDecayIndex)

	decayLOESS = loess(voltageDecay~timeDecay, span=0.001);
	voltageDecayLOESS = predict(decayLOESS)
	voltageDecayLOESShead = head(voltageDecayLOESS, timeEndNoiseDecayIndex)

	maxTime = timeDecayHead[which.max(voltageDecayHead)]
	timeShift = timeDarkPeak - maxTime

	fitfunfun = function(cLinear, cBias, cDecay, cSlow, cSlow2, cNoiseDecay) {
		timeCenteredOnPeak = timeDecayHead - maxTime
		time = timeShift + maxTime + (1+0.5*sin(cSlow)^2)*timeCenteredOnPeak + cSlow2*timeCenteredOnPeak^2
		timeTransformedDarkNoise = darkCEvoltageDecayLOESSfun(time)
		noise = cLinear * timeTransformedDarkNoise + cNoiseDecay * time * timeTransformedDarkNoise + exp(cBias) * exp(-time/cDecay)
		nonZeros = as.logical(timeTransformedDarkNoise)
		return(list(noise=noise, time=time, nonZeros=nonZeros))
	}
	fitfun = function(params){
		output = fitfunfun(params[1], params[2], params[3], params[4], params[5], params[6])
		differences = (output$noise[output$nonZeros] - voltageDecayLOESShead[output$nonZeros])^2
		result = sum(differences)

		#plot(timeDecayHead, voltageDecayHead, col="green", pch="+")
		#lines(timeDecayHead, output$noise)
		#if(first){
		#	Sys.sleep(3)
		#	first <<- FALSE
		#}
		return(result)
	}
	if(!exists("startList")){
		startList = c(1, coef(expfitCE)["C"], coef(expfitCE)["D"], pi/2, 0, 0);
	}
	first <<- TRUE
	startList[2] = coef(expfitCE)["C"]
	startList[3] = coef(expfitCE)["D"]

	fitNoise = optim(startList, fitfun)#lower = c(0, -Inf, 0, 0, -Inf, -Inf), upper = c(10, 10, 1, 10, Inf, Inf))
	message("*****************fitNoise*****************")
	print(fitNoise)
	# convergence is good when convergence = 0
	if(fitNoise$convergence){

	fitfun0 = function(params){
		output = fitfunfun(params[1], params[2], params[3], fitNoise$par[4], fitNoise$par[5], fitNoise$par[6])
		differences = (output$noise[output$nonZeros] - voltageDecayLOESShead[output$nonZeros])^2
		result = sum(differences)

		#plot(timeDecayHead, voltageDecayHead, col="red", pch="+")
		#lines(timeDecayHead, output$noise)
		return(result)
	}

		fitNoise0 = optim(head(startList,3), fitfun0)
		message("*****************fitNoise0*****************")
		print(fitNoise0)
		#startList <- c(fitNoise0$par, c(startList[4], startList[], startList[]))

	fitfun1 = function(params){
		output = fitfunfun(fitNoise$par[1], fitNoise$par[2], fitNoise$par[3], params[1], params[2], fitNoise$par[6])
		differences = (output$noise[output$nonZeros] - voltageDecayLOESShead[output$nonZeros])^2
		result = sum(differences)

		#plot(timeDecayHead, voltageDecayHead, col="orange", pch="+")
		#lines(timeDecayHead, output$noise)
		return(result)
	}

		fitNoise1 = optim(c(0,0), fitfun1)#c(fitNoise$par[4], fitNoise$par[5]), fitfun1)
		message("*****************fitNoise1*****************")
		print(fitNoise1)
		#startList <- c(fitNoise0$par, fitNoise1$par, )
	fitfun2 = function(params){
		output = fitfunfun(params[1], params[2], params[3], fitNoise$par[4], fitNoise$par[5], params[4])
		differences = (output$noise[output$nonZeros] - voltageDecayLOESShead[output$nonZeros])^2
		result = sum(differences)

		#plot(timeDecayHead, voltageDecayHead, col="orange", pch="+")
		#lines(timeDecayHead, output$noise)
		return(result)
	}

		fitNoise2 = optim(c(fitNoise0$par[1], fitNoise0$par[2], fitNoise0$par[3], fitNoise$par[6]), fitfun2)
		message("*****************fitNoise2*****************")
		print(fitNoise2)
		startList <- c(fitNoise2$par[1], fitNoise2$par[2], fitNoise2$par[3], fitNoise1$par[1], fitNoise1$par[2], fitNoise2$par[4])
		fitNoise = optim(startList, fitfun)
		message("*****************fitNoise*****************")
		print(fitNoise)
	}
	startList <<- fitNoise$par
	message("HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH")

	fitOutput = fitfunfun(fitNoise$par[1], -Inf, Inf, fitNoise$par[4], fitNoise$par[5], fitNoise$par[6]);
	noiseProfileFun = approxfun(timeDecayHead, fitOutput$noise, method="linear", 0, 0)

	noiseProfile = noiseProfileFun(timeDecay)
	voltageMinusNoise = voltageDecay - noiseProfile

	#just for plotting
	fitOutput2 = fitfunfun(fitNoise$par[1], fitNoise$par[2], fitNoise$par[3], fitNoise$par[4], fitNoise$par[5], fitNoise$par[6]);
	noiseProfileFun2 = approxfun(timeDecayHead, fitOutput2$noise, method="linear", 0, 0)
	noiseProfile2 = noiseProfileFun2(timeDecay)
	#plot(timeDecay, noiseProfile, xlim=c(timeStart,noiseEndTime), log="x",col="red", type="l")
	#lines(timeDecay,voltageDecay, pch="+")
	#Sys.sleep(1)

	currentMinusNoise <- voltageMinusNoise/50
	chargeMinusNoise <- cumsum(currentMinusNoise)*deltaT
	totalchargeMinusNoise=mean(chargeMinusNoise[round(length(chargeMinusNoise)*0.9):round(length(chargeMinusNoise)*0.95)])
	totalchargedensityMinusNoise=totalchargeMinusNoise/0.09


	b<-strsplit(x, "_")
	c<-unlist(b)
	c2 <- c[grepl("mV",c)]
	d<-as.numeric(sub("mV.*", "", c2))
        outputChargeDensityCE <- t(c(d, totalchargedensityMinusNoise));
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
	lines(timeDecay, voltageMinusNoise, col="green")
	lines(timeDecay, noiseProfile2, col="red")
	lines(darkCEtimeDecay, darkCEvoltageDecay, col="orange")

	ylim_charge=c(min(chargeNoBaseline, charge)/0.09, max(chargeMinusNoise, charge)/0.09)
	par(new=TRUE)
	plot(timeDecay,chargeMinusNoise/0.09, type="l", col="green", xaxt="n",yaxt="n",xlab="",ylab="", xlim=xlim, ylim=ylim_charge, log="x")
	abline(h=0,col="red")
	eaxis(4,col.ticks="red",col.axis="red", col="red", cex.axis=1.5)
	text(xlim[2]*0.75,ylim_charge[2]*0.9,labels=bquote(.(signif(totalchargedensityMinusNoise,3))~"C/cm"^"2"),cex=2,col="green")
	text(xlim[2]*0.75,ylim_charge[2]*0.8,labels=bquote(.(signif(totalchargedensity,3))~"C/cm"^"2"),cex=2,col="orange")
	text(xlim[2]*0.75,ylim_charge[2]*0.7,labels=bquote(.(signif(totalchargedensityNoBaseline,3))~"C/cm"^"2"),cex=2,col="red")

	lines(mydata[[x]]$time,chargeNoBaseline/0.09, type="l", col="orange")
	lines(mydata[[x]]$time, charge/0.09, type="l", col="red")

	graphics.off()
#reset the plotting margins
par(op)


write.table(outputChargeDensityCE, file=file.path(cedir,"outputChargeDensityCE.txt"), append=TRUE, col.names=F, row.names=F, quote=F);
})
}