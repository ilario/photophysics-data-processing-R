ceIntegrateExp <- function(cedir="ce")
{
library(sfsmisc)
library(robustbase)
library(minpack.lm)

timeStartMinimum = 3e-8

options(error=function() { traceback(2); if(!interactive()) quit("no", status = 1, runLast = FALSE) })

rm("expfitCE1")
rm("expfitCE2")
rm("expfitCE3")
rm("startListIntegrateExp")

print("CE: INTEGRATING")
files <- list.files(path=cedir, pattern="^CE.*\\.txt.table$");
mydata <- lapply(file.path(cedir,files), read.table, header=FALSE, col.names=c("time","voltage"));
files <- sub(".txt.table","",files);
names(mydata) <- files;
write.table(t(c("Voc","ChargeDensityCE")), file=file.path(cedir,"outputChargeDensityCE.txt"), append=FALSE, col.names=F, row.names=F);

print(files[1])

len<-length(mydata[[files[1]]]$voltage)



trashfornullmessages <- lapply(files, function(x) {
	message(x);
	time = mydata[[x]]$time
	voltage = mydata[[x]]$voltage
	deltaT = time[2] - time[1]
	startVoltage <- mean(head(voltage, 600))
	endVoltage <- mean(tail(voltage, 600))
	baseline <- seq(startVoltage, endVoltage, length.out=len)
	voltage2 <- voltage - baseline
	
	current <- voltage2/50
	charge <- cumsum(current)*deltaT

	charge=charge-charge[match(0,time)]
	totalcharge=mean(charge[round(length(charge)*0.9):round(length(charge)*0.95)])
	totalchargedensity=totalcharge/0.09

	voltageNoBaseline <- mydata[[x]]$voltage
	currentNoBaseline <- voltageNoBaseline/50
	chargeNoBaseline <- cumsum(currentNoBaseline)*deltaT

	chargeNoBaseline=chargeNoBaseline-chargeNoBaseline[mydata[[x]]$time == 0]
	totalchargeNoBaseline=mean(chargeNoBaseline[round(length(chargeNoBaseline)*0.9):round(length(chargeNoBaseline)*0.95)])
	totalchargedensityNoBaseline=totalchargeNoBaseline/0.09

timeStartMinimumIndex = which.min(abs(time - timeStartMinimum))
maxVoltageAfterStartMinimumIndex = which.max(tail(voltage, -timeStartMinimumIndex))
voltageForSearchingZero = voltage[timeStartMinimumIndex : (timeStartMinimumIndex+maxVoltageAfterStartMinimumIndex)]
voltageChangeOfSign = sign(head(voltageForSearchingZero, -1)) + sign(tail(voltageForSearchingZero, -1))
voltageChangeOfSignIndexes = which(voltageChangeOfSign == 0)
voltageChangeOfSignIndexesLast = tail(voltageChangeOfSignIndexes,1)
timeStartIndex = timeStartMinimumIndex + voltageChangeOfSignIndexesLast
timeStart = time[timeStartIndex]

	timeDecay = tail(time, -timeStartIndex)
	voltageDecay = tail(voltage2, -timeStartIndex)

	startCvalue = log(quantile(voltageDecay,0.999, names=FALSE))
	expfitCE1 <- nlsLM(voltageDecay ~ exp(C)*exp(-timeDecay/D), start=list(C=startCvalue,D=3e-7))

	#tryCatch({
		previousC = coef(expfitCE1)["C"]
		previousD = coef(expfitCE1)["D"]
		expfitCE2 <- nlrob(voltageDecay ~ exp(C) * exp(-timeDecay / D), start=list(C=previousC,D=previousD), data=data.frame(voltageDecay =voltageDecay, timeDecay=timeDecay))#, method="MM")#, control=nlrob.control("MM"))
	#}, error=function(e) cat("Failed monoexponential robust fit", e$message, "\n"))
	coefC = coef(expfitCE2)["C"]
	coefD = coef(expfitCE2)["D"]
	voltageIntegralExp = function(x){
		exp(coefC) * coefD * (exp(-timeStart/coefD) - exp(-x/coefD))
	}

	if(expfitCE2$status != "converged"){
		if(!exists("startListIntegrateExp")){
			newC = coefC/5
			newD1 = asin(sqrt((coefD*3 - 5e-8)*1e5))
			newD2 = asin(sqrt((coefD/3 - 5e-8)*1e5))
			startListIntegrateExp=list(C1=newC, C2=newC, D1=newD1, D2=newD2)
		}
		tryCatch({
			expfitCE3 <- nlrob(voltageDecay ~ exp(C1) * exp(-timeDecay / (5e-8+1e-5*sin(D1)^2)) + exp(C2) * exp(-timeDecay / (5e-8+1e-5*sin(D2)^2)), start=startListIntegrateExp, data=data.frame(voltageDecay =voltageDecay, timeDecay=timeDecay))
			coefC1 = coef(expfitCE3)["C1"]
			coefC2 = coef(expfitCE3)["C2"]
			coefD1 = coef(expfitCE3)["D1"]
			coefD2 = coef(expfitCE3)["D2"]
			coefD1complete = 5e-8+1e-5*sin(coefD1)^2
			coefD2complete = 5e-8+1e-5*sin(coefD2)^2
			startListIntegrateExp <<- expfitCE3$coefficients
			#if(expfitCE3$status == "converged"){
			voltageIntegralExp = function(x){
				exp(coefC1) * coefD1complete * (exp(-timeStart/coefD1complete) - exp(-x/coefD1complete)) + exp(coefC2) * coefD2complete * (exp(-timeStart/coefD2complete) - exp(-x/coefD2complete))
			#}
			}
		}, error=function(e) cat("Failed biexponential robust fit", e$message, "\n"))
	}
	if(interactive()){
		plot(timeDecay, voltageDecay, log="x")
		lines(timeDecay,predict(expfitCE1),col="red")
		lines(timeDecay,predict(expfitCE2),col="orange")
		if(exists("expfitCE3")){
			lines(timeDecay,predict(expfitCE3),col="green")
		}
		Sys.sleep(3)
	}

	chargeIntegratedExp <- voltageIntegralExp(timeDecay)/50
	totalchargeIntegratedExp = voltageIntegralExp(Inf)/50
	totalchargedensityIntegratedExp=totalchargeIntegratedExp/0.09


	b<-strsplit(x, "_")
	c<-unlist(b)
	c2 <- c[grepl("mV",c)]
	d<-as.numeric(sub("mV.*", "", c2))
        outputChargeDensityCE <- t(c(d, totalchargedensityIntegratedExp));
	write.table(outputChargeDensityCE, file=file.path(cedir,"outputChargeDensityCE.txt"), append=TRUE, col.names=F, row.names=F, quote=F);

	png(file.path(cedir,paste(x, ".png", sep="")), width=image_width, height=image_height)
	op <- par(mar = c(5,7,4,8.5) + 0.1) ## default is c(5,4,4,2) + 0.1
	
	xlim=c(timeStart, tail(mydata[[x]]$time,1))
	plot(mydata[[x]],type="l", ylab="", xlab="", xaxt="n", xlim=xlim, cex.axis=1.5, log="x")#, yaxt="n"
	title(ylab="Voltage (V)", cex.lab=2, line=4)
	title(xlab="Time (s)", cex.lab=2, line=3.5)
	mtext(bquote("Collected Charge Density (C/cm"^"2"*")"), cex=2, side=4,line=7,col="red")
	eaxis(side=1, cex.axis=1.5)
	#eaxis(side=2, cex.axis=1.5)
	
	lines(mydata[[x]]$time, baseline, col="blue")
	if(exists("expfitCE3")){# && expfitCE3$status == "converged"){
		lines(timeDecay, predict(expfitCE3, newdata=data.frame(timeDecay=timeDecay)), col="orange")
	}else{
		lines(timeDecay, predict(expfitCE2, newdata=data.frame(timeDecay=timeDecay)), col="green")
	}
		
	ylim_charge=c(min(chargeNoBaseline, charge)/0.09, max(charge, chargeIntegratedExp)/0.09)
	par(new=TRUE)
	plot(timeDecay,chargeIntegratedExp/0.09, type="l", xaxt="n",yaxt="n",xlab="",ylab="", xlim=xlim, ylim=ylim_charge, col="red", log="x")
	#abline(h=0,col="red")
	eaxis(4,col.ticks="red",col.axis="red", col="red", cex.axis=1.5)
	text(xlim[2]*0.5,ylim_charge[2]*0.9,labels=bquote(.(signif(totalchargedensityIntegratedExp,3))~"C/cm"^"2"),cex=2,col="red")
	#text(xlim[2]*0.75,ylim_charge[2]*0.8,labels=bquote(.(signif(totalchargedensity,3))~"C/cm"^"2"),cex=2,col="orange")
	#text(xlim[2]*0.75,ylim_charge[2]*0.7,labels=bquote(.(signif(totalchargedensityNoBaseline,3))~"C/cm"^"2"),cex=2,col="red")

	#lines(mydata[[x]]$time,chargeNoBaseline/0.09, type="l", col="orange")
	#lines(mydata[[x]]$time, charge/0.09, type="l")

	graphics.off()
#reset the plotting margins
par(op)


write.table(outputChargeDensityCE, file=file.path(cedir,"outputChargeDensityCE.txt"), append=TRUE, col.names=F, row.names=F, quote=F);
})
}
