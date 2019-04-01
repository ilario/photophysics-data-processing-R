ceIntegrateExp <- function(cedir="ce")
{
  library(sfsmisc)
  library(robustbase)
  library(minpack.lm)
  library(RColorBrewer)
  
  rm("maximumDecayTime")
  
  timeStartMinimum = 3e-8
  #maximumDecayTime = 1e-4 # leave commented for setting the last time point as maximumDecayTime
  minimumDecayTime = 5e-8
  maximumDeltaV = 2
  minimumDeltaV = 1e-4
  
  options(error=function() { traceback(2); if(!interactive()) quit("no", status = 1, runLast = FALSE) })
  mycolors=brewer.pal(8,"Dark2")
  
  logdownsampling <- function(data, s=0.0003)
  {
    i <- 1
    t <- 1
    means <- c()
    length <- length(data)
    while(i < length)
    {
      by <- floor(10^(i*s))-1
      means[t] <- mean(data[i:(i+by)])
      i <- i+by+1
      t <- t+1
    }
    means=means[!is.na(means)]
    return ( 
      values = means  
    )
  }
  
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
  write.table(t(c("Voc","ChargeCE")), file=file.path(cedir,"outputChargeCE.txt"), append=FALSE, col.names=F, row.names=F);
  
  print(files[1])
  
  len<-length(mydata[[files[1]]]$voltage)
  
  
  intoRangeDV <- function(number){
    number1 = min(number, maximumDeltaV)
    number2 = max(number1, minimumDeltaV)
    return(number2)
  }
  
  trashfornullmessages <- lapply(files, function(x) {
    message(x);
    time = mydata[[x]]$time
    voltage = mydata[[x]]$voltage
    if(!exists("maximumDecayTime")){
      maximumDecayTime = tail(time,1)
    }
    deltaT = time[2] - time[1]
    
    intoRangeT <- function(number){
      number1 = min(number, maximumDecayTime)
      number2 = max(number1, minimumDecayTime)
      return(number2)
    }
    
    startVoltage <- mean(head(voltage, 600))
    endVoltage <- mean(tail(voltage, 600))
    baseline <- seq(startVoltage, endVoltage, length.out=len)
    voltage2 <- voltage - baseline
    
    current <- voltage2/50
    charge <- cumsum(current)*deltaT
    
    charge=charge-charge[match(0,time)]
    totalcharge=mean(charge[round(length(charge)*0.9):round(length(charge)*0.95)])
    totalchargedensity=totalcharge/cellArea
    
    voltageNoBaseline <- mydata[[x]]$voltage
    currentNoBaseline <- voltageNoBaseline/50
    chargeNoBaseline <- cumsum(currentNoBaseline)*deltaT
    
    chargeNoBaseline=chargeNoBaseline-chargeNoBaseline[mydata[[x]]$time == 0]
    totalchargeNoBaseline=mean(chargeNoBaseline[round(length(chargeNoBaseline)*0.9):round(length(chargeNoBaseline)*0.95)])
    totalchargedensityNoBaseline=totalchargeNoBaseline/cellArea
    
    timeStartMinimumIndex = which.min(abs(time - timeStartMinimum))
    maxVoltageAfterStartMinimumIndex = which.max(tail(voltage, -timeStartMinimumIndex))
    voltageForSearchingZero = voltage[timeStartMinimumIndex : (timeStartMinimumIndex+maxVoltageAfterStartMinimumIndex)]
    voltageChangeOfSign = sign(head(voltageForSearchingZero, -1)) + sign(tail(voltageForSearchingZero, -1))
    voltageChangeOfSignIndexes = which(voltageChangeOfSign == 0)
    voltageChangeOfSignIndexesLast = tail(voltageChangeOfSignIndexes,1)
    if(!length(voltageChangeOfSignIndexesLast)){voltageChangeOfSignIndexesLast = 0}
    timeStartIndex = timeStartMinimumIndex + voltageChangeOfSignIndexesLast
    timeStart = time[timeStartIndex]
    
    timeDecay = tail(time, -timeStartIndex)
    voltageDecay = tail(voltage2, -timeStartIndex)
    
    startCvalue = quantile(voltageDecay,0.999, names=FALSE)
    
    nlsLMsuccess = FALSE
    tryCatch({
      expfitCE1 <- nlsLM(voltageDecay ~ C*exp(-timeDecay/D), start=list(C=intoRangeDV(startCvalue),D=3e-7))
      nlsLMsuccess = TRUE
    }, error=function(e) cat("Failed first exponential fit, it is likely just noise, returning a point at (0,0)", e$message, "\n"))
    if(nlsLMsuccess){
      
      coefC = coef(expfitCE1)["C"]
      coefD = coef(expfitCE1)["D"]
      tryCatch({
        expfitCE2 <- nlrob(voltageDecay ~ C * exp(-timeDecay / D), start=list(C=intoRangeDV(coefC),D=intoRangeT(coefD)), data=data.frame(voltageDecay =voltageDecay, timeDecay=timeDecay), lower=c(C=minimumDeltaV, D=minimumDecayTime), upper=c(C=maximumDeltaV, D=maximumDecayTime), algorithm = "port")
        coefC = coef(expfitCE2)["C"]
        coefD = coef(expfitCE2)["D"]
      }, error=function(e) cat("Failed monoexponential robust fit", e$message, "\n"))
      
      voltageIntegralExp = function(x){
        coefC * coefD * (exp(-timeStart/coefD) - exp(-x/coefD))
      }
      
      #if(!exists("expfitCE2") || expfitCE2$status != "converged"){
      if(!exists("startListIntegrateExp")){
        newC = coefC/5
        newD1 = coefD*3
        newD2 = coefD/3
        startListIntegrateExp=list(C1=intoRangeDV(newC), C2=intoRangeDV(newC), D1=intoRangeT(newD1), D2=intoRangeT(newD2))
      }
      tryCatch({
        expfitCE3 <- nlrob(voltageDecay ~ C1 * exp(-timeDecay / D1) + C2 * exp(-timeDecay / D2), start=startListIntegrateExp, data=data.frame(voltageDecay =voltageDecay, timeDecay=timeDecay), lower=c(C1=minimumDeltaV, C2=minimumDeltaV, D1=minimumDecayTime, D2=minimumDecayTime), upper=c(C1=maximumDeltaV, C2=maximumDeltaV, D1=maximumDecayTime, D2=maximumDecayTime), algorithm = "port")
        coefC1 = coef(expfitCE3)["C1"]
        coefC2 = coef(expfitCE3)["C2"]
        coefD1 = coef(expfitCE3)["D1"]
        coefD2 = coef(expfitCE3)["D2"]
        #if(expfitCE3$status == "converged"){
        if(signif(coefC1,1) == minimumDeltaV || signif(coefC1,1) == maximumDeltaV || signif(coefC2,1) == minimumDeltaV || signif(coefC2,1) == maximumDeltaV || signif(coefD1,1) == minimumDecayTime || signif(coefD1,1) == maximumDecayTime || signif(coefD2,1) == minimumDecayTime || signif(coefD2,1) == maximumDecayTime){
          rm("expfitCE3")
          print("Biexponential robust fit coefficients are at the bounds, likely it failed, removing")
        } else {
          startListIntegrateExp <<- expfitCE3$coefficients
          voltageIntegralExp = function(x){
            coefC1 * coefD1 * (exp(-timeStart/coefD1) - exp(-x/coefD1)) + coefC2 * coefD2 * (exp(-timeStart/coefD2) - exp(-x/coefD2))
          }
        }
      }, error=function(e) cat("Failed biexponential robust fit", e$message, "\n"))
      #}
      if(interactive()){
        plot(timeDecay, voltageDecay, log="x")
        lines(timeDecay,predict(expfitCE1),col="red")
        if(exists("expfitCE2")){
          lines(timeDecay,predict(expfitCE2),col="green")
        }
        if(exists("expfitCE3")){
          lines(timeDecay,predict(expfitCE3),col="orange")
        }
        Sys.sleep(1)
      }
      
      chargeIntegratedExp <- voltageIntegralExp(timeDecay)/50
      totalchargeIntegratedExp = voltageIntegralExp(Inf)/50
      totalchargedensityIntegratedExp=totalchargeIntegratedExp/cellArea
      
      b<-strsplit(x, "_")
      c<-unlist(b)
      c2 <- c[grepl("mV",c)]
      d<-as.numeric(sub("mV.*", "", c2))
      
      if(output_pdf){
        pdf(file.path(cedir,paste(x, ".pdf", sep="")), width=image_midpdf_width, height=image_midpdf_height, pointsize=7)
      }else{
        png(file.path(cedir,paste(x, ".png", sep="")), width=image_width, height=image_height)
      }
      op <- par(mar = c(5,8.5,2,8.5) + 0.1) ## default is c(5,4,4,2) + 0.1
      
      xlim=c(5e-8, tail(mydata[[x]]$time,1))
      
      timeDownsampled = logdownsampling(mydata[[x]]$time)
      #plot(mydata[[x]],type="l", ylab="", xlab="", xaxt="n", xlim=xlim, cex.axis=1.4, log="x", yaxt="n", panel.first=c(lines(mydata[[x]]$time, baseline, col="gray70")))
      plot(timeDownsampled,logdownsampling(mydata[[x]]$voltage),type="l", ylab="", xlab="", xaxt="n", xlim=xlim, cex.axis=1.4, log="x", yaxt="n", panel.first=c(lines(timeDownsampled, logdownsampling(baseline), col="gray70")))
      
      title(ylab="Voltage (V)", cex.lab=1.7, line=6)
      title(xlab="Time (s)", cex.lab=1.7, line=3.5)
      mtext(bquote("Collected Charge Density (C/cm"^"2"*")"), cex=1.7, side=4,line=7,col=mycolors[3])
      eaxis(side=1, cex.axis=1.4, n.axp=1)
      yaxt <- round(axTicks(2), digits=10) # without round, the zero can get printed as 2E-17
      eaxis(side=2, cex.axis=1.4, labels = pretty10exp(yaxt))
      
      #ylim_charge=c(min(chargeNoBaseline, charge)/cellArea, max(charge, chargeIntegratedExp)/cellArea)
      
      timeDecayDownsampled = logdownsampling(timeDecay)
      if(exists("expfitCE3")){# && expfitCE3$status == "converged"){
        lines(timeDecayDownsampled, predict(expfitCE3, newdata=data.frame(timeDecay=timeDecayDownsampled)), col=mycolors[1])
        mtext(paste("Tau1 =", signif(coefD1,3),"s"), side=3, line=-10, cex=1.7, adj=0.95)
        mtext(paste("Tau2 =", signif(coefD2,3),"s"), side=3, line=-12, cex=1.7, adj=0.95)
        legendtext_fit = "Bi-Exp fitting"
      }else if (exists("expfitCE2")){
        lines(timeDecayDownsampled, predict(expfitCE2, newdata=data.frame(timeDecay=timeDecayDownsampled)), col=mycolors[1])
        mtext(paste("Tau =", signif(coefD,3),"s"), side=3, line=-10, cex=1.7, adj=0.95)
        legendtext_fit = "Exp fitting"
      } else {
        lines(timeDecayDownsampled, predict(expfitCE1, newdata=data.frame(timeDecay=timeDecayDownsampled)), col=mycolors[1])
        mtext(paste("Tau =", signif(coefD,3),"s"), side=3, line=-10, cex=1.7, adj=0.95)
        legendtext_fit = "Exp fitting"
      }	
      
      par(new=TRUE)
      plot(timeDecayDownsampled,logdownsampling(chargeIntegratedExp)/cellArea, type="l", xaxt="n",yaxt="n",xlab="",ylab="", xlim=xlim,  col=mycolors[3], log="x")#ylim=ylim_charge,
      #abline(h=0,col="red")
      eaxis(4,col.ticks=mycolors[3],col.axis=mycolors[3], col=mycolors[3], cex.axis=1.4)
      #text(xlim[2]*0.5,ylim_charge[2]*0.9,labels=bquote(.(signif(totalchargedensityIntegratedExp,3))~"C/cm"^"2"),cex=2,col="red")
      mtext(bquote(.(signif(totalchargedensityIntegratedExp,3))~"C/cm"^"2"), side=3, line=-8, cex=1.7, adj=0.95)
      
      
      #text(xlim[2]*0.75,ylim_charge[2]*0.8,labels=bquote(.(signif(totalchargedensity,3))~"C/cm"^"2"),cex=2,col="orange")
      #text(xlim[2]*0.75,ylim_charge[2]*0.7,labels=bquote(.(signif(totalchargedensityNoBaseline,3))~"C/cm"^"2"),cex=2,col="red")
      
      #lines(mydata[[x]]$time,chargeNoBaseline/cellArea, type="l", col="orange")
      #lines(mydata[[x]]$time, charge/cellArea, type="l")
      
      #legendtext = c("Signal", legendtext_fit,"Integrated charge")
      #legend(x="bottomright",inset=0,legendtext,col=c("black",mycolors[c(1,3)]), cex=1.5, lwd=4, bty="n")
      legendtext = legendtext_fit
      legend(x="bottomright",inset=0,legendtext,col=mycolors[1], cex=1.5, lwd=4, bty="n")

      graphics.off()
      #reset the plotting margins
      par(op)
      
      outputChargeDensityCE <- t(c(d, totalchargedensityIntegratedExp));
      outputChargeCE <- t(c(d, totalchargeIntegratedExp));
    }else{
      outputChargeDensityCE <- t(c(0, 0));
      outputChargeCE <- t(c(0, 0));
    }
    
    write.table(outputChargeDensityCE, file=file.path(cedir,"outputChargeDensityCE.txt"), append=TRUE, col.names=F, row.names=F, quote=F);
    write.table(outputChargeCE, file=file.path(cedir,"outputChargeCE.txt"), append=TRUE, col.names=F, row.names=F, quote=F);
  })
}
