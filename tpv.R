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

tpv <- function(tpvdir="tpv")
{
  print("TPV: FITTING")
  library(robustbase)
  library(sfsmisc)
  library(RColorBrewer)
  library(gplots)

  doBiexpFit=T
  robust=T
  doPlots=T
  plotHist2d=T
  logy=F
  logx=F
  residuals=F
  thresholdBiexp=10
  thresholdRobustBiexp=100
  forbidNegativeDecays=T
  noiseTime=5e-8
  #debugDeltaV info will appear in biexp logx graphics
  debugDeltaV=F
  
  if(doPlots){
    if(!exists("output_pdf")){stop("images width and height variables must be set, via limits_for_graphics.R")}
    mycolors=brewer.pal(8,"Dark2")
    mycolors=mycolors[c(2,4,6,7)] # remove green and blue cause they're not visible over YlGnBu
    mycolors2_temp=brewer.pal(8,"YlGnBu")
    mycolors2_get = colorRampPalette(mycolors2_temp[3:8])
    #mycolors2 = c("#ffffff", mycolors2_get(15))# white background
    mycolors2 = mycolors2_get(16)
  }
  
  options(error=function() { traceback(2); if(!interactive()) quit("no", status = 1, runLast = FALSE) })
  
  files <- list.files(path=tpvdir, pattern="^TPV.*\\.txt.table$");
  mydata <- lapply(file.path(tpvdir,files), read.table, header=FALSE, col.names=c("time","voltage"));
  files <- sub(".txt.table","",files);
  names(mydata) <- files;
  ## output for importing
  write.table(t(c("file","Voc","A","T","T.error")), file=file.path(tpvdir,"output-monoexp.txt"), append=FALSE, col.names=F, row.names=F);
  write.table(t(c("file","Voc","deltaV")), file=file.path(tpvdir,"outputDeltaV.txt"), append=FALSE, col.names=F, row.names=F);
  write.table(t(c("file","Voc","deltaV")), file=file.path(tpvdir,"outputDeltaVloess.txt"), append=FALSE, col.names=F, row.names=F);
  write.table(t(c("file","Voc","deltaV")), file=file.path(tpvdir,"outputDeltaVfirstPoints.txt"), append=FALSE, col.names=F, row.names=F);
  write.table(t(c("file","Voc","deltaV")), file=file.path(tpvdir,"outputDeltaVmonoexp.txt"), append=FALSE, col.names=F, row.names=F);
  if(robust){
    write.table(t(c("file","Voc","A","T","T.error")), file=file.path(tpvdir,"output-robustmonoexp.txt"), append=FALSE, col.names=F, row.names=F);
  }
  if(doBiexpFit){
    write.table(t(c("file","Voc","A1","T1","T1.error","A2","T2","T2.error")), file=file.path(tpvdir,"output-biexp.txt"), append=FALSE, col.names=F, row.names=F);
    write.table(t(c("file","Voc","deltaV")), file=file.path(tpvdir,"outputDeltaVbiexp.txt"), append=FALSE, col.names=F, row.names=F);
    if(robust){
      write.table(t(c("file","Voc","A1","T1","T1.error","A2","T2","T2.error")), file=file.path(tpvdir,"output-robustbiexp.txt"), append=FALSE, col.names=F, row.names=F);
    }
  }
  
  trashfornullmessages <- lapply(files, function(x) {
    message(x);	
    #header =  read.table(file.path(tpvdir,paste(x,".txt",sep="")), skip=3, header=FALSE, nrows=3)$V2
    
    peaktime <- mydata[[x]]$time[which.max(mydata[[x]]$voltage)];
    starttime <- head(mydata[[x]]$time, n=1)
    voltage2 <- subset(mydata[[x]], time >= peaktime, select=voltage); 
    time2 <- subset(mydata[[x]], time >= peaktime, select=time);
    temp <- data.frame(time2, voltage2);
    voltage_nonfit <- subset(mydata[[x]], time < peaktime, select=voltage); 
    time_nonfit <- subset(mydata[[x]], time < peaktime, select=time);
    temp_nonfit <- data.frame(time_nonfit, voltage_nonfit);
    endtime = tail(temp$time, n=1)
    startingvoltage <- mean(mydata[[x]]$voltage[0:200]);
    maxVoltageIndex <- which.max(mydata[[x]]$voltage)
    deltavoltage <- mydata[[x]]$voltage[maxVoltageIndex] - startingvoltage;
    outputDeltaV <- t(c(x, startingvoltage, deltavoltage));
    write.table(outputDeltaV, file=file.path(tpvdir,"outputDeltaV.txt"), append=TRUE, col.names=F, row.names=F, quote=F);
    
    expectedresult <- 0.00002;
    A=startingvoltage; B=deltavoltage; 
    C=expectedresult; slowdecay <- C; fastdecay <- C; #in caso non venisse effettuato il monoexp esplorativo
    tempsubset <- subset(mydata[[x]], time >= peaktime & time < endtime/20, select=c(time,voltage));
    tempsubset2 <- subset(mydata[[x]], time > endtime/120, select=c(time,voltage));
    lo2 <- loess(tempsubset2$voltage~tempsubset2$time, span=0.02);
    deltaVloess <- max(predict(lo2))-startingvoltage
    outputDeltaVloess <- t(c(x, startingvoltage, deltaVloess));
    write.table(outputDeltaVloess, file=file.path(tpvdir,"outputDeltaVloess.txt"), append=TRUE, col.names=F, row.names=F, quote=F);
    
    afterPeakNoiseTimePoints = time2 < peaktime + noiseTime
    afterPeakNoiseVoltagePoints = voltage2[afterPeakNoiseTimePoints]
    print(paste("DeltaVfirstPoints averaging on", length(afterPeakNoiseVoltagePoints), "points"))
    deltaVfirstPoints <- mean(afterPeakNoiseVoltagePoints) - startingvoltage
    outputDeltaVfirstPoints <- t(c(x, startingvoltage, deltaVfirstPoints));
    write.table(outputDeltaVfirstPoints, file=file.path(tpvdir,"outputDeltaVfirstPoints.txt"), append=TRUE, col.names=F, row.names=F, quote=F);
    
    tryCatch({
      print("EvaluationMonoexp/Fit: Performing");
      fit <- nls(voltage ~ cbind(1, exp(-time/C)), start=list(C=expectedresult*runif(1,0.5,2)),trace=F,data=temp,alg="plinear");
      capture.output(summary(fit), file=file.path(tpvdir,paste(x, "-fit", sep="")),  append=TRUE);
      C <- coef(fit)["C"];
      C.err <- summary(fit)$coefficients["C",2];
      slowdecay <- C; fastdecay <- C; #in caso non venisse effettuato il biexp
    }, error=function(e) cat("EvaluationMonoexp: Error ", e$message, "\n"));
    
    
    
    if(doBiexpFit){
      biexpsuccess = 0
      print("Biexp/Fit: Performing")
      step <- 0.5;
      removedPoints <- 0
      while(!biexpsuccess & removedPoints < thresholdBiexp){
        for(rand in (step*seq(1, 10, by=step))){
          if(forbidNegativeDecays){while(mean(temp[1:10,2]) < mean(temp[10:20,2])) {print("Biexp/Fit: Removing a point for clearing the curve - rough"); temp <- temp[-1,]}
            while(temp[1,2] < mean(temp[2:4,2]) | temp[1,2] < mean(temp[10:15,2])) {print("Biexp/Fit: Removing a point for clearing the curve - fine"); temp <- temp[-1,]}
            while(temp[1,2] < mean(temp[5:10,2]) + 0.003*deltavoltage) {print("Biexp/Fit: Removing a point for clearing the curve - flat part"); temp <- temp[-1,]}}
          tryCatch({
            if(exists("CFromPreviousDecay") && rand == step){
              print("Biexp/Fit: Trying with previous parameters")
              Crand <- CFromPreviousDecay*(2^runif(1,-step,step))
              Frand <- FFromPreviousDecay*(2^runif(1,-step,step))
            } else {
              Crand <- C*(3^runif(1,-rand,rand))*2;
              Frand <- max(endtime/500*(3^runif(1,-rand,rand)), endtime/1500)#C*(3^runif(1,-rand,rand))/100;
            }
            fit2 <- nls(voltage ~ cbind(1, exp(-time/C), exp(-time/F)), start=list(C=Crand, F=Frand),trace=F,data=temp,alg="plinear");
            CFromPreviousDecay <<- coef(fit2)["C"] 
            FFromPreviousDecay <<- coef(fit2)["F"]
            A2 <- coef(fit2)[".lin1"]; A2.err <- summary(fit2)$coefficients[".lin1",2];
            #print(paste(startingvoltage, "=>", coef(fit2)[".lin1"], "/// ... =>", coef(fit2)[".lin2"], "///", Crand, "=>", coef(fit2)["C"], "/// ... =>", coef(fit2)[".lin3"], "--", Frand, "=>", coef(fit2)["F"]));
            if(coef(fit2)["C"] > coef(fit2)["F"]) {slowampl <- coef(fit2)[".lin2"]; slowdecay <- coef(fit2)["C"]; fastampl <- coef(fit2)[".lin3"]; fastdecay <- coef(fit2)["F"]; slowampl.err <- summary(fit2)$coefficients[".lin2",2]; slowdecay.err <- summary(fit2)$coefficients["C",2]; fastampl.err <- summary(fit2)$coefficients[".lin3",2]; fastdecay.err <- summary(fit2)$coefficients["F",2]} else {fastampl <- coef(fit2)[".lin2"]; fastdecay <- coef(fit2)["C"]; slowampl <- coef(fit2)[".lin3"]; slowdecay <- coef(fit2)["F"]; fastampl.err <- summary(fit2)$coefficients[".lin2",2]; fastdecay.err <- summary(fit2)$coefficients["C",2]; slowampl.err <- summary(fit2)$coefficients[".lin3",2]; slowdecay.err <- summary(fit2)$coefficients["F",2]};
            biexpNext=0
            if(fastdecay < endtime/1500){print("Biexp/Fit: Bad result, a decay is too fast"); biexpNext=1}
            if(slowampl < 0){print("Biexp/Fit: Bad result, slow decay is negative"); biexpNext=1}
            if(biexpNext){next}
            if(forbidNegativeDecays){if(coef(fit2)[".lin2"] < 0 | coef(fit2)[".lin3"] < 0){print("Biexp/Fit: Bad result"); next}};
            print(paste("C starting: ", Crand, ", C fitted: ", CFromPreviousDecay))
            print(paste("F starting: ", Frand, ", F fitted: ", FFromPreviousDecay))
            capture.output(summary(fit2), file=file.path(tpvdir,paste(x, "-fit", sep="")),  append=TRUE);
            
            if(doPlots){
              print("Biexp/Plot/Linear: Performing");
              tryCatch({
                png(file.path(tpvdir,paste(x, "-biexp.png", sep="")), width = image_width, height = image_height);
                plot(mydata[[x]]$time, mydata[[x]]$voltage, xlab = "Time (s)", ylab = "Voltage (V)", pch=".", col="yellow");
                points(temp, pch=".");
                lines(tempsubset2$time, predict(lo2), col="black", lwd=3);
                lines(temp$time, predict(fit2), col="aquamarine4", lwd=2);
                segments(starttime, A2, 0, A2, col="aquamarine4")
                # TO BE FIXED: REMOVE HEADER VARIABLE
                #lines(temp$time+header[3]/10*header[1], A2 + slowampl*exp(-temp$time/slowdecay)+deltavoltage/10, col="green");
                mtext(paste("Tau1 =", signif(slowdecay,digits=4), "s"), side=3, line=-5, adj=NA, col="green", cex=2);
                # TO BE FIXED: REMOVE HEADER VARIABLE
                # lines(temp$time+header[3]/10*header[1], A2 + fastampl*exp(-temp$time/fastdecay)+deltavoltage/10, col="blue");
                mtext(paste("Tau2 =", signif(fastdecay,digits=4), "s"), side=3, line=-7, adj=NA, col="blue", cex=2);
                graphics.off();
              }, error=function(e) cat("Biexp/Plot/Linear: Error ", e$message, "\n"));
            }
            if(logy & doPlots){			print("Biexp/Plot/Log: Performing");
              tryCatch({
                temp_logy_biexp = subset(temp, time <= slowdecay*5)
                png(file.path(tpvdir,paste(x, "-biexp-log.png", sep="")), width = image_width, height = image_height);
                plot(temp_logy_biexp$time, temp_logy_biexp$voltage - A2, xlab = "Time (s)", ylab = paste("Log(Voltage (V) -", A2, "V)"), pch=".", log="y", ylim=c(max(mean(tail(temp_logy_biexp$voltage, 10))-A2, deltavoltage/1e5), deltavoltage));
                #			points(tempsubset$time, tempsubset$voltage - A2, pch=".", col="yellow");
                lines(tempsubset2$time, predict(lo2) - A2, col='black', lwd=2);
                lines(temp$time,predict(fit2) - A2, col="aquamarine4", lwd=3);
                mtext(paste("T1 =", signif(slowdecay,digits=4), "s"), side=3, line=-5, adj=NA, col="green", cex=2);
                mtext(paste("T2 =", signif(fastdecay,digits=4), "s"), side=3, line=-7, adj=NA, col="blue", cex=2);
                graphics.off();	
              }, error=function(e) cat("Biexp/Plot/Log: Error ", e$message, "\n"));
            }
            if(logx & doPlots){
              print("Biexp/Plot/LogX: Performing");
              tryCatch({
                
                png(file.path(tpvdir,paste(x, "-biexp-logx.png", sep="")), width = image_width, height = image_height);
                plot(temp$time, temp$voltage, xlab = "Log(Time (s))", ylab = "Voltage (V)", pch=".", log="x", xlim=c(max(temp$time[1],-temp$time[1]),endtime));
                #lines(temp$time, predict(lo2), col='black', lwd=3);
                lines(temp$time,predict(fit2), col="aquamarine4", lwd=2);
                mtext(paste("Tau1 =", signif(slowdecay,digits=4), "\u00b1", signif(slowdecay.err,digits=4), "s"), side=3, line=-5, adj=NA, col="green", cex=2);
                mtext(paste("Tau2 =", signif(fastdecay,digits=4), "\u00b1", signif(fastdecay.err,digits=4), "s"), side=3, line=-7, adj=NA, col="blue", cex=2);
                
                if(debugDeltaV){
                  abline(h=startingvoltage + deltaVfirstPoints, lwd=5)
                  abline(h=startingvoltage + deltaVloess, col="red")
                  abline(h=startingvoltage + deltavoltage, col="green")
                }
                
                graphics.off();	
              }, error=function(e) cat("Biexp/Plot/LogX: Error ", e$message, "\n"));
            }
            if(residuals & doPlots){
              tryCatch({
                print("Biexp/Plot/Residuals: Performing");
                png(file.path(tpvdir,paste(x, "-biexp-residuals.png", sep="")), width = image_width, height = image_height);
                yresidualfit2 <- temp$voltage-predict(fit2,newdata=data.frame(time=temp$time))
                plot(temp$time, yresidualfit2, ylim=c(quantile(yresidualfit2,0.001),quantile(yresidualfit2,0.999)), xlab = "Time (s)", ylab = "Voltage (V)", pch=".", log="x", col=rgb(0,0,0, (((endtime*8)-temp$time)/(endtime*8))^10));
                abline(h=0, col="red");
                mtext(paste("Tau1 =", signif(slowdecay,digits=4), "\u00b1", signif(slowdecay.err,digits=4), "s"), side=3, line=-5, adj=NA, col="green", cex=2);
                mtext(paste("Tau2 =", signif(fastdecay,digits=4), "\u00b1", signif(fastdecay.err,digits=4), "s"), side=3, line=-7, adj=NA, col="blue", cex=2);
                graphics.off();
              }, error=function(e) cat("Biexp/Plot/Residuals: Error ", e$message, "\n"));
            }
            ## output for importing
            output <- t(c(x, A2, slowampl, slowdecay, slowdecay.err, fastampl, fastdecay, fastdecay.err));
            write.table(output, file=file.path(tpvdir,"output-biexp.txt"), append=TRUE, col.names=F, row.names=F, quote=F);
            
            #		outputDeltaVbiexp <- t(c(x, startingvoltage, predict(fit2,newdata=data.frame(time=0))-startingvoltage));
            outputDeltaVbiexp <- t(c(x, A2, slowampl+fastampl));
            write.table(outputDeltaVbiexp, file=file.path(tpvdir,"outputDeltaVbiexp.txt"), append=TRUE, col.names=F, row.names=F, quote=F);
            
            biexpsuccess = 1;
            break
          }, error=function(e) cat("Biexp/Fit: Error ", e$message, "\n"));
        };
        if(biexpsuccess) {break}
        else {
          removedPoints=removedPoints+1
          temp <- tail(temp,-removedPoints*2)
          print(paste("Biexp/Fit: Removed", removedPoints*2, "points, fitting", nrow(temp), "datapoints"));
          step <- min((step+0.5),5)
        }
      }
    }
    
    
    
    
    
    
    
    
    
    for(rand in seq(1, 10, by=0.5)){
      tryCatch({
        print("Monoexp/Fit: Performing");
        
        if(doBiexpFit){
          Crand=sign(slowdecay*fastdecay)*sqrt(abs(slowdecay*fastdecay))*(2^runif(1,-rand,rand))
        }else{
          Crand=C*(2^runif(1,-rand,rand))
        }
        
        fit <- nls(voltage ~ cbind(1, exp(-time/C)), start=list(C=Crand),trace=F,data=temp,alg="plinear");
        capture.output(summary(fit), file=file.path(tpvdir,paste(x, "-fit", sep="")),  append=TRUE);
        A <- coef(fit)[".lin1"]; B <- coef(fit)[".lin2"]; C <- coef(fit)["C"];
        C.err <- summary(fit)$coefficients["C",2];
        
        if(doPlots){
          print("Monoexp/Plot/Linear: Performing");
	  xlim_monoexp = c(head(temp_nonfit$time,1), tail(temp$time,1))
	  ylim_monoexp = c(min(temp$voltage), max(temp$voltage))
          dev.monoexp = plotMonoexpStart(in.tpvdir=tpvdir, in.samplename=x, in.data.nonfit=temp_nonfit, in.data.fit=temp, in.data.loess=tempsubset2, in.loess_fit=lo2, in.suffix="monoexp", in.log="", in.xlim=xlim_monoexp, in.ylim=ylim_monoexp, in.color=mycolors2, in.plotHist2d=plotHist2d)
          plotMonoexpAddline(in.data.nonfit=temp_nonfit, in.data.fit=temp, in.monoexp_fit=fit, in.color=mycolors[1], in.mtext="Exp", in.dev=dev.monoexp)
        }
        if(logy & doPlots){
          print("Monoexp/Plot/Log: Performing");
          dev.monoexp.logy = plotMonoexpStart(in.tpvdir=tpvdir, in.samplename=x, in.data.nonfit=temp_nonfit, in.data.fit=temp, in.data.loess=tempsubset2, in.loess_fit=lo2, in.yshift=-A, in.suffix="monoexp-log", in.log="y", in.xlim=c(0, C*5), in.ylim=c(deltavoltage/1e5, deltavoltage), in.image_width=image_width, in.image_height=image_height, in.plotHist2d=F)
          plotMonoexpAddline(in.data.nonfit=temp_nonfit, in.data.fit=temp, in.monoexp_fit=fit, in.yshift=-A, in.dev=dev.monoexp.logy)
        }
        if(logx & doPlots){
          print("Monoexp/Plot/LogX: Performing");
          dev.monoexp.logx = plotMonoexpStart(in.tpvdir=tpvdir, in.samplename=x, in.data.nonfit=temp_nonfit, in.data.fit=temp, in.data.loess=tempsubset2, in.loess_fit=lo2, in.suffix="monoexp-logx", in.log="x", in.xlim=NULL, in.ylim=NULL, in.plotHist2d=F)
          plotMonoexpAddline(in.data.nonfit=temp_nonfit, in.data.fit=temp, in.monoexp_fit=fit, in.dev=dev.monoexp.logx)
        }
        if(residuals & doPlots){
          print("Monoexp/Plot/Residuals: Performing");
          png(file.path(tpvdir,paste(x, "-monoexp-residuals.png", sep="")), width = image_width, height = image_height);
          yresidualfit <- temp$voltage-predict(fit,newdata=data.frame(time=temp$time))
          plot(temp$time, yresidualfit, ylim=c(quantile(yresidualfit,0.001),quantile(yresidualfit,0.999)), xlab = "Time (s)", ylab = "Voltage (V)", pch=".", log="x", col=rgb(0,0,0, 0.5));
          abline(h=0, col="red");
          mtext(paste("T =", signif(C,digits=4), "s"), side=3, line=-5, adj=NA, col="red", cex=2);
          graphics.off();
        }	
        outputmonoexp <- t(c(x, A, B, C, C.err));
        write.table(outputmonoexp, file=file.path(tpvdir,"output-monoexp.txt"), append=TRUE, col.names=F, row.names=F, quote=F);
        
        #		outputDeltaVmonoexp <- t(c(x, startingvoltage, predict(fit,newdata=data.frame(time=0))-startingvoltage));
        outputDeltaVmonoexp <- t(c(x, A, B));
        write.table(outputDeltaVmonoexp, file=file.path(tpvdir,"outputDeltaVmonoexp.txt"), append=TRUE, col.names=F, row.names=F, quote=F);
        
        if(doBiexpFit){
          tryCatch({
            capture.output(anova(fit,fit2), file=file.path(tpvdir,paste(x, "-fit", sep="")),  append=TRUE);
          }, error=function(e) cat("Anova comparison: Error ", e$message, "\n"));
        }
        break;
      }, error=function(e) cat("Monoexp: Error ", e$message, "\n"));
    }
    
    if(robust){
      for(rand in seq(0, 10, by=0.5)){
        tryCatch({
          print("RobustMonoexp/Fit: Performing");
          
          CRrand=C*(2^runif(1,-rand,rand))
          fitR <- nlrob(voltage ~ A + B*exp(-time/C), start=list(A=A,B=B,C=CRrand),trace=F,data=temp);
          capture.output(summary(fitR), file=file.path(tpvdir,paste(x, "-fit", sep="")),  append=TRUE);
          AR <- coef(fitR)["A"]; BR <- coef(fitR)["B"]; 
          CR <- coef(fitR)["C"];
          CR.err <- summary(fitR)$coefficients["C",2];
          
          if(doPlots){
            print("RobustMonoexp/Plot/Linear: Performing");
            #plotMonoexpStart(in.tpvdir=tpvdir, in.samplename=x, in.data.nonfit=temp_nonfit, in.data.fit=temp, in.data.loess=tempsubset2, in.loess_fit=lo2, in.suffix="robustmonoexp", in.log="", in.xlim=NULL, in.ylim=NULL)
            plotMonoexpAddline(in.data.nonfit=temp_nonfit, in.data.fit=temp, in.monoexp_fit=fitR, in.color=mycolors[2], in.mtext="Robust Exp", in.mtextline=-7, in.dev=dev.monoexp)
          }
          if(logy & doPlots){
            print("RobustMonoexp/Plot/Log: Performing");
            #plot_monoexp(in.tpvdir=tpvdir, in.samplename=x, in.data.nonfit=temp_nonfit, in.data.fit=temp, in.data.loess=tempsubset2, in.monoexp_fit=fitR, in.loess_fit=lo2, in.yshift=-AR, in.suffix="robustmonoexp-log", in.log="y", in.xlim=c(0, CR*5), in.ylim=c(deltavoltage/1e5, deltavoltage), in.image_width=image_width, in.image_height=image_height)
            plotMonoexpAddline(in.data.nonfit=temp_nonfit, in.data.fit=temp, in.monoexp_fit=fitR, in.yshift=-AR, in.dev=dev.monoexp.logy)
            
          }
          if(logx & doPlots){
            print("RobustMonoexp/Plot/LogX: Performing");
            #plot_monoexp(in.tpvdir=tpvdir, in.samplename=x, in.data.nonfit=temp_nonfit, in.data.fit=temp, in.data.loess=tempsubset2, in.monoexp_fit=fitR, in.loess_fit=lo2, in.suffix="robustmonoexp-logx", in.log="x", in.xlim=NULL, in.ylim=NULL)
            plotMonoexpAddline(in.data.nonfit=temp_nonfit, in.data.fit=temp, in.monoexp_fit=fitR, in.dev=dev.monoexp.logx)
          }
          if(residuals & doPlots){
            print("RobustMonoexp/Plot/Residuals: Performing");
            png(file.path(tpvdir,paste(x, "-robustmonoexp-residuals.png", sep="")), width = image_width, height = image_height);
            yresidualfitR <- temp$voltage-predict(fitR,newdata=data.frame(time=temp$time))
            plot(temp$time, yresidualfitR, ylim=c(quantile(yresidualfitR,0.001),quantile(yresidualfitR,0.999)), xlab = "Time (s)", ylab = "Voltage (V)", pch=".", log="x", col=rgb(0,0,0, 0.5));
            abline(h=0, col="red");
            mtext(paste("Tau =", signif(CR,digits=4), "\u00b1", signif(CR.err,digits=4), "s"), side=3, line=-5, adj=NA, col="red", cex=2);
            graphics.off();
          }	
          outputmonoexp <- t(c(x, AR, BR, CR, CR.err));
          write.table(outputmonoexp, file=file.path(tpvdir,"output-robustmonoexp.txt"), append=TRUE, col.names=F, row.names=F, quote=F);
          
          tryCatch({
            capture.output(anova(fit,fitR), file=file.path(tpvdir,paste(x, "-fit", sep="")),  append=TRUE);
          }, error=function(e) cat("Anova comparison: Error ", e$message, "\n"));
          break;
        }, error=function(e) cat("RobustMonoexp: Error ", e$message, "\n"));
      }
    }
    
    if(robust & doBiexpFit){
      print("RobustBiexp/Fit: Performing")
      robustbiexpsuccess = 0
      tries=0
      while(!robustbiexpsuccess & tries < thresholdRobustBiexp){
        tries=tries+1
        tryCatch({
          fit2R <- nlrob(voltage ~ A2R + C2Rampl*exp(-time/C2R) + F2Rampl*exp(-time/F2R), start=list(A2R=A2, C2Rampl=slowampl, C2R=slowdecay, F2Rampl=fastampl, F2R=fastdecay),trace=F,data=temp);
          A2R <- coef(fit2R)["A2R"]; A2R.err <- summary(fit2R)$coefficients["A2R",2];
          if(coef(fit2R)["C2Rampl"] < 0 | coef(fit2R)["F2Rampl"] < 0){print("RobustBiexp/Fit: Bad result");}
          Rslowampl<-coef(fit2R)["C2Rampl"]
          Rslowampl.err<-summary(fit2R)$coefficients["C2Rampl",2]
          Rfastampl<-coef(fit2R)["F2Rampl"]
          Rfastampl.err<-summary(fit2R)$coefficients["F2Rampl",2]
          Rslowdecay<-coef(fit2R)["C2R"]
          Rslowdecay.err<-summary(fit2R)$coefficients["C2R",2]
          Rfastdecay<-coef(fit2R)["F2R"]
          Rfastdecay.err<-summary(fit2R)$coefficients["F2R",2]
          
          capture.output(summary(fit2R), file=file.path(tpvdir,paste(x, "-fit", sep="")),  append=TRUE);
          if(doPlots){
            print("RobustBiexp/Plot/Linear: Performing");
            tryCatch({
              png(file.path(tpvdir,paste(x, "-robustbiexp.png", sep="")), width = image_width, height = image_height);
              plot(mydata[[x]]$time, mydata[[x]]$voltage, xlab = "Time (s)", ylab = "Voltage (V)", pch=".", col="yellow");
              points(temp, pch=".");
              lines(tempsubset2$time, predict(lo2), col="black", lwd=3);
              lines(temp$time, predict(fit2R), col="aquamarine4", lwd=2);
              # TO BE FIXED: REMOVE HEADER VARIABLE
              #lines(temp$time+header[3]/10*header[1], A2R + Rslowampl*exp(-temp$time/Rslowdecay)+deltavoltage/10, col="green");
              mtext(paste("Tau1 =", signif(Rslowdecay,digits=4), "\u00b1", signif(Rslowdecay.err,digits=4), "s"), side=3, line=-5, adj=NA, col="green", cex=2);
              # TO BE FIXED: REMOVE HEADER VARIABLE
              #lines(temp$time+header[3]/10*header[1], A2R + Rfastampl*exp(-temp$time/Rfastdecay)+deltavoltage/10, col="blue");
              mtext(paste("Tau2 =", signif(Rfastdecay,digits=4), "\u00b1", signif(Rfastdecay.err,digits=4), "s"), side=3, line=-7, adj=NA, col="blue", cex=2);
              graphics.off();
            }, error=function(e) cat("RobustBiexp/Plot/Linear: Error ", e$message, "\n"));
          }
          if(logy & doPlots){			print("RobustBiexp/Plot/Log: Performing");
            tryCatch({
              Rquitelowerpointbiexp <- max((predict(fit2R, newdata = data.frame(time=tail(tempsubset$time, n=1)))-A2R)/5,min(subset(tempsubset, voltage > A2R, select=voltage)-A2R));
              Rhigherpointbiexp <- max(predict(fit2R, newdata = data.frame(time=0))-A2R,head(tempsubset$voltage,n=1)-A2R);
              png(file.path(tpvdir,paste(x, "-robustbiexp-log.png", sep="")), width = image_width, height = image_height);
              plot(tempsubset$time, tempsubset$voltage - A2R, xlab = "Time (s)", ylab = paste("Log(Voltage (V) -", A2R, "V)"), pch=".", log="y", ylim=c(Rquitelowerpointbiexp,Rhigherpointbiexp), col="yellow");
              points(temp$time, temp$voltage - A2R, pch=".");
              lines(tempsubset2$time, predict(lo2) - A2R, col='black', lwd=3);
              lines(temp$time,predict(fit2R) - A2R, col="aquamarine4", lwd=2);
              mtext(paste("Tau1 =", signif(Rslowdecay,digits=4), "\u00b1", signif(Rslowdecay.err,digits=4), "s"), side=3, line=-5, adj=NA, col="green", cex=2);
              mtext(paste("Tau2 =", signif(Rfastdecay,digits=4), "\u00b1", signif(Rfastdecay.err,digits=4), "s"), side=3, line=-7, adj=NA, col="blue", cex=2);
              graphics.off();	
            }, error=function(e) cat("RobustBiexp/Plot/Log: Error ", e$message, "\n"));
          }
          if(logx & doPlots){
            print("RobustBiexp/Plot/LogX: Performing");
            tryCatch({
              png(file.path(tpvdir,paste(x, "-robustbiexp-logx.png", sep="")), width = image_width, height = image_height);
              plot(temp$time, temp$voltage,  xlab = "Log(Time (s))", ylab = "Voltage (V)", pch=".", log="x", xlim=c(max(temp$time[1],-temp$time[1]),endtime));
              lines(temp$time,predict(fit2R), col="aquamarine4", lwd=2);
              mtext(paste("Tau1 =", signif(Rslowdecay,digits=4), "\u00b1", signif(Rslowdecay.err,digits=4), "s"), side=3, line=-5, adj=NA, col="green", cex=2);
              mtext(paste("Tau2 =", signif(Rfastdecay,digits=4), "\u00b1", signif(Rfastdecay.err,digits=4), "s"), side=3, line=-7, adj=NA, col="blue", cex=2);
              graphics.off();	
            }, error=function(e) cat("RobustBiexp/Plot/LogX: Error ", e$message, "\n"));
          }
          if(residuals & doPlots){
            tryCatch({
              print("RobustBiexp/Plot/Residuals: Performing");
              png(file.path(tpvdir,paste(x, "-robustbiexp-residuals.png", sep="")), width = image_width, height = image_height);
              yresidualfit2R <- temp$voltage-predict(fit2R,newdata=data.frame(time=temp$time))
              plot(temp$time, yresidualfit2R, ylim=c(quantile(yresidualfit2R,0.001),quantile(yresidualfit2R,0.999)), xlab = "Time (s)", ylab = "Voltage (V)", pch=".", log="x", col=rgb(0,0,0, (((endtime*8)-temp$time)/(endtime*8))^10));
              abline(h=0, col="red");
              mtext(paste("Tau1 =", signif(Rslowdecay,digits=4), "\u00b1", signif(Rslowdecay.err,digits=4), "s"), side=3, line=-5, adj=NA, col="green", cex=2);
              mtext(paste("Tau2 =", signif(Rfastdecay,digits=4), "\u00b1", signif(Rfastdecay.err,digits=4), "s"), side=3, line=-7, adj=NA, col="blue", cex=2);
              graphics.off();
            }, error=function(e) cat("RobustBiexp/Plot/Residuals: Error ", e$message, "\n"));
          }
          ## output for importing
          outputrobust <- t(c(x, A2R, Rslowampl, Rslowdecay, Rslowdecay.err, Rfastampl, Rfastdecay, Rfastdecay.err));
          write.table(outputrobust, file=file.path(tpvdir,"output-robustbiexp.txt"), append=TRUE, col.names=F, row.names=F, quote=F);
          
          robustbiexpsuccess = 1;
          break
        }, error=function(e) cat("RobustBiexp/Fit: Error ", e$message, "\n"));
        if(robustbiexpsuccess) {break}
        else {print("RobustBiexp/Fit: Removing a point for helping the fitting");
          temp <- temp[-1,]
        }
      }
    }
    if(doPlots){
      graphics.off(); # close all the graphical devices
      #reset the plotting margins
      #par(op)
    }
  })
  
}



plotMonoexpStart <- function(in.tpvdir="tpv", in.samplename, in.data.nonfit, in.data.fit, in.data.loess, in.loess_fit, in.yshift=0, in.suffix="monoexp", in.log="", in.xlim=NULL, in.ylim=NULL, in.output_pdf=F, in.color="black", in.plotHist2d=F){
  if(output_pdf){
    pdf(file.path(in.tpvdir,paste(in.samplename, "-", in.suffix, ".pdf", sep="")), width = image_bigpdf_width, height = image_bigpdf_height, pointsize=7);
  }else{
    png(file.path(in.tpvdir,paste(in.samplename, "-", in.suffix, ".png", sep="")), width = image_width, height = image_height);
  }
  op <- par(mar = c(5,8,4,2) + 0.1) ## default is c(5,4,4,2) + 0.1
  if(in.plotHist2d){
    steplength = function(x){sign(length(x))*(length(x)+10)}
    df = data.frame(c(in.data.nonfit$time, in.data.fit$time), c(in.data.nonfit$voltage, in.data.fit$voltage) + in.yshift)
    names(df) = c("time", "voltage")
    df = subset(df, voltage > in.ylim[1] & voltage < in.ylim[2])
    h2 = hist2d(df, col=in.color, xlab="", ylab="", cex.lab=1.7, cex.axis=1.4, xaxt="n", yaxt="n", FUN=steplength, panel.first=c(abline(h=0, col="gray80"),abline(v=0, col="gray80")))#xlim=in.xlim, ylim=in.ylim,
  }else{
    plot(in.data.fit$time, in.data.fit$voltage + in.yshift, xlab = "", ylab = "", pch=".", cex.lab=1.7, cex.axis=1.4, xaxt="n", yaxt="n", xlim=in.xlim, ylim=in.ylim, log=in.log, col="gray30", panel.first=c(abline(h=0, col="gray80"),abline(v=0, col="gray80")));
    points(in.data.nonfit$time, in.data.nonfit$voltage + in.yshift, pch=".", col="yellow");
    lines(in.data.loess$time, predict(in.loess_fit) + in.yshift, col='black', lwd=2);
  }
  eaxis(side=1, cex.axis=1.4)
  eaxis(side=2, cex.axis=1.4)
  #line is for introducing more space between label and axis
  title(ylab = "Voltage (V)", cex.lab = 1.7, line = 5)
  title(xlab = "Time (s)", cex.lab = 1.7, line = 3.5)
  return(dev.cur())
}

plotMonoexpAddline <- function(in.data.nonfit, in.data.fit, in.monoexp_fit, in.yshift=0, in.color="red", in.mtext="", in.mtextline=-5, in.dev){
  dev.set(in.dev)
  # print(paste("Plotting on", dev.cur()))
  xarray = seq(head(in.data.fit$time,1), tail(in.data.fit$time,1), length.out=100)
  ydata = predict(in.monoexp_fit, newdata = data.frame(time = xarray))
  lines(xarray, ydata + in.yshift, col=in.color, lwd=2);
  if(is.na(coef(in.monoexp_fit)["A"])){
    intercept=coef(in.monoexp_fit)[".lin1"]
  } else {
    intercept=coef(in.monoexp_fit)["A"]
  }
  segments(head(in.data.nonfit$time, 1), intercept + in.yshift, 0, intercept + in.yshift, col=in.color, lwd=2)
  mtext(paste(in.mtext, "T =", signif(coef(in.monoexp_fit)["C"],digits=4), "s"), side=3, line=in.mtextline, adj=0.95, col=in.color, cex=1.7);
}
