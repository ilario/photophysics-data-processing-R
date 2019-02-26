ce <- function(cedir="ce")
{
  library(sfsmisc)
  library(RColorBrewer)
  
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
  
  print("CE: INTEGRATING")
  files <- list.files(path=cedir, pattern="^CE.*\\.txt.table$");
  mydata <- lapply(file.path(cedir,files), read.table, header=FALSE, col.names=c("time","voltage"));
  files <- sub(".txt.table","",files);
  names(mydata) <- files;
  write.table(t(c("Voc","ChargeDensityCE")), file=file.path(cedir,"outputChargeDensityCE.txt"), append=FALSE, col.names=F, row.names=F);
  write.table(t(c("Voc","ChargeCE")), file=file.path(cedir,"outputChargeCE.txt"), append=FALSE, col.names=F, row.names=F);
  
  print(files[1])
  darkCEvoltage = mydata[[files[1]]]$voltage;
  
  trashfornullmessages <- lapply(files, function(x) {
    message(x);
    len<-length(mydata[[x]]$voltage)
    startVoltage <- mean(mydata[[x]]$voltage[1:600])
    endVoltage <- mean(mydata[[x]]$voltage[(len-600):len])
    baseline <- seq(startVoltage, endVoltage, length.out=len)
    voltage2 <- mydata[[x]]$voltage - baseline
    
    current <- voltage2/50
    charge <- cumsum(current)*(mydata[[x]]$time[2]-mydata[[x]]$time[1])
    
    charge=charge-charge[match(0,mydata[[x]]$time)]
    totalcharge=mean(charge[round(length(charge)*0.9):round(length(charge)*0.95)])
    totalchargedensity=totalcharge/cellArea
    
    #voltagezero <- mydata[[x]]$voltage
    #currentzero <- voltagezero/50
    #chargezero <- cumsum(currentzero)*(mydata[[x]]$time[2]-mydata[[x]]$time[1])
    #chargezero=chargezero-chargezero[match(0,mydata[[x]]$time)]
    
    b<-strsplit(x, "_")
    c<-unlist(b)
    c2 <- c[grepl("mV",c)]
    d<-as.numeric(sub("mV.*", "", c2))
    outputChargeDensityCE <- t(c(d, totalchargedensity));
    outputChargeCE <- t(c(d, totalcharge));
    write.table(outputChargeDensityCE, file=file.path(cedir,"outputChargeDensityCE.txt"), append=TRUE, col.names=F, row.names=F, quote=F);
    
    if(output_pdf){
      pdf(file.path(cedir,paste(x, ".pdf", sep="")), width=image_bigpdf_width, height=image_bigpdf_height, pointsize=7)
    }else{
      png(file.path(cedir,paste(x, ".png", sep="")), width=image_width, height=image_height)
    }
    op <- par(mar = c(5,8.5,2,8.5) + 0.1) ## default is c(5,4,4,2) + 0.1
    
    xlim=c(1e-7,tail(mydata[[x]]$time,1))
    #plot(mydata[[x]],type="l", ylab="", xlab="", xaxt="n", yaxt="n", xlim=xlim, log="x", cex.axis=1.4, panel.first=c(lines(mydata[[x]]$time, baseline, col="gray70")))
    timeDownsampled = logdownsampling(mydata[[x]]$time)
    plot(timeDownsampled,logdownsampling(mydata[[x]]$voltage),type="l", ylab="", xlab="", xaxt="n", yaxt="n", xlim=xlim, log="x", cex.axis=1.4, panel.first=c(lines(mydata[[x]]$time, baseline, col="gray70")))
    
    title(ylab="Voltage (V)", cex.lab=1.7, line=6)
    title(xlab="Time (s)", cex.lab=1.7, line=3.5)
    mtext(bquote("Collected Charge Density (C/cm"^"2"*")"), cex=1.7, side=4,line=7,col=mycolors[3])#"red")
    eaxis(side=1, cex.axis=1.4)
    yaxt <- round(axTicks(2), digits=10) # without round, the zero can get printed as 2E-17
    eaxis(side=2, cex.axis=1.4, labels = pretty10exp(yaxt))
    
    #lines(mydata[[x]]$time, voltagezero - darkCEvoltage, col=mycolors[1])#"blue")
    
    #ylim_charge=c(min(chargezero, charge)/cellArea, max(chargezero, charge)/cellArea)
    par(new=TRUE)
    plot(timeDownsampled,logdownsampling(charge)/cellArea, type="l", col=mycolors[3], xaxt="n",yaxt="n",xlab="",ylab="", xlim=xlim, log="x")#ylim=ylim_charge,
    #abline(h=0,col=mycolors[3])
    eaxis(4,col.ticks=mycolors[3],col.axis=mycolors[3], col=mycolors[3], cex.axis=1.4)
    #text(xlim[2]*0.9,totalchargedensity*0.9,labels=bquote(.(signif(totalchargedensity,3))~"C/cm"^"2"),cex=1.7, adj=1)
    mtext(bquote(.(signif(totalchargedensity,3))~"C/cm"^"2"), side=3, line=-10, cex=1.7, adj=0.95)
    
    
    #par(new=TRUE)
    #plot(mydata[[x]]$time,chargezero/cellArea, type="l", col=mycolors[4], xaxt="n",yaxt="n",xlab="",ylab="", xlim=xlim, ylim=ylim_charge, log="x")
    
    #legendtext = c("Signal","Integrated charge")
    #legend(x="bottomright",inset=0,legendtext,col=c("black", mycolors[3]), cex=1.5, lwd=4, bty="n")
    
    graphics.off()
    #reset the plotting margins
    par(op)
    
    
    write.table(outputChargeDensityCE, file=file.path(cedir,"outputChargeDensityCE.txt"), append=TRUE, col.names=F, row.names=F, quote=F);
    write.table(outputChargeCE, file=file.path(cedir,"outputChargeCE.txt"), append=TRUE, col.names=F, row.names=F, quote=F);
  })
}
