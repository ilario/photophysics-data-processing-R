height=1200
width=1600

allfiles = list.files()
pngfiles = list.files(pattern='png')
list=allfiles[!allfiles %in% pngfiles] 

lapply(list, function(x){print(x); 
  data=read.csv(x,header=F)
  names(data)=c("mode","time","voltage","current")
  data$current=data$current*1000/0.09
  data$power=data$voltage*data$current
  endscan=tail(data[data$mode==1,],1)
  data2 = data[1:as.numeric(rownames(endscan)),]
  mpp=tail(data2[data2$mode==0,],1)
  
  names(endscan)
  print(mpp)
  
  png(paste(x,".png",sep=""), height=height, width=width)
  par(mar = c(5, 4, 5, 4) + 0.3)
  plot(NULL, xlim=range(data$time), ylim=range(data$voltage), xlab="Time (s)", ylab="Voltage (V)", cex.axis=2, cex.lab=2)
  points(data$time, data$voltage, pch=".")
  points(mpp$time, mpp$voltage, cex=3, pch=16)
  par(new=T)
  range=range(data$current, data$power)
  plot(data$time, data$current, pch=".", axes = FALSE, bty = "n", xlab = "", ylab = "", col="red", ylim=range)
  par(new=T)
  plot(data$time, data$power, pch=".", axes = FALSE, bty = "n", xlab = "", ylab = "", col="green", ylim=range)
  axis(side=4, at = pretty(range), cex.axis=2)
  mtext("Current density (mA/cm2, red), PCE (%, green)", side=4, line=3, cex=2)
  points(mpp$time,mpp$current, cex=3, pch=16, col="red")
  points(mpp$time,mpp$power, cex=3, pch=16, col="green")
  graphics.off()}
)

#png(paste(file,"-justdwell.png",sep=""), height=height, width=width)
#par(mar = c(5, 4, 5, 4) + 0.3)
#data=data[data$mode==0,]
#plot(data$time, data$voltage, pch=".", xlab="Voltage (V)", ylab="Time (s)")
#par(new=T)
#range=range(data$current,data$power)
#plot(data$time, data$current, pch=".", axes = FALSE, bty = "n", xlab = "", ylab = "", col="red", ylim=range)
#par(new=T)
#plot(data$time, data$power, pch=".", axes = FALSE, bty = "n", xlab = "", ylab = "", col="green", ylim=range)
#axis(side=4, at = pretty(range))
#mtext("Current density (mA/cm2, red), PCE (%, green)", side=4, line=3)
#graphics.off()
