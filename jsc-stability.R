files <- list.files(path=".", pattern="^ig.*\\.txt$");
mydata <- lapply(files, read.table, header=T)
names(mydata) <- files

h2=new.env(hash=TRUE)
h2[["ig101-1514-1-jsc.txt"]]=1
h2[["ig101-1553-3-jsc.txt"]]=1
h2[["ig101-1555-1-jsc.txt"]]=2
h2[["ig101-1567-2-jsc.txt"]]=2
colors=c("red","blue")
i=0
png("Jsc-stability.png", height=600, width=600)
par(mar=c(5,6,2,2))
plot(NULL,xlim=c(0,80),ylim=c(-0.003,0), xlab="Illumination Time (s)", ylab=bquote("J"["SC"]*" (mA)"), cex.axis=1.5, cex.lab=1.5)
lapply(ls(h2), function(x){
  lines(mydata[[x]]$TIME, mydata[[x]]$CURR, lwd=2, col=colors[h2[[x]]])
  i<<-i+1
  legend(x="bottomright",inset=0.05,c("spiro-OMeTAD", "TAE4"), col=colors, lwd=4, cex=1.4)
})
graphics.off()
