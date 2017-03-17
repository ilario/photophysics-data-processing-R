files <- list.files(path=".", pattern="^ig.*\\.txt$");
mydata <- lapply(files, read.table, header=T)
names(mydata) <- files
h=new.env(hash=TRUE)
h[["ig79-c25-2-redLED80mA.txt"]]=1
h[["ig79-c25-2-blueLED80mA.txt"]]=2

colors=c("red","blue")
#lty=c(3,3,3,2,2,2,1,1,1)
i=0
png("Voc-stability.png", height=600, width=600)
plot(NULL,xlim=c(0,80),ylim=c(0,1.1), xlab="Illumination Time (s)", ylab="Voc (V)", cex.axis=1.5, cex.lab=1.5)
lapply(ls(h), function(x){
	       	       lines(mydata[[x]]$TIME, mydata[[x]]$VOLT, lwd=2, col=colors[h[[x]]])
		       	       i<<-i+1
		       legend(x="bottomright",inset=0.05,c("red LED", "blue LED"), col=colors, lwd=4, cex=1.4)
})
graphics.off()

