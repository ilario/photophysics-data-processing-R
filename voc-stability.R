files <- list.files(path=".", pattern="^ig.*\\.txt$");
mydata <- lapply(files, read.table, header=T)
names(mydata) <- files
h=new.env(hash=TRUE)
h[["ig92-1356-1week-1.txt"]]=1
h[["ig92-1370-1week-3.txt"]]=2
h[["ig92-1395-1week-4.txt"]]=1
h[["ig92-1396-1week-2.txt"]]=2
h[["ig92-1427-1week-1.txt"]]=2
h[["ig92-1454-1week-2.txt"]]=2
h[["ig92-1456-1week-2.txt"]]=1
h[["ig92-1467-1week-3.txt"]]=1
h[["ig92-1370-1week-1.txt"]]=2
h[["ig92-1370-1weekbis-3.txt"]]=2
colors=c("red","blue")
i=0
png("Voc-stability.png", height=600, width=600)
par(mar=c(5,6,2,2))
plot(NULL,xlim=c(0,80),ylim=c(0,1.1), xlab="Illumination Time (s)", ylab=bquote("V"["OC"]*" (V)"), cex.axis=1.5, cex.lab=1.5)
lapply(ls(h), function(x){
	       deltas = mydata[[x]]$VOLT[seq(2,length(mydata[[x]]$VOLT))] - mydata[[x]]$VOLT[seq(1,length(mydata[[x]]$VOLT)-1)]
	       jump = which.max(abs(deltas))
	       	       lines(mydata[[x]]$TIME-mydata[[x]]$TIME[jump], mydata[[x]]$VOLT, lwd=2, col=colors[h[[x]]])
		       	       i<<-i+1
		       legend(x="bottomright",inset=0.05,c("red LED", "blue LED"), col=colors, lwd=4, cex=1.4)
})
graphics.off()

