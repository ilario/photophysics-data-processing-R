ce.files <- list.files(path="ce", pattern="^CE.*\\.txt.table$");
tpv.files  <- list.files(path="tpv", pattern="^TPV.*\\.txt.table$");
tpc.sun.files <- list.files(path="tpc", pattern="^TPC.*sun.*\\.txt.table$");
tpc.dark.files <- list.files(path="tpc", pattern="^TPC.*dark.*\\.txt.table$");

ce.file <- tail(ce.files, n=1)
tpv.file <- tail(tpv.files, n=1)
tpc.sun.file <- tail(tpc.sun.files, n=1)
tpc.dark.file <- tail(tpc.dark.files, n=1)

ce <- read.table(paste("ce/",ce.file,sep=""),header=F)
tpc.sun <- read.table(paste("tpc/",tpc.sun.file,sep=""),header=F)
tpc.dark <- read.table(paste("tpc/",tpc.dark.file,sep=""),header=F)
tpv <- read.table(paste("tpv/",tpv.file,sep=""),header=F)
tpc.sun$V2 <- -tpc.sun$V2
tpc.dark$V2 <- -tpc.dark$V2

ce.max <- max(ce$V2)
tpc.sun.max <- max(tpc.sun$V2)
tpc.dark.max <- max(tpc.dark$V2)
tpc.sun.min <- mean(tpc.sun$V2[1:50])
tpc.dark.min <- mean(tpc.dark$V2[1:50])
tpv.max <- max(tpv$V2)
ce.min <- mean(ce$V2[1:50])
tpv.min <- mean(tpv$V2[1:50])
main <- gsub(".txt.table","",gsub("CE_","",ce.file))
main <- strsplit(main, "_")[[1]][1]
png(paste("tpc_vs_tpv_vs_ce-", main, ".png", sep=""), width=600, heigh=600)
plot(tpv$V1, (tpv$V2-tpv.min)/(tpv.max-tpv.min), xlim=c(-5e-7,0.5e-5), ylim=c(-0.1,1), type="l", ylab="Normalized Voltage", xlab="Time (s)", main=paste(main, "TPC dark vs TPC 1sun vs TPV vs CE"), col="blue", yaxt='n', cex.axis=1.4, cex.lab=1.4)
lines(ce$V1, (ce$V2-ce.min)/(ce.max-ce.min), lwd=2, col="green")
lines(tpc.sun$V1, (tpc.sun$V2-tpc.sun.min)/(max(tpc.sun.max-tpc.sun.min, tpc.dark.max-tpc.dark.min)), lwd=1, col="red")
lines(tpc.dark$V1, (tpc.dark$V2-tpc.dark.min)/(max(tpc.sun.max-tpc.sun.min, tpc.dark.max-tpc.dark.min)), lwd=1, col="black", lty=1)
legend(x="topright",inset=0.05, c("TPC 1sun","TPC dark", "TPV 1sun", "CE 1sun"), lwd=4, col=c("red","black","blue","green"), cex=2)
graphics.off()


