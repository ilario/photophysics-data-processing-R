options(error=function() { traceback(2); if(!interactive()) quit("no", status = 1, runLast = FALSE) })

library(robustbase)

files <- list.files(path=".", pattern="\\.dat$");
lapply(files, function(x) {print(x);
	original <- read.csv(x, skip=6, header=F);

	#downsampling by 10 times
	downsampled <- head(cbind(round(matrix(original$V1, nrow=10)[1,]),colSums(matrix(original$V2, nrow=10))/10),-1);
	a = data.frame(position=downsampled[,1], height=downsampled[,2])

	length <- round((tail(a[,1], n=1)-a[1,1])/1e4)*10
 
	substrateToZero = data.frame(position=a$position, height=a$height-quantile(a$height,0.02))
	quantiles = quantile(substrateToZero$height, c(0.5,0.95))
        mid <- quantiles[1]
        high <- quantiles[2]
	if(output_pdf){
		pdf(paste(x,".pdf",sep=""), width=image_bigpdf_width, height=image_bigpdf_height, pointsize=7)
	}else{
		png(paste(x,".png",sep=""), width=image_width, height=image_height)
	}
		par(mar = c(4,5,1,1) + 0.1) ## default is c(5,4,4,2) + 0.1
	        plot(substrateToZero, type="l", ylab="Height (nm)", xlab="", xaxt="n", bty="n", ylim=c(0,min(1500,high*1.2)), lwd=2, cex.lab=1.7, cex.axis=1.4, panel.first=c(abline(h=0, col="gray80"), if(mid < 0.9*high){abline(h=mid, col="gray80")}, abline(h=high, col="gray80")))
        	mtext(side = 1, text = paste(length,"nm"), line=1, cex=1.7)
        	text(length*1000*0.9, high, paste(round(high),"nm"), col="red", cex=1.7, font=2)
        	if(mid < 0.9*high){text(length*1000*0.9, mid, paste(round(mid),"nm"), col="red", cex=1.7, font=2)}
		abline(h=0, col="gray50")
	dev.off()

	heightAverage = mean(a$height)
	aSubstrate = subset(a, height < heightAverage/2)
	fitParabolaSub = lmrob(height ~ poly(position, 3), data=aSubstrate)
	flattenedSub = data.frame(position=a$position, height=a$height - predict(fitParabolaSub, newdata=data.frame(position=a$position)))
	substrateToZeroSub = data.frame(position=flattenedSub$position, height=flattenedSub$height-quantile(flattenedSub$height,0.03))
	quantilesSub = quantile(substrateToZeroSub$height, c(0.5,0.95))
        midSub <- quantilesSub[1]
        highSub <- quantilesSub[2]
	if(output_pdf){
		pdf(paste(x,"-flat.pdf",sep=""), width=image_bigpdf_width, height=image_bigpdf_height, pointsize=7)
	}else{
		png(paste(x,"-flat.png",sep=""), width=image_width, height=image_height)
	}
		par(mar = c(4,5,1,1) + 0.1)
	        plot(substrateToZeroSub, type="l", ylab="Height (nm)", xlab="", xaxt="n", bty="n", ylim=c(0,min(1500,highSub*1.2)), lwd=2, cex.lab=1.7, cex.axis=1.4, panel.first=c(abline(h=0, col="gray80"), if(midSub < 0.9*highSub){abline(h=midSub, col="gray80")}, abline(h=highSub, col="gray80")))
        	mtext(side = 1, text = paste(length,"nm"), line=1, cex=1.7)
        	text(length*1000*0.9, highSub, paste(round(highSub),"nm"), col="red", cex=1.7, font=2)
        	if(midSub < 0.9*highSub){text(length*1000*0.9, midSub, paste(round(midSub),"nm"), col="red", cex=1.7, font=2)}
		abline(h=0, col="gray50")
	dev.off()

	aSurface = subset(a, height > heightAverage/2)
	fitParabolaSur = lmrob(height ~ poly(position, 3), data=aSurface)
	flattenedSur = data.frame(position=a$position, height=a$height - predict(fitParabolaSur, newdata=data.frame(position=a$position)))
	substrateToZeroSur = data.frame(position=flattenedSur$position, height=flattenedSur$height-quantile(flattenedSur$height,0.03))
	quantilesSur = quantile(substrateToZeroSur$height, c(0.5,0.95))
        midSur <- quantilesSur[1]
        highSur <- quantilesSur[2]
	if(output_pdf){
		pdf(paste(x,"-flat2.pdf",sep=""), width=image_bigpdf_width, height=image_bigpdf_height, pointsize=7)
	}else{
		png(paste(x,"-flat2.png",sep=""), width=image_width, height=image_height)
	}
		par(mar = c(4,5,1,1) + 0.1)
	        plot(substrateToZeroSur, type="l", ylab="Height (nm)", xlab="", xaxt="n", bty="n", ylim=c(0,min(1500,highSur*1.2)), lwd=2, cex.lab=1.7, cex.axis=1.4, panel.first=c(abline(h=0, col="gray80"), if(midSur < 0.9*highSur){abline(h=midSur, col="gray80")}, abline(h=highSur, col="gray80")))
        	mtext(side = 1, text = paste(length,"nm"), line=1, cex=1.7)
        	text(length*1000*0.9, highSur, paste(round(highSur),"nm"), col="red", cex=2, font=2)
		if(midSur < 0.9*highSur){text(length*1000*0.9, midSur, paste(round(midSur),"nm"), col="red", cex=2, font=2)}
		abline(h=0, col="gray50")
	dev.off()
})
