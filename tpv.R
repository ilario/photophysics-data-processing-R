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
library(robustbase)

powerlaw=0
robust=0
logy=0
logx=1
residuals=1
thresholdBiexp=70
thresholdRobustBiexp=100

files <- list.files(path=tpvdir, pattern="^TPV.*\\.txt.table$");
mydata <- lapply(file.path(tpvdir,files), read.table, header=FALSE, col.names=c("time","voltage"));
files <- sub(".txt.table","",files);
names(mydata) <- files;
## output for importing
write.table(t(c("file","Voc","A1","T1","T1.error","A2","T2","T2.error")), file=file.path(tpvdir,"output-biexp.txt"), append=FALSE, col.names=F, row.names=F);
write.table(t(c("file","Voc","A","T","T.error")), file=file.path(tpvdir,"output-monoexp.txt"), append=FALSE, col.names=F, row.names=F);
write.table(t(c("file","Voc","deltaV")), file=file.path(tpvdir,"outputDeltaV.txt"), append=FALSE, col.names=F, row.names=F);
write.table(t(c("file","Voc","A","T","T.error")), file=file.path(tpvdir,"output-robustmonoexp.txt"), append=FALSE, col.names=F, row.names=F);
write.table(t(c("file","Voc","A1","T1","T1.error","A2","T2","T2.error")), file=file.path(tpvdir,"output-robustbiexp.txt"), append=FALSE, col.names=F, row.names=F);




trashfornullmessages <- lapply(files, function(x) {
	message(x);
#	if(!file.exists(paste(x, "-biexp-log.png", sep="")) | !file.exists(paste(x, "-monoexp-log.png", sep=""))){
		#names(mydata[[x]]) <- c("time","voltage")
	
		#workaround
		header = as.numeric(system(paste("head -6 '", file.path(tpvdir,paste(x,".txt",sep="")), "' | tail -3|sed 's/\r$//' | cut -f2 -d' '", sep=""), intern = TRUE))


		indexpeaktime <- which.max(mydata[[x]]$voltage);
		peaktime <- mydata[[x]]$time[indexpeaktime];
		voltage2 <- subset(mydata[[x]], time >= peaktime, select=voltage); 
		time2 <- subset(mydata[[x]], time >= peaktime, select=time);
		temp <- data.frame(time2, voltage2);
		#names(temp) <- c("time","voltage")
		maxtemp=max(temp$time)
		startingvoltage <- mean(mydata[[x]]$voltage[0:10]);
		outputDeltaV <- t(c(x, startingvoltage, max(temp$voltage)-startingvoltage));
		write.table(outputDeltaV, file=file.path(tpvdir,"outputDeltaV.txt"), append=TRUE, col.names=F, row.names=F, quote=F);
		deltavoltage <- max(mydata[[x]]$voltage) - startingvoltage;
		expectedresult <- 0.00002;
		A=startingvoltage; B=deltavoltage; 
		C=expectedresult; slowdecay <- C; fastdecay <- C; #in caso non venisse effettuato il monoexp esplorativo
		tempsubset <- subset(mydata[[x]], time >= peaktime & time < tail(temp$time, n=1)/20, select=c(time,voltage));
		tempsubset2 <- subset(mydata[[x]], time > tail(temp$time, n=1)/120, select=c(time,voltage));
		lo2 <- loess(tempsubset2$voltage~tempsubset2$time, span=0.02);
		tryCatch({
		print("EvaluationMonoexp/Fit: Performing");
		fit <- nls(voltage ~ cbind(1, exp(-time/C)), start=list(C=expectedresult*runif(1,0.5,2)),trace=F,data=temp,alg="plinear");
		capture.output(summary(fit), file=file.path(tpvdir,paste(x, "-fit", sep="")),  append=TRUE);
		C <- coef(fit)["C"];
		C.err <- summary(fit)$coefficients["C",2];
		slowdecay <- C; fastdecay <- C; #in caso non venisse effettuato il biexp
		}, error=function(e) {print("EvaluationMonoexp: Error")});
		
		biexpsuccess = 0
		print("Biexp/Fit: Performing")
		step <- 0.5;
		removedPoints <- 0
		while(!biexpsuccess & removedPoints < thresholdBiexp){
		for(rand in (step*seq(1, 10, by=step))){
		while(mean(temp[1:10,2]) < mean(temp[10:20,2])) {print("Biexp/Fit: Removing a point for clearing the curve - rough"); temp <- temp[-1,]}
		while(temp[1,2] < mean(temp[2:4,2]) | temp[1,2] < mean(temp[10:15,2])) {print("Biexp/Fit: Removing a point for clearing the curve - fine"); temp <- temp[-1,]}
tryCatch({
		Crand <- C*(2^runif(1,-rand,rand))*10;
		Frand <- C*(2^runif(1,-rand,rand))/10;
		fit2 <- nls(voltage ~ cbind(1, exp(-time/C), exp(-time/F)), start=list(C=Crand, F=Frand),trace=F,data=temp,alg="plinear");
			A2 <- coef(fit2)[".lin1"]; A2.err <- summary(fit2)$coefficients[".lin1",2];
			if(coef(fit2)[".lin2"] < 0 | coef(fit2)[".lin3"] < 0 | coef(fit2)["C"] < 0 | coef(fit2)["F"] < 0){print("Biexp/Fit: Bad result"); next};
			#print(paste(startingvoltage, "=>", coef(fit2)[".lin1"], "/// ... =>", coef(fit2)[".lin2"], "///", Crand, "=>", coef(fit2)["C"], "/// ... =>", coef(fit2)[".lin3"], "--", Frand, "=>", coef(fit2)["F"]));
			if(coef(fit2)["C"] > coef(fit2)["F"]) {slowampl <- coef(fit2)[".lin2"]; slowdecay <- coef(fit2)["C"]; fastampl <- coef(fit2)[".lin3"]; fastdecay <- coef(fit2)["F"]; slowampl.err <- summary(fit2)$coefficients[".lin2",2]; slowdecay.err <- summary(fit2)$coefficients["C",2]; fastampl.err <- summary(fit2)$coefficients[".lin3",2]; fastdecay.err <- summary(fit2)$coefficients["F",2]} else {fastampl <- coef(fit2)[".lin2"]; fastdecay <- coef(fit2)["C"]; slowampl <- coef(fit2)[".lin3"]; slowdecay <- coef(fit2)["F"]; fastampl.err <- summary(fit2)$coefficients[".lin2",2]; fastdecay.err <- summary(fit2)$coefficients["C",2]; slowampl.err <- summary(fit2)$coefficients[".lin3",2]; slowdecay.err <- summary(fit2)$coefficients["F",2]};
			capture.output(summary(fit2), file=file.path(tpvdir,paste(x, "-fit", sep="")),  append=TRUE);
			
		outputDeltaVbiexp <- t(c(x, startingvoltage, predict(fit2,newdata=data.frame(time=0))-startingvoltage));
		write.table(outputDeltaVbiexp, file=file.path(tpvdir,"outputDeltaVbiexp.txt"), append=TRUE, col.names=F, row.names=F, quote=F);

			print("Biexp/Plot/Linear: Performing");
		tryCatch({
			png(file.path(tpvdir,paste(x, "-biexp.png", sep="")), width = 1280, height = 800);
			plot(mydata[[x]]$time, mydata[[x]]$voltage, main=paste(x,"biexp"), xlab = "Time (s)", ylab = "Voltage (V)", pch=3, col="yellow");
			points(temp, pch=3);
			lines(tempsubset2$time, predict(lo2), col="black", lwd=3);
			lines(temp$time, predict(fit2), col="aquamarine4", lwd=2);
			lines(temp$time+header[3]/10*header[1], A2 + slowampl*exp(-temp$time/slowdecay)+deltavoltage/10, col="green");
			mtext(paste("Tau1 =", signif(slowdecay,digits=4), "\u00b1", signif(slowdecay.err,digits=4), "s"), side=3, line=-5, adj=NA, col="green", cex=2);
			lines(temp$time+header[3]/10*header[1], A2 + fastampl*exp(-temp$time/fastdecay)+deltavoltage/10, col="blue");
			mtext(paste("Tau2 =", signif(fastdecay,digits=4), "\u00b1", signif(fastdecay.err,digits=4), "s"), side=3, line=-7, adj=NA, col="blue", cex=2);
			graphics.off();
		}, error=function(e) print("Biexp/Plot/Linear: Error"));
if(logy){			print("Biexp/Plot/Log: Performing");
		tryCatch({
			quitelowerpointbiexp <- max((predict(fit2, newdata = data.frame(time=tail(tempsubset$time, n=1)))-A2)/5,min(subset(tempsubset, voltage > A2, select=voltage)-A2));
			higherpointbiexp <- max(predict(fit2, newdata = data.frame(time=0))-A2,head(tempsubset$voltage,n=1)-A2);
			png(file.path(tpvdir,paste(x, "-biexp-log.png", sep="")), width = 1280, height = 800);
			plot(tempsubset$time, tempsubset$voltage - A2, main=paste(x,"biexp LOG"), xlab = "Time (s)", ylab = paste("Log(Voltage (V) -", A2, "V)"), pch=3, log="y", ylim=c(quitelowerpointbiexp,higherpointbiexp), col="yellow");
			points(temp$time, temp$voltage - A2, pch=3);
			lines(tempsubset2$time, predict(lo2) - A2, col='black', lwd=3);
			lines(temp$time,predict(fit2) - A2, col="aquamarine4", lwd=2);
			mtext(paste("Tau1 =", signif(slowdecay,digits=4), "\u00b1", signif(slowdecay.err,digits=4), "s"), side=3, line=-5, adj=NA, col="green", cex=2);
			mtext(paste("Tau2 =", signif(fastdecay,digits=4), "\u00b1", signif(fastdecay.err,digits=4), "s"), side=3, line=-7, adj=NA, col="blue", cex=2);
			graphics.off();	
		}, error=function(e) print("Biexp/Plot/Log: Error"));
}
if(logx){
			print("Biexp/Plot/LogX: Performing");
		tryCatch({

			png(file.path(tpvdir,paste(x, "-biexp-logx.png", sep="")), width = 1280, height = 800);
			plot(temp$time, temp$voltage, main=paste(x,"biexp LogX"), xlab = "Log(Time (s))", ylab = "Voltage (V)", pch=3, log="x", xlim=c(max(temp$time[1],-temp$time[1]),tail(temp$time,1)));
			#lines(temp$time, predict(lo2), col='black', lwd=3);
			lines(temp$time,predict(fit2), col="aquamarine4", lwd=2);
			mtext(paste("Tau1 =", signif(slowdecay,digits=4), "\u00b1", signif(slowdecay.err,digits=4), "s"), side=3, line=-5, adj=NA, col="green", cex=2);
			mtext(paste("Tau2 =", signif(fastdecay,digits=4), "\u00b1", signif(fastdecay.err,digits=4), "s"), side=3, line=-7, adj=NA, col="blue", cex=2);
			graphics.off();	
		}, error=function(e) print("Biexp/Plot/LogX: Error"));
}
if(residuals){
		tryCatch({
			print("Biexp/Plot/Residuals: Performing");
			png(file.path(tpvdir,paste(x, "-biexp-residuals.png", sep="")), width = 1280, height = 800);
			yresidualfit2 <- temp$voltage-predict(fit2,newdata=data.frame(time=temp$time))
			plot(temp$time, yresidualfit2, ylim=c(quantile(yresidualfit2,0.001),quantile(yresidualfit2,0.999)), main=paste(x,"residuals"), xlab = "Time (s)", ylab = "Voltage (V)", pch=3, log="x", col=rgb(0,0,0, (((maxtemp*8)-temp$time)/(maxtemp*8))^10));
			abline(h=0, col="red");
			mtext(paste("Tau1 =", signif(slowdecay,digits=4), "\u00b1", signif(slowdecay.err,digits=4), "s"), side=3, line=-5, adj=NA, col="green", cex=2);
			mtext(paste("Tau2 =", signif(fastdecay,digits=4), "\u00b1", signif(fastdecay.err,digits=4), "s"), side=3, line=-7, adj=NA, col="blue", cex=2);
			graphics.off();
		}, error=function(e) print("Biexp/Plot/Residuals: Error"));
}
## output for importing
		output <- t(c(x, A2, slowampl, slowdecay, slowdecay.err, fastampl, fastdecay, fastdecay.err));
		write.table(output, file=file.path(tpvdir,"output-biexp.txt"), append=TRUE, col.names=F, row.names=F, quote=F);
		
		biexpsuccess = 1;
		break
}, error=function(e) print("Biexp/Fit: Error"));
		};
		if(biexpsuccess) {break}
		else {
		removedPoints=removedPoints+1
		print(paste("Biexp/Fit: Removing", removedPoints, "points for helping the fitting"));
		temp <- temp[-removedPoints,]
		step <- min((step+0.5),5)
		}
		}
		
#		if(biexpsuccess) {
		for(rand in seq(1, 10, by=0.5)){
tryCatch({
		print("Monoexp/Fit: Performing");

		Crand=sqrt(slowdecay*fastdecay)*(2^runif(1,-rand,rand))
		fit <- nls(voltage ~ cbind(1, exp(-time/C)), start=list(C=Crand),trace=F,data=temp,alg="plinear");
		capture.output(summary(fit), file=file.path(tpvdir,paste(x, "-fit", sep="")),  append=TRUE);
		A <- coef(fit)[".lin1"]; B <- coef(fit)[".lin2"]; C <- coef(fit)["C"];
		C.err <- summary(fit)$coefficients["C",2];
		quitelowerpointmonoexp <- max((predict(fit, newdata = data.frame(time=tail(tempsubset$time, n=1)))-A)/5,min(subset(tempsubset, voltage > A, select=voltage)-A));
		higherpointmonoexp <- max(predict(fit, newdata = data.frame(time=0))-A,head(tempsubset$voltage,n=1)-A);

		outputDeltaVmonoexp <- t(c(x, startingvoltage, predict(fit,newdata=data.frame(time=0))-startingvoltage));
		write.table(outputDeltaVmonoexp, file=file.path(tpvdir,"outputDeltaVmonoexp.txt"), append=TRUE, col.names=F, row.names=F, quote=F);
		
		print("Monoexp/Plot/Linear: Performing");
		png(file.path(tpvdir,paste(x, "-monoexp.png", sep="")), width = 1280, height = 800);
		plot(mydata[[x]]$time, mydata[[x]]$voltage, main=paste(x,"monoexp"), xlab = "Time (s)", ylab = "Voltage (V)", pch=3, col="yellow");
		points(temp, pch=3);
		lines(tempsubset2$time, predict(lo2), col='black', lwd=3);
		lines(temp$time, predict(fit), col="red", lwd=2);
		mtext(paste("Tau =", signif(C,digits=4), "\u00b1", signif(C.err,digits=4), "s"), side=3, line=-5, adj=NA, col="red", cex=2);
		graphics.off();
if(logy){
		print("Monoexp/Plot/Log: Performing");
		png(file.path(tpvdir,paste(x, "-monoexp-log.png", sep="")), width = 1280, height = 800);
		plot(tempsubset$time, tempsubset$voltage - A, main=paste(x,"monoexp LOG"), xlab = "Time (s)", ylab = paste("Log(Voltage (V) -", A, "V)"), pch=3, log="y", ylim=c(quitelowerpointmonoexp,higherpointmonoexp), col="yellow");
		points(temp$time, temp$voltage - A, pch=3);
		lines(tempsubset2$time, predict(lo2) - A, col='black', lwd=3);
		lines(temp$time,predict(fit) - A, col="red", lwd=2);
		mtext(paste("Tau =", signif(C,digits=4), "\u00b1", signif(C.err,digits=4), "s"), side=3, line=-5, adj=NA, col="red", cex=2);
		graphics.off();
}
if(logx){
		print("Monoexp/Plot/LogX: Performing");
		png(file.path(tpvdir,paste(x, "-monoexp-logx.png", sep="")), width = 1280, height = 800);
		plot(temp$time, temp$voltage, main=paste(x,"monoexp LogX"), xlab = "Log(Time (s))", ylab = "Voltage (V)", pch=3, log="x");
		#lines(temp$time, predict(lo2), col='black', lwd=3);
		lines(temp$time,predict(fit), col="red", lwd=2);
		mtext(paste("Tau =", signif(C,digits=4), "\u00b1", signif(C.err,digits=4), "s"), side=3, line=-5, adj=NA, col="red", cex=2);
		graphics.off();
}
if(residuals){
		print("Monoexp/Plot/Residuals: Performing");
		png(file.path(tpvdir,paste(x, "-monoexp-residuals.png", sep="")), width = 1280, height = 800);
		yresidualfit <- temp$voltage-predict(fit,newdata=data.frame(time=temp$time))
                plot(temp$time, yresidualfit, ylim=c(quantile(yresidualfit,0.001),quantile(yresidualfit,0.999)), main=paste(x,"residuals"), xlab = "Time (s)", ylab = "Voltage (V)", pch=3, log="x", col=rgb(0,0,0, 0.5));
		abline(h=0, col="red");
		mtext(paste("Tau =", signif(C,digits=4), "\u00b1", signif(C.err,digits=4), "s"), side=3, line=-5, adj=NA, col="red", cex=2);
		graphics.off();
}	
                outputmonoexp <- t(c(x, A, B, C, C.err));
                write.table(outputmonoexp, file=file.path(tpvdir,"output-monoexp.txt"), append=TRUE, col.names=F, row.names=F, quote=F);

tryCatch({
			capture.output(anova(fit,fit2), file=file.path(tpvdir,paste(x, "-fit", sep="")),  append=TRUE);
}, error=function(e) print("Anova comparison: Error"));
		break;
}, error=function(e) print("Monoexp: Error"));
		}

if(powerlaw){
		for(rand in seq(1, 10, by=0.5)){
		tryCatch({
		print("PowerMonoexp/Fit: Performing");

		Grand=sqrt(slowdecay*fastdecay)*(2^runif(1,-rand,rand))
		Hrand=1000*(2^runif(1,-rand,rand))
		fit3 <- nls(voltage ~ cbind(1, exp(-time/G), 1/(1+H*time)), start=list(G=Grand, H=Hrand),trace=F,data=temp,alg="plinear");
		A3 <- coef(fit3)[".lin1"]; A3.err <- summary(fit3)$coefficients[".lin1",2];
		Gampl <- coef(fit3)[".lin2"]; Gampl.err <- summary(fit3)$coefficients[".lin2",2];
	       	G <- coef(fit3)["G"];
		G.err <- summary(fit3)$coefficients["G",2];
		Hampl <- coef(fit3)[".lin3"]; Hampl.err <- summary(fit3)$coefficients[".lin3",2];
	       	H <- coef(fit3)["H"];
		H.err <- summary(fit3)$coefficients["H",2];
		capture.output(summary(fit3), file=file.path(tpvdir,paste(x, "-fit", sep="")),  append=TRUE);
		K <- H/G;
		K.err <- K*(H.err/H + G.err/G)

		print("PowerMonoexp/Plot/Linear: Performing");
		tryCatch({
			png(file.path(tpvdir,paste(x, "-powermonoexp.png", sep="")), width = 1280, height = 800);
			plot(mydata[[x]]$time, mydata[[x]]$voltage, main=paste(x,"power + monoexp"), xlab = "Time (s)", ylab = "Voltage (V)", pch=3, col="yellow");
			points(temp, pch=3);
			lines(tempsubset2$time, predict(lo2), col="black", lwd=3);
			lines(temp$time, predict(fit3), col="red", lwd=2);
			lines(temp$time+header[3]/10*header[1], A3 + Gampl*exp(-temp$time/G)+deltavoltage/10, col="green");
			mtext(paste("Tau =", signif(G,digits=4), "\u00b1", signif(G.err,digits=2), "s"), side=3, line=-5, adj=NA, col="green", cex=2);
			lines(temp$time+header[3]/10*header[1], A3 + Hampl/(1+H*temp$time)+deltavoltage/10, col="blue");
			mtext(paste("K =", signif(K,digits=4), "\u00b1", signif(K.err,digits=2), "s"), side=3, line=-7, adj=NA, col="blue", cex=2);
			graphics.off();
		}, error=function(e) print("PowerMonoexp/Plot/Linear: Error"));
if(logy){
			print("PowerMonoexp/Plot/Log: Performing");
		tryCatch({
			png(file.path(tpvdir,paste(x, "-powermonoexp-log.png", sep="")), width = 1280, height = 800);
			plot(tempsubset$time, tempsubset$voltage - A3, main=paste(x,"power + monoexp LOG"), xlab = "Time (s)", ylab = paste("Log(Voltage (V) -", A3, "V)"), pch=3, log="y", ylim=c(quitelowerpointbiexp,higherpointbiexp), col="yellow");
			points(temp$time, temp$voltage - A3, pch=3);
			lines(tempsubset2$time, predict(lo2) - A3, col='black', lwd=3);
			lines(temp$time,predict(fit3) - A3, col="red", lwd=2);
			mtext(paste("Tau =", signif(G,digits=4), "\u00b1", signif(G.err,digits=2), "s"), side=3, line=-5, adj=NA, col="green", cex=2);
			mtext(paste("K =", signif(K,digits=4), "\u00b1", signif(K.err,digits=2), "s"), side=3, line=-7, adj=NA, col="blue", cex=2);
			graphics.off();	
		}, error=function(e) print("PowerMonoexp/Plot/Log: Error"));
}
if(residuals){
		tryCatch({
			print("PowerMonoexp/Plot/Residuals: Performing");
			png(file.path(tpvdir,paste(x, "-powermonoexp-residuals.png", sep="")), width = 1280, height = 800);
			yresidualfit3 <- temp$voltage-predict(fit3,newdata=data.frame(time=temp$time))
                        plot(temp$time, yresidualfit3, ylim=c(quantile(yresidualfit3,0.001),quantile(yresidualfit3,0.999)), main=paste(x,"power + monoexp residuals"), xlab = "Time (s)", ylab = "Voltage (V)", pch=3, log="x", col=rgb(0,0,0, (((maxtemp*8)-temp$time)/(maxtemp*8))^10));
			abline(h=0, col="red");
			mtext(paste("Tau =", signif(G,digits=4), "\u00b1", signif(G.err,digits=2), "s"), side=3, line=-5, adj=NA, col="green", cex=2);
			mtext(paste("K =", signif(K,digits=4), "\u00b1", signif(K.err,digits=2), "s"), side=3, line=-7, adj=NA, col="blue", cex=2);
			graphics.off();
		}, error=function(e) print("PowerMonoexp/Plot/Residuals: Error"));
}
		break
		}, error=function(e) print("PowerMonoexp: Error"))}
}


if(robust){
		for(rand in seq(0, 10, by=0.5)){
		#tryCatch({
		print("RobustMonoexp/Fit: Performing");

		CRrand=C*(2^runif(1,-rand,rand))
		#fitR <- nlrob(voltage ~ cbind(1, exp(-time/CR)), start=list(CR=CRrand),trace=F,data=temp,alg="plinear");
		fitR <- nlrob(voltage ~ AR + BR*exp(-time/CR), start=list(AR=A,BR=B,CR=CRrand),trace=F,data=temp);
		capture.output(summary(fitR), file=file.path(tpvdir,paste(x, "-fit", sep="")),  append=TRUE);
#		AR <- coef(fitR)[".lin1"]; BR <- coef(fitR)[".lin2"]; 
		AR <- coef(fitR)["AR"]; BR <- coef(fitR)["BR"]; 
		CR <- coef(fitR)["CR"];
		CR.err <- summary(fitR)$coefficients["CR",2];
		quitelowerpointmonoexp <- max((predict(fitR, newdata = data.frame(time=tail(tempsubset$time, n=1)))-AR)/5,min(subset(tempsubset, voltage > AR, select=voltage)-AR));
		higherpointmonoexp <- max(predict(fitR, newdata = data.frame(time=0))-AR,head(tempsubset$voltage,n=1)-AR);
		print("RobustMonoexp/Plot/Linear: Performing");
		png(file.path(tpvdir,paste(x, "-robustmonoexp.png", sep="")), width = 1280, height = 800);
		plot(mydata[[x]]$time, mydata[[x]]$voltage, main=paste(x,"robust monoexp"), xlab = "Time (s)", ylab = "Voltage (V)", pch=3, col="yellow");
		points(temp, pch=3);
		lines(tempsubset2$time, predict(lo2), col='black', lwd=3);
		lines(temp$time, predict(fitR), col="red", lwd=2);
		mtext(paste("Tau =", signif(CR,digits=4), "\u00b1", signif(CR.err,digits=4), "s"), side=3, line=-5, adj=NA, col="red", cex=2);
		graphics.off();
if(logy){
		print("RobustMonoexp/Plot/Log: Performing");
		png(file.path(tpvdir,paste(x, "-robustmonoexp-log.png", sep="")), width = 1280, height = 800);
		plot(tempsubset$time, tempsubset$voltage - AR, main=paste(x,"robust monoexp LOG"), xlab = "Time (s)", ylab = paste("Log(Voltage (V) -", AR, "V)"), pch=3, log="y", ylim=c(quitelowerpointmonoexp,higherpointmonoexp), col="yellow");
		points(temp$time, temp$voltage - AR, pch=3);
		lines(tempsubset2$time, predict(lo2) - AR, col='black', lwd=3);
		lines(temp$time,predict(fitR) - AR, col="red", lwd=2);
		mtext(paste("Tau =", signif(CR,digits=4), "\u00b1", signif(CR.err,digits=4), "s"), side=3, line=-5, adj=NA, col="red", cex=2);
		graphics.off();
}
if(logx){
		print("RobustMonoexp/Plot/LogX: Performing");
		png(file.path(tpvdir,paste(x, "-robustmonoexp-logx.png", sep="")), width = 1280, height = 800);
		plot(temp$time, temp$voltage, main=paste(x,"robust monoexp LogX"), xlab = "Log(Time (s))", ylab = "Voltage (V)", pch=3, log="x");
		#lines(temp$time, predict(lo2), col='black', lwd=3);
		lines(temp$time,predict(fitR), col="red", lwd=2);
		mtext(paste("Tau =", signif(CR,digits=4), "\u00b1", signif(CR.err,digits=4), "s"), side=3, line=-5, adj=NA, col="red", cex=2);
		graphics.off();
}
if(residuals){
		print("RobustMonoexp/Plot/Residuals: Performing");
		png(file.path(tpvdir,paste(x, "-robustmonoexp-residuals.png", sep="")), width = 1280, height = 800);
		yresidualfitR <- temp$voltage-predict(fitR,newdata=data.frame(time=temp$time))
                plot(temp$time, yresidualfitR, ylim=c(quantile(yresidualfitR,0.001),quantile(yresidualfitR,0.999)), main=paste(x,"robust residuals"), xlab = "Time (s)", ylab = "Voltage (V)", pch=3, log="x", col=rgb(0,0,0, 0.5));
		abline(h=0, col="red");
		mtext(paste("Tau =", signif(CR,digits=4), "\u00b1", signif(CR.err,digits=4), "s"), side=3, line=-5, adj=NA, col="red", cex=2);
		graphics.off();
}	
                outputmonoexp <- t(c(x, AR, BR, CR, CR.err));
                write.table(outputmonoexp, file=file.path(tpvdir,"output-robustmonoexp.txt"), append=TRUE, col.names=F, row.names=F, quote=F);

		tryCatch({
			capture.output(anova(fit,fitR), file=file.path(tpvdir,paste(x, "-fit", sep="")),  append=TRUE);
		}, error=function(e) print("Anova comparison: Error"));
		break;
#	}, error=function(e) print("RobustMonoexp: Error"));
		}
}




if(powerlaw){if(robust){
		for(rand in seq(0, 10, by=0.5)){
		tryCatch({
		print("RobustPowerMonoexp/Fit: Performing");

		Irand=G*(2^runif(1,-rand,rand))
		Lrand=H*(2^runif(1,-rand,rand))
		fit3r <- nlrob(voltage ~ A3r + Iampl*exp(-time/I) + Lampl/(1+L*time), start=list(A3r=A3, Iampl=Gampl, I=Irand, Lampl=Hampl, L=Lrand),trace=F,data=temp); #na.action?
		#A3r <- coef(fit3r)[".lin1"]; A3r.err <- summary(fit3r)$coefficients[".lin1",2];
		#Iampl <- coef(fit3r)[".lin2"]; Iampl.err <- summary(fit3r)$coefficients[".lin2",2];
		A3r <- coef(fit3r)["A3r"]; A3r.err <- summary(fit3r)$coefficients["A3r",2];
		Iampl <- coef(fit3r)["Iampl"]; Iampl.err <- summary(fit3r)$coefficients["Iampl",2];
	       	I <- coef(fit3r)["I"];
		I.err <- summary(fit3r)$coefficients["I",2];
		#Lampl <- coef(fit3r)[".lin3"]; Lampl.err <- summary(fit3r)$coefficients[".lin3",2];
		Lampl <- coef(fit3r)["Lampl"]; Lampl.err <- summary(fit3r)$coefficients["Lampl",2];
	       	L <- coef(fit3r)["L"];
		L.err <- summary(fit3r)$coefficients["L",2];
		capture.output(summary(fit3r), file=file.path(tpvdir,paste(x, "-fit", sep="")),  append=TRUE);
		Kr <- I/L;
		Kr.err <- Kr*(I.err/I + L.err/L)

		print("RobustPowerMonoexp/Plot/Linear: Performing");
		tryCatch({
			png(file.path(tpvdir,paste(x, "-robustpowermonoexp.png", sep="")), width = 1280, height = 800);
			plot(mydata[[x]]$time, mydata[[x]]$voltage, main=paste(x,"robust power + monoexp"), xlab = "Time (s)", ylab = "Voltage (V)", pch=3, col="yellow");
			points(temp, pch=3);
			lines(tempsubset2$time, predict(lo2), col="black", lwd=3);
			lines(temp$time, predict(fit3r), col="red", lwd=2);
			lines(temp$time+header[3]/10*header[1], A3r + Iampl*exp(-temp$time/I)+deltavoltage/10, col="green");
			mtext(paste("Tau =", signif(I,digits=4), "\u00b1", signif(I.err,digits=2), "s"), side=3, line=-5, adj=NA, col="green", cex=2);
			lines(temp$time+header[3]/10*header[1], A3r + Lampl/(1+L*temp$time)+deltavoltage/10, col="blue");
			mtext(paste("K =", signif(Kr,digits=4), "\u00b1", signif(Kr.err,digits=2), "s"), side=3, line=-7, adj=NA, col="blue", cex=2);
			graphics.off();
		}, error=function(e) print("RobustPowerMonoexp/Plot/Linear: Error"));

if(logy){
			print("RobustPowerMonoexp/Plot/Log: Performing");
		tryCatch({
			png(file.path(tpvdir,paste(x, "-robustpowermonoexp-log.png", sep="")), width = 1280, height = 800);
			plot(tempsubset$time, tempsubset$voltage - A3r, main=paste(x,"robust power + monoexp LOG"), xlab = "Time (s)", ylab = paste("Log(Voltage (V) -", A3r, "V)"), pch=3, log="y", ylim=c(quitelowerpointbiexp,higherpointbiexp), col="yellow");
			points(temp$time, temp$voltage - A3r, pch=3);
			lines(tempsubset2$time, predict(lo2) - A3r, col='black', lwd=3);
			lines(temp$time,predict(fit3r) - A3r, col="red", lwd=2);
			mtext(paste("Tau =", signif(I,digits=4), "\u00b1", signif(I.err,digits=2), "s"), side=3, line=-5, adj=NA, col="green", cex=2);
			mtext(paste("K =", signif(Kr,digits=4), "\u00b1", signif(Kr.err,digits=2), "s"), side=3, line=-7, adj=NA, col="blue", cex=2);
			graphics.off();	
		}, error=function(e) print("RobustPowerMonoexp/Plot/Log: Error"));
}
if(residuals){
		tryCatch({
			print("RobustPowerMonoexp/Plot/Residuals: Performing");
			png(file.path(tpvdir,paste(x, "-robustpowermonoexp-residuals.png", sep="")), width = 1280, height = 800);
			yresidualfit3r <- temp$voltage-predict(fit3r,newdata=data.frame(time=temp$time))
                        plot(temp$time, yresidualfit3r, ylim=c(quantile(yresidualfit3r,0.001),quantile(yresidualfit3r,0.999)), main=paste(x,"robust power + monoexp residuals"), xlab = "Time (s)", ylab = "Voltage (V)", pch=3, log="x", col=rgb(0,0,0, (((maxtemp*8)-temp$time)/(maxtemp*8))^10));
			abline(h=0, col="red");
			mtext(paste("Tau =", signif(I,digits=4), "\u00b1", signif(I.err,digits=2), "s"), side=3, line=-5, adj=NA, col="green", cex=2);
			mtext(paste("K =", signif(Kr,digits=4), "\u00b1", signif(Kr.err,digits=2), "s"), side=3, line=-7, adj=NA, col="blue", cex=2);
			graphics.off();
		}, error=function(e) print("RobustPowerMonoexp/Plot/Residuals: Error"));
}
		break
		}, error=function(e) print("RobustPowerMonoexp: Error"))
		}
		}}



if(robust){
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

			print("RobustBiexp/Plot/Linear: Performing");
		tryCatch({
			png(file.path(tpvdir,paste(x, "-robustbiexp.png", sep="")), width = 1280, height = 800);
			plot(mydata[[x]]$time, mydata[[x]]$voltage, main=paste(x,"robust biexp"), xlab = "Time (s)", ylab = "Voltage (V)", pch=3, col="yellow");
			points(temp, pch=3);
			lines(tempsubset2$time, predict(lo2), col="black", lwd=3);
			lines(temp$time, predict(fit2R), col="aquamarine4", lwd=2);
			lines(temp$time+header[3]/10*header[1], A2R + Rslowampl*exp(-temp$time/Rslowdecay)+deltavoltage/10, col="green");
			mtext(paste("Tau1 =", signif(Rslowdecay,digits=4), "\u00b1", signif(Rslowdecay.err,digits=4), "s"), side=3, line=-5, adj=NA, col="green", cex=2);
			lines(temp$time+header[3]/10*header[1], A2R + Rfastampl*exp(-temp$time/Rfastdecay)+deltavoltage/10, col="blue");
			mtext(paste("Tau2 =", signif(Rfastdecay,digits=4), "\u00b1", signif(Rfastdecay.err,digits=4), "s"), side=3, line=-7, adj=NA, col="blue", cex=2);
			graphics.off();
		}, error=function(e) print("RobustBiexp/Plot/Linear: Error"));
if(logy){			print("RobustBiexp/Plot/Log: Performing");
		tryCatch({
			Rquitelowerpointbiexp <- max((predict(fit2R, newdata = data.frame(time=tail(tempsubset$time, n=1)))-A2R)/5,min(subset(tempsubset, voltage > A2R, select=voltage)-A2R));
			Rhigherpointbiexp <- max(predict(fit2R, newdata = data.frame(time=0))-A2R,head(tempsubset$voltage,n=1)-A2R);
			png(file.path(tpvdir,paste(x, "-robustbiexp-log.png", sep="")), width = 1280, height = 800);
			plot(tempsubset$time, tempsubset$voltage - A2R, main=paste(x,"robust biexp LOG"), xlab = "Time (s)", ylab = paste("Log(Voltage (V) -", A2R, "V)"), pch=3, log="y", ylim=c(Rquitelowerpointbiexp,Rhigherpointbiexp), col="yellow");
			points(temp$time, temp$voltage - A2R, pch=3);
			lines(tempsubset2$time, predict(lo2) - A2R, col='black', lwd=3);
			lines(temp$time,predict(fit2R) - A2R, col="aquamarine4", lwd=2);
			mtext(paste("Tau1 =", signif(Rslowdecay,digits=4), "\u00b1", signif(Rslowdecay.err,digits=4), "s"), side=3, line=-5, adj=NA, col="green", cex=2);
			mtext(paste("Tau2 =", signif(Rfastdecay,digits=4), "\u00b1", signif(Rfastdecay.err,digits=4), "s"), side=3, line=-7, adj=NA, col="blue", cex=2);
			graphics.off();	
		}, error=function(e) print("RobustBiexp/Plot/Log: Error"));
}
if(logx){
			print("RobustBiexp/Plot/LogX: Performing");
		tryCatch({
			png(file.path(tpvdir,paste(x, "-robustbiexp-logx.png", sep="")), width = 1280, height = 800);
			plot(temp$time, temp$voltage, main=paste(x,"robust biexp LogX"), xlab = "Log(Time (s))", ylab = "Voltage (V)", pch=3, log="x", xlim=c(max(temp$time[1],-temp$time[1]),tail(temp$time,1)));
			lines(temp$time,predict(fit2R), col="aquamarine4", lwd=2);
			mtext(paste("Tau1 =", signif(Rslowdecay,digits=4), "\u00b1", signif(Rslowdecay.err,digits=4), "s"), side=3, line=-5, adj=NA, col="green", cex=2);
			mtext(paste("Tau2 =", signif(Rfastdecay,digits=4), "\u00b1", signif(Rfastdecay.err,digits=4), "s"), side=3, line=-7, adj=NA, col="blue", cex=2);
			graphics.off();	
		}, error=function(e) print("RobustBiexp/Plot/LogX: Error"));
}
if(residuals){
		tryCatch({
			print("RobustBiexp/Plot/Residuals: Performing");
			png(file.path(tpvdir,paste(x, "-robustbiexp-residuals.png", sep="")), width = 1280, height = 800);
			yresidualfit2R <- temp$voltage-predict(fit2R,newdata=data.frame(time=temp$time))
                        plot(temp$time, yresidualfit2R, ylim=c(quantile(yresidualfit2R,0.001),quantile(yresidualfit2R,0.999)), main=paste(x,"robust residuals"), xlab = "Time (s)", ylab = "Voltage (V)", pch=3, log="x", col=rgb(0,0,0, (((maxtemp*8)-temp$time)/(maxtemp*8))^10));
			abline(h=0, col="red");
			mtext(paste("Tau1 =", signif(Rslowdecay,digits=4), "\u00b1", signif(Rslowdecay.err,digits=4), "s"), side=3, line=-5, adj=NA, col="green", cex=2);
			mtext(paste("Tau2 =", signif(Rfastdecay,digits=4), "\u00b1", signif(Rfastdecay.err,digits=4), "s"), side=3, line=-7, adj=NA, col="blue", cex=2);
			graphics.off();
		}, error=function(e) print("RobustBiexp/Plot/Residuals: Error"));
}
## output for importing
		outputrobust <- t(c(x, A2R, Rslowampl, Rslowdecay, Rslowdecay.err, Rfastampl, Rfastdecay, Rfastdecay.err));
		write.table(outputrobust, file=file.path(tpvdir,"output-robustbiexp.txt"), append=TRUE, col.names=F, row.names=F, quote=F);
		
		robustbiexpsuccess = 1;
		break
	}, error=function(e) print("RobustBiexp/Fit: Error"));
		if(robustbiexpsuccess) {break}
		else {print("RobustBiexp/Fit: Removing a point for helping the fitting");
		temp <- temp[-1,]
		}
		}
#}




#		}
	}})
}
