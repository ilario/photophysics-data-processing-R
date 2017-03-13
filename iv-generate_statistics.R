name="ig"

variable1="PCBM"
v1g1=c()
v1g2=c()
v1g3=c()
v1g4=c()
v1all=c(v1g1, v1g2, v1g3, v1g4)

#variable2="Vacuum"
#v2g1=c()
#v2g2=c()
#v2g3=c()
#v2all=c(v2g1, v2g2, v2g3)
#v2g1=c()
#v2g2=all[! all %in% v2g1]

a <- read.table("table", col.names = c("device","diode","reverse","Jsc","Voc","FF","PCE","notes"),fill = TRUE, na.strings=c("","NA"))
levels(a$notes) <- c(levels(a$notes), "NANA")
a[is.na(a)] <- "NANA"

colors<-c("blue","blue","darkviolet","darkviolet","darkred","darkred")

b <- data.frame(a, temp1 = factor(a$device, v1all))
colnames(b)[ncol(b)] <- variable1
levels(b[[variable1]]) <- list(`` = v1g1, `` = v1g2, `` = v1g3, `` = v1g4);
#z <- data.frame(b, temp2 = factor(a$device, v2all))#(b$notes, v2all, exclude=NULL))
#colnames(z)[ncol(z)] <- variable2
#levels(z[[variable2]]) <- list(`` = v2g1,`` = v2g2, ``=v2g3);
z <- b

d <- data.frame(z, rev = factor(z$reverse, c(0,1)))
levels(d$rev) <- list('fwd' = 0, 'rev' = 1)


png(paste(name, "-", variable1, ".png", sep=""));
#pdf(paste(name, "-", variable1, ".pdf", sep=""), width=10, height=10);
par(mar=c(8,6,1,1), mgp = c(4, 1, 0))
boxplot(PCE~droplevels(interaction(rev,PCBM)), d#[d$HTM=="",]
	, cex.axis=1.5, cex.lab=1.5, border=colors, las=2, ylab="PCE (%)")#,main=paste(name, variable1))#,xaxt="n")
#axis(side=1,labels=c("",""), at=seq(1,), las=2, cex.axis=2)
graphics.off()

png(paste(name, "-", variable1, "-Jsc.png", sep=""));
#pdf(paste(name, "-", variable1, ".pdf", sep=""), width=10, height=10);
par(mar=c(8,6,1,1), mgp = c(4, 1, 0))
boxplot(Jsc~droplevels(interaction(rev,PCBM)), d, cex.axis=1.5, cex.lab=1.5, border=colors, las=2, ylab=bquote("J"["sc"]~"(mA/cm"^"2"*")"))#,main=paste(name, variable1))#,xaxt="n")
#axis(side=1,labels=c("",""), at=seq(1,), las=2, cex.axis=2)
graphics.off()

png(paste(name, "-", variable1, "-Voc.png", sep=""));
#pdf(paste(name, "-", variable1, ".pdf", sep=""), width=10, height=10);
par(mar=c(8,6,1,1), mgp = c(4, 1, 0))
boxplot(Voc~droplevels(interaction(rev,PCBM)), d, cex.axis=1.5, cex.lab=1.5, border=colors, las=2, ylab=bquote("V"["oc"]~"(V)"))#,main=paste(name, variable1))#,xaxt="n")
#axis(side=1,labels=c("",""), at=seq(1,), las=2, cex.axis=2)
graphics.off()

png(paste(name, "-", variable1, "-FF.png", sep=""));
#pdf(paste(name, "-", variable1, ".pdf", sep=""), width=10, height=10);
par(mar=c(8,6,1,1), mgp = c(4, 1, 0))
boxplot(FF~droplevels(interaction(rev,PCBM)), d, cex.axis=1.5, cex.lab=1.5, border=colors, las=2, ylab="FF")#,main=paste(name, variable1))#,xaxt="n")
#axis(side=1,labels=c("",""), at=seq(1,), las=2, cex.axis=2)
graphics.off()



if(exists("variable2")){
png(paste(name, "-", variable2, ".png", sep=""));
#pdf(paste(name, "-", variable2, ".pdf", sep=""), width=10, height=10);
par(mar=c(8,6,1,1), mgp = c(4, 1, 0))
boxplot(PCE~droplevels(interaction(rev,Vacuum)), d#[d$HTM=="",]
	, cex.axis=1.5, cex.lab=1.5, border=colors, las=2, ylab="PCE (%)")#,main=paste(name, variable2))#,xaxt="n")
#axis(side=1,labels=c("",""), at=seq(1,), las=2, cex.axis=2)
graphics.off()

png(paste(name, "-", variable2, "-Jsc.png", sep=""));
#pdf(paste(name, "-", variable2, ".pdf", sep=""), width=10, height=10);
par(mar=c(8,6,1,1), mgp = c(4, 1, 0))
boxplot(Jsc~droplevels(interaction(rev,Vacuum)), d, cex.axis=1.5, cex.lab=1.5, border=colors, las=2, ylab=bquote("J"["sc"]~"(mA/cm"^"2"*")"))#,main=paste(name, variable2))#,xaxt="n")
#axis(side=1,labels=c("",""), at=seq(1,), las=2, cex.axis=2)
graphics.off()

png(paste(name, "-", variable2, "-Voc.png", sep=""));
#pdf(paste(name, "-", variable2, ".pdf", sep=""), width=10, height=10);
par(mar=c(8,6,1,1), mgp = c(4, 1, 0))
boxplot(Voc~droplevels(interaction(rev,Vacuum)), d, cex.axis=1.5, cex.lab=1.5, border=colors, las=2, ylab=bquote("V"["oc"]~"(V)"))#,main=paste(name, variable2))#,xaxt="n")
#axis(side=1,labels=c("",""), at=seq(1,), las=2, cex.axis=2)
graphics.off()

png(paste(name, "-", variable2, "-FF.png", sep=""));
#pdf(paste(name, "-", variable2, ".pdf", sep=""), width=10, height=10);
par(mar=c(8,6,1,1), mgp = c(4, 1, 0))
boxplot(FF~droplevels(interaction(rev,Vacuum)), d, cex.axis=1.5, cex.lab=1.5, border=colors, las=2, ylab="FF")#,main=paste(name, variable2))#,xaxt="n")
#axis(side=1,labels=c("",""), at=seq(1,), las=2, cex.axis=2)
graphics.off()





png(paste(name, "-", variable1, "-", variable2, ".png", sep=""));
#pdf(paste(name, "-", variable2, ".pdf", sep=""), width=10, height=10);
par(mar=c(8,6,1,1), mgp = c(4, 1, 0))
boxplot(PCE~droplevels(interaction(rev,PCBM,Vacuum)), d#[d$HTM=="",]
	, cex.axis=1.5, cex.lab=1.5, border=colors, las=2, ylab="PCE (%)")#,main=paste(name, variable1, variable2))#,xaxt="n")
#axis(side=1,labels=c("",""), at=seq(1,), las=2, cex.axis=2)
graphics.off()

png(paste(name,  "-", variable1,"-", variable2, "-Jsc.png", sep=""));
#pdf(paste(name, "-", variable2, ".pdf", sep=""), width=10, height=10);
par(mar=c(8,6,1,1), mgp = c(4, 1, 0))
boxplot(Jsc~droplevels(interaction(rev,PCBM,Vacuum)), d, cex.axis=1.5, cex.lab=1.5, border=colors, las=2, ylab=bquote("J"["sc"]~"(mA/cm"^"2"*")"))#,main=paste(name, variable1, variable2))#,xaxt="n")
#axis(side=1,labels=c("",""), at=seq(1,), las=2, cex.axis=2)
graphics.off()

png(paste(name,  "-", variable1,"-", variable2, "-Voc.png", sep=""));
#pdf(paste(name, "-", variable2, ".pdf", sep=""), width=10, height=10);
par(mar=c(8,6,1,1), mgp = c(4, 1, 0))
boxplot(Voc~droplevels(interaction(rev,PCBM,Vacuum)), d, cex.axis=1.5, cex.lab=1.5, border=colors, las=2, ylab=bquote("V"["oc"]~"(V)"))#,main=paste(name, variable1, variable2))#,xaxt="n")
#axis(side=1,labels=c("",""), at=seq(1,), las=2, cex.axis=2)
graphics.off()

png(paste(name,  "-", variable1,"-", variable2, "-FF.png", sep=""));
#pdf(paste(name, "-", variable2, ".pdf", sep=""), width=10, height=10);
par(mar=c(8,6,1,1), mgp = c(4, 1, 0))
boxplot(FF~droplevels(interaction(rev,PCBM,Vacuum)), d, cex.axis=1.5, cex.lab=1.5, border=colors, las=2, ylab="FF")#,main=paste(name, variable1, variable2))#,xaxt="n")
#axis(side=1,labels=c("",""), at=seq(1,), las=2, cex.axis=2)
graphics.off()
}
