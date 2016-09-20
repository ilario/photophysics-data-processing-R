name="ig63"

variable1="TiO2"
v1g1=c()
v1g2=c()
v1all=c(v1g1, v1g2)

variable2="Vacuum"
v2g1=c()
v2g2=c()
v2g3=c()
v2all=c(v2g1, v2g2, v2g3)
#v2g1=c()
#v2g2=all[! all %in% v2g1]

a <- read.table("table", col.names = c("device","diode","reverse","Jsc","Voc","FF","PCE","notes"),fill = TRUE, na.strings=c("","NA"))
levels(a$notes) <- c(levels(a$notes), "NANA")
a[is.na(a)] <- "NANA"

colors<-c("blue","blue","darkviolet","darkviolet","darkred","darkred")

b <- data.frame(a, temp1 = factor(a$device, v1all))
colnames(b)[ncol(b)] <- variable1
levels(b[[variable1]]) <- list(`single` = v1g1, `double` = v1g2);
z <- data.frame(b, temp2 = factor(a$device, v2all))#(b$notes, v2all, exclude=NULL))
colnames(z)[ncol(z)] <- variable2
levels(z[[variable2]]) <- list(`PhCl` = v2g1,`good` = v2g2, `bad`=v2g3);
#z <- b

d <- data.frame(z, rev = factor(z$reverse, c(0,1)))
levels(d$rev) <- list('fwd' = 0, 'rev' = 1)


png(paste(name, "-", variable1, ".png", sep=""));
#pdf(paste(name, "-", variable1, ".pdf", sep=""), width=10, height=10);
par(mar=c(13,5,1,1))
boxplot(PCE~droplevels(interaction(rev,TiO2)), d#[d$HTM=="",]
	, cex.axis=1.5, cex.lab=1.5, border=colors, las=2, ylab="PCE",main=paste(name, variable1))#,xaxt="n")
#axis(side=1,labels=c("",""), at=seq(1,), las=2, cex.axis=2)
graphics.off()

png(paste(name, "-", variable1, "-Jsc.png", sep=""));
#pdf(paste(name, "-", variable1, ".pdf", sep=""), width=10, height=10);
par(mar=c(13,5,1,1))
boxplot(Jsc~droplevels(interaction(rev,TiO2)), d, cex.axis=1.5, cex.lab=1.5, border=colors, las=2, ylab="Jsc",main=paste(name, variable1))#,xaxt="n")
#axis(side=1,labels=c("",""), at=seq(1,), las=2, cex.axis=2)
graphics.off()

png(paste(name, "-", variable1, "-Voc.png", sep=""));
#pdf(paste(name, "-", variable1, ".pdf", sep=""), width=10, height=10);
par(mar=c(13,5,1,1))
boxplot(Voc~droplevels(interaction(rev,TiO2)), d, cex.axis=1.5, cex.lab=1.5, border=colors, las=2, ylab="Voc",main=paste(name, variable1))#,xaxt="n")
#axis(side=1,labels=c("",""), at=seq(1,), las=2, cex.axis=2)
graphics.off()

png(paste(name, "-", variable1, "-FF.png", sep=""));
#pdf(paste(name, "-", variable1, ".pdf", sep=""), width=10, height=10);
par(mar=c(13,5,1,1))
boxplot(FF~droplevels(interaction(rev,TiO2)), d, cex.axis=1.5, cex.lab=1.5, border=colors, las=2, ylab="FF",main=paste(name, variable1))#,xaxt="n")
#axis(side=1,labels=c("",""), at=seq(1,), las=2, cex.axis=2)
graphics.off()



if(exists("variable2")){
png(paste(name, "-", variable2, ".png", sep=""));
#pdf(paste(name, "-", variable2, ".pdf", sep=""), width=10, height=10);
par(mar=c(13,5,1,1))
boxplot(PCE~droplevels(interaction(rev,Vacuum)), d#[d$HTM=="",]
	, cex.axis=1.5, cex.lab=1.5, border=colors, las=2, ylab="PCE",main=paste(name, variable2))#,xaxt="n")
#axis(side=1,labels=c("",""), at=seq(1,), las=2, cex.axis=2)
graphics.off()

png(paste(name, "-", variable2, "-Jsc.png", sep=""));
#pdf(paste(name, "-", variable2, ".pdf", sep=""), width=10, height=10);
par(mar=c(13,5,1,1))
boxplot(Jsc~droplevels(interaction(rev,Vacuum)), d, cex.axis=1.5, cex.lab=1.5, border=colors, las=2, ylab="Jsc",main=paste(name, variable2))#,xaxt="n")
#axis(side=1,labels=c("",""), at=seq(1,), las=2, cex.axis=2)
graphics.off()

png(paste(name, "-", variable2, "-Voc.png", sep=""));
#pdf(paste(name, "-", variable2, ".pdf", sep=""), width=10, height=10);
par(mar=c(13,5,1,1))
boxplot(Voc~droplevels(interaction(rev,Vacuum)), d, cex.axis=1.5, cex.lab=1.5, border=colors, las=2, ylab="Voc",main=paste(name, variable2))#,xaxt="n")
#axis(side=1,labels=c("",""), at=seq(1,), las=2, cex.axis=2)
graphics.off()

png(paste(name, "-", variable2, "-FF.png", sep=""));
#pdf(paste(name, "-", variable2, ".pdf", sep=""), width=10, height=10);
par(mar=c(13,5,1,1))
boxplot(FF~droplevels(interaction(rev,Vacuum)), d, cex.axis=1.5, cex.lab=1.5, border=colors, las=2, ylab="FF",main=paste(name, variable2))#,xaxt="n")
#axis(side=1,labels=c("",""), at=seq(1,), las=2, cex.axis=2)
graphics.off()





png(paste(name, "-", variable1, "-", variable2, ".png", sep=""));
#pdf(paste(name, "-", variable2, ".pdf", sep=""), width=10, height=10);
par(mar=c(13,5,1,1))
boxplot(PCE~droplevels(interaction(rev,TiO2,Vacuum)), d#[d$HTM=="",]
	, cex.axis=1.5, cex.lab=1.5, border=colors, las=2, ylab="PCE",main=paste(name, variable1, variable2))#,xaxt="n")
#axis(side=1,labels=c("",""), at=seq(1,), las=2, cex.axis=2)
graphics.off()

png(paste(name,  "-", variable1,"-", variable2, "-Jsc.png", sep=""));
#pdf(paste(name, "-", variable2, ".pdf", sep=""), width=10, height=10);
par(mar=c(13,5,1,1))
boxplot(Jsc~droplevels(interaction(rev,TiO2,Vacuum)), d, cex.axis=1.5, cex.lab=1.5, border=colors, las=2, ylab="Jsc",main=paste(name, variable1, variable2))#,xaxt="n")
#axis(side=1,labels=c("",""), at=seq(1,), las=2, cex.axis=2)
graphics.off()

png(paste(name,  "-", variable1,"-", variable2, "-Voc.png", sep=""));
#pdf(paste(name, "-", variable2, ".pdf", sep=""), width=10, height=10);
par(mar=c(13,5,1,1))
boxplot(Voc~droplevels(interaction(rev,TiO2,Vacuum)), d, cex.axis=1.5, cex.lab=1.5, border=colors, las=2, ylab="Voc",main=paste(name, variable1, variable2))#,xaxt="n")
#axis(side=1,labels=c("",""), at=seq(1,), las=2, cex.axis=2)
graphics.off()

png(paste(name,  "-", variable1,"-", variable2, "-FF.png", sep=""));
#pdf(paste(name, "-", variable2, ".pdf", sep=""), width=10, height=10);
par(mar=c(13,5,1,1))
boxplot(FF~droplevels(interaction(rev,TiO2,Vacuum)), d, cex.axis=1.5, cex.lab=1.5, border=colors, las=2, ylab="FF",main=paste(name, variable1, variable2))#,xaxt="n")
#axis(side=1,labels=c("",""), at=seq(1,), las=2, cex.axis=2)
graphics.off()
}
