name="ig"

variable1="HTM"
v1g1=c()
v1g2=c()
v1g3=c()
all=c(v1g1, v1g2, v1g3)
variable2=""
v2g1=""
v2g2=""
all2=c(v2g1, v2g2)
#v2g1=c()
#v2g2=all[! all %in% v2g1]

a <- read.table("table", col.names = c("device","diode","reverse","Jsc","Voc","FF","PCE","notes"),fill = TRUE, na.strings=c("","NA"))
levels(a$notes) <- c(levels(a$notes), "NANA")
a[is.na(a)] <- "NANA"

colors<-c("blue","blue","darkviolet","darkviolet","darkred","darkred")

b <- data.frame(a, temp1 = factor(a$device, all))
colnames(b)[ncol(b)] <- variable1
levels(b[[variable1]]) <- list(`` = v1g1, `` = v1g2, `` = v1g3);
z <- data.frame(b, temp2 = factor(b$notes, all2, exclude=NULL))
colnames(z)[ncol(z)] <- variable2
levels(z[[variable2]]) <- list(`` = v2g1,`` = v2g2);
#z <- b

d <- data.frame(z, rev = factor(z$reverse, c(0,1)))
levels(d$rev) <- list('fwd' = 0, 'rev' = 1)


png(paste(name, "-", variable1, ".png", sep=""));
#pdf(paste(name, "-", variable1, ".pdf", sep=""), width=10, height=10);
par(mar=c(13,5,1,1))
boxplot(PCE~droplevels(interaction(rev,)), d[d$HTM=="",], cex.axis=1.5, cex.lab=1.5, border=colors, las=2, ylab="PCE",main=paste(name, variable1))#,xaxt="n")
#axis(side=1,labels=c("",""), at=seq(1,), las=2, cex.axis=2)
graphics.off()


png(paste(name, "-", variable1, "-Jsc.png", sep=""));
#pdf(paste(name, "-", variable1, ".pdf", sep=""), width=10, height=10);
par(mar=c(13,5,1,1))
boxplot(Jsc~droplevels(interaction(rev,)), d[d$HTM=="",], cex.axis=1.5, cex.lab=1.5, border=colors, las=2, ylab="PCE",main=paste(name, variable1))#,xaxt="n")
#axis(side=1,labels=c("",""), at=seq(1,), las=2, cex.axis=2)
graphics.off()

png(paste(name, "-", variable1, "-Voc.png", sep=""));
#pdf(paste(name, "-", variable1, ".pdf", sep=""), width=10, height=10);
par(mar=c(13,5,1,1))
boxplot(Voc~droplevels(interaction(rev,)), d[d$HTM=="",], cex.axis=1.5, cex.lab=1.5, border=colors, las=2, ylab="PCE",main=paste(name, variable1))#,xaxt="n")
#axis(side=1,labels=c("",""), at=seq(1,), las=2, cex.axis=2)
graphics.off()

png(paste(name, "-", variable1, "-FF.png", sep=""));
#pdf(paste(name, "-", variable1, ".pdf", sep=""), width=10, height=10);
par(mar=c(13,5,1,1))
boxplot(FF~droplevels(interaction(rev,)), d[d$HTM=="",], cex.axis=1.5, cex.lab=1.5, border=colors, las=2, ylab="PCE",main=paste(name, variable1))#,xaxt="n")
#axis(side=1,labels=c("",""), at=seq(1,), las=2, cex.axis=2)
graphics.off()



if(exists("variable2")){
}
