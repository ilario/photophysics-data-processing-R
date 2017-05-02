name="ig"

var1=""
v1g1=c()
v1g2=c()
v1all=c(v1g1, v1g2)
list1=list(`` = v1g1, `` = v1g2)

var2=""
v2g1=c()
v2g2=c()#v1all[! v1all %in% v2g1]
v2all=c(v2g1, v2g2)
list2=list(`` = v2g1, `` = v2g2)

#var3=""
#v3g1=c()
#v3g2=c()#v1all[! v1all %in% v3g1]
#v3all=c(v3g1, v3g2)
#list3=list(`` = v3g1, `` = v3g2)


a <- read.table("table", col.names = c("device","diode","reverse","Jsc","Voc","FF","PCE","notes"),fill = TRUE, na.strings=c("","NA"))
levels(a$notes) <- c(levels(a$notes), "NANA")
a[is.na(a)] <- "NANA"

colors<-c("blue","blue","darkviolet","darkviolet","darkred","darkred")

a1 <- data.frame(a, temp1 = factor(a$device, v1all))
colnames(a1)[ncol(a1)] <- var1
levels(a1[[var1]]) <- list1;
z <- a1
if(exists("var2")){
	a2 <- data.frame(a1, temp2 = factor(a$device, v2all))#(b$notes, v2all, exclude=NULL))
	colnames(a2)[ncol(a2)] <- var2
	levels(a2[[var2]]) <- list2;
	z <- a2
	if(exists("var3")){
		a3 <- data.frame(a2, temp3 = factor(a$device, v3all))#(b$notes, v2all, exclude=NULL))
		colnames(a3)[ncol(a3)] <- var3
		levels(a3[[var3]]) <- list3;
		z <- a3
		}
	}

d <- data.frame(z, rev = factor(z$reverse, c(0,1)))
levels(d$rev) <- list('fwd' = 0, 'rev' = 1)

single <- function(x,mar,mgp){
	png(paste(name, "-", x, ".png", sep=""));
	par(mar=mar, mgp = mgp)
	boxplot(PCE~droplevels(interaction(rev,get(x))), d, cex.axis=1.5, cex.lab=1.5, border=colors, las=2, ylab="PCE (%)")
	graphics.off()

	png(paste(name, "-", x, "-Jsc.png", sep=""));
	par(mar=mar, mgp = mgp)
	boxplot(Jsc~droplevels(interaction(rev,get(x))), d, cex.axis=1.5, cex.lab=1.5, border=colors, las=2, ylab=bquote("J"["sc"]~"(mA/cm"^"2"*")"))
	graphics.off()

	png(paste(name, "-", x, "-Voc.png", sep=""));
	par(mar=mar, mgp = mgp)
	boxplot(Voc~droplevels(interaction(rev,get(x))), d, cex.axis=1.5, cex.lab=1.5, border=colors, las=2, ylab=bquote("V"["oc"]~"(V)"))
	graphics.off()

	png(paste(name, "-", x, "-FF.png", sep=""));
	par(mar=mar, mgp = mgp)
	boxplot(FF~droplevels(interaction(rev,get(x))), d, cex.axis=1.5, cex.lab=1.5, border=colors, las=2, ylab="FF")
	graphics.off()
}

single(x=var1, mar=c(8,6,1,1), mgp=c(4,1,0))
if(exists("var2")){single(x=var2, mar=c(8,6,1,1), mgp=c(4,1,0))
if(exists("var3")){single(x=var3, mar=c(8,6,1,1), mgp=c(4,1,0))}}

double <- function(x,y,mar,mgp){
	png(paste(name, "-", x, "-", y, ".png", sep=""));
	par(mar=mar, mgp = mgp)
	boxplot(PCE~droplevels(interaction(rev,get(x),get(y))), d, cex.axis=1.5, cex.lab=1.5, border=colors, las=2, ylab="PCE (%)")
	graphics.off()

	png(paste(name,  "-", x,"-", y, "-Jsc.png", sep=""));
	par(mar=mar, mgp = mgp)
	boxplot(Jsc~droplevels(interaction(rev,get(x),get(y))), d, cex.axis=1.5, cex.lab=1.5, border=colors, las=2, ylab=bquote("J"["sc"]~"(mA/cm"^"2"*")"))
	graphics.off()

	png(paste(name,  "-", x,"-", y, "-Voc.png", sep=""));
	par(mar=mar, mgp = mgp)
	boxplot(Voc~droplevels(interaction(rev,get(x),get(y))), d, cex.axis=1.5, cex.lab=1.5, border=colors, las=2, ylab=bquote("V"["oc"]~"(V)"))
	graphics.off()

	png(paste(name,  "-", x,"-", y, "-FF.png", sep=""));
	par(mar=mar, mgp = mgp)
	boxplot(FF~droplevels(interaction(rev,get(x),get(y))), d, cex.axis=1.5, cex.lab=1.5, border=colors, las=2, ylab="FF")
	graphics.off()
	}

if(exists("var2")){
	double(x=var1, y=var2, mar=c(12,6,1,1), mgp=c(4,1,0))
	if(exists("var3")){
		double(x=var1, y=var3, mar=c(12,6,1,1), mgp=c(4,1,0))
		double(x=var2, y=var3, mar=c(12,6,1,1), mgp=c(4,1,0))
		}
	}
