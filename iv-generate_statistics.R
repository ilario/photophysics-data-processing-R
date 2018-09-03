#name="ig"

#var1=""
#v1g1=c()
#v1g2=c()
#v1all=c(v1g1, v1g2)
#list1=list(`` = v1g1, `` = v1g2)

#var2=""
#v2g1=c()
#v2g2=c()#v1all[! v1all %in% v2g1]
#v2all=c(v2g1, v2g2)
#list2=list(`` = v2g1, `` = v2g2)

#var3=""
#v3g1=c()
#v3g2=c()#v1all[! v1all %in% v3g1]
#v3all=c(v3g1, v3g2)
#list3=list(`` = v3g1, `` = v3g2)


a <- read.table("table", col.names = c("device","diode","reverse","Jsc","Voc","FF","PCE","comments","notes"),fill = TRUE, na.strings=c("","NA"))
levels(a$notes) <- c(levels(a$notes), "NANA")
a[is.na(a)] <- "NANA"

colors<-c("blue","blue","darkviolet","darkviolet","darkred","darkred", "darkgreen", "darkgreen")

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

## Add an alpha value to a colour https://gist.github.com/mages/5339689
add.alpha <- function(col, alpha=1){
	if(missing(col))
		stop("Please provide a vector of colours.")
	apply(sapply(col, col2rgb)/255, 2,
		function(x)
			rgb(x[1], x[2], x[3], alpha=alpha))
}

signif.mean <- function(values, digits){
	signif(mean(values), digits=digits)
}
signif.stdev <- function(values, digits){
	signif(sd(values), digits=digits)
}

single <- function(x,mar,mgp.y,mgp.x){
	print(aggregate(d$PCE, by=list(interaction(d$rev, d[[x]])), FUN=length))

	png(paste(name, "-", x, ".png", sep=""), width=1200, height=600);
	par(mar=mar)
	boxplot(PCE~droplevels(interaction(rev,get(x),sep="  ")), d, cex.axis=1.5, border=colors, las=2, outline=FALSE, ylim=c(min(d$PCE),max(d$PCE)))
	stripchart(PCE ~ droplevels(interaction(rev,get(x))), vertical = TRUE, data = d, method = "jitter", add = TRUE, pch = 20, col = add.alpha(colors,0.3), cex=2)
	title(ylab="PCE (%)", mgp=mgp.y, cex.lab=1.5)
	title(xlab=gsub("_", " ", x), mgp=mgp.x, cex.lab=1.5)
	print("PCE")
	print(aggregate(d$PCE, by=list(interaction(d$rev, d[[x]])), FUN=signif.mean, digits=3))
	print("PCE-stdev")
	print(aggregate(d$PCE, by=list(interaction(d$rev, d[[x]])), FUN=signif.stdev, digits=2))
	graphics.off()

	png(paste(name, "-", x, "-Jsc.png", sep=""), width=1200, height=600);
	par(mar=mar)
	boxplot(Jsc~droplevels(interaction(rev,get(x),sep="  ")), d, cex.axis=1.5, border=colors, las=2, outline=FALSE, ylim=c(min(d$Jsc),max(d$Jsc)))
	stripchart(Jsc ~ droplevels(interaction(rev,get(x))), vertical = TRUE, data = d, method = "jitter", add = TRUE, pch = 20, col = add.alpha(colors,0.3), cex=2)
	title(ylab=bquote("J"["sc"]~"(mA/cm"^"2"*")"), mgp=mgp.y, cex.lab=1.5)
	title(xlab=gsub("_", " ", x), mgp=mgp.x, cex.lab=1.5)
	print("Jsc")
	print(aggregate(d$Jsc, by=list(interaction(d$rev, d[[x]])), FUN=signif.mean, digits=3))
	print("Jsc-stdev")
	print(aggregate(d$Jsc, by=list(interaction(d$rev, d[[x]])), FUN=signif.stdev, digits=2))
	graphics.off()

	png(paste(name, "-", x, "-Voc.png", sep=""), width=1200, height=600);
	par(mar=mar)
	boxplot(Voc~droplevels(interaction(rev,get(x),sep="  ")), d, cex.axis=1.5, border=colors, las=2, outline=FALSE, ylim=c(min(d$Voc),max(d$Voc)))
	stripchart(Voc ~ droplevels(interaction(rev,get(x))), vertical = TRUE, data = d, method = "jitter", add = TRUE, pch = 20, col = add.alpha(colors,0.3), cex=2)
	title(ylab=bquote("V"["oc"]~"(V)"), mgp=mgp.y, cex.lab=1.5)
	title(xlab=gsub("_", " ", x), mgp=mgp.x, cex.lab=1.5)
	print("Voc")
	print(aggregate(d$Voc, by=list(interaction(d$rev, d[[x]])), FUN=signif.mean, digits=3))
	print("Voc-stdev")
	print(aggregate(d$Voc, by=list(interaction(d$rev, d[[x]])), FUN=signif.stdev, digits=2))
	graphics.off()

	png(paste(name, "-", x, "-FF.png", sep=""), width=1200, height=600);
	par(mar=mar)
	boxplot(FF~droplevels(interaction(rev,get(x),sep="  ")), d, cex.axis=1.5, border=colors, las=2, outline=FALSE, ylim=c(min(d$FF),max(d$FF)))
	stripchart(FF ~ droplevels(interaction(rev,get(x))), vertical = TRUE, data = d, method = "jitter", add = TRUE, pch = 20, col = add.alpha(colors,0.3), cex=2)
	title(ylab="FF", mgp=mgp.y, cex.lab=1.5)
	title(xlab=gsub("_", " ", x), mgp=mgp.x, cex.lab=1.5)
	print("FF")
	print(aggregate(d$FF, by=list(interaction(d$rev, d[[x]])), FUN=signif.mean, digits=3))
	print("FF-stdev")
	print(aggregate(d$FF, by=list(interaction(d$rev, d[[x]])), FUN=signif.stdev, digits=2))
	graphics.off()
}

single(x=var1, mar=c(12,6,1,1), mgp.y=c(4,1,0), mgp.x=c(8,1,0))
if(exists("var2")){single(x=var2, mar=c(12,6,1,1), mgp.y=c(4,1,0), mgp.x=c(8,1,0))
if(exists("var3")){single(x=var3, mar=c(12,6,1,1), mgp.y=c(4,1,0), mgp.x=c(8,1,0))}}

double <- function(x,y,mar,mgp.y,mgp.x){
	png(paste(name, "-", x, "-", y, ".png", sep=""), width=1200, height=600);
	par(mar=mar)
	boxplot(PCE~droplevels(interaction(rev,get(x),get(y),sep="  ")), d, cex.axis=1.5, border=colors, las=2, outline=FALSE, ylim=c(min(d$PCE),max(d$PCE)))
	stripchart(PCE ~ droplevels(interaction(rev,get(x),get(y))), vertical = TRUE, data = d, method = "jitter", add = TRUE, pch = 20, col = add.alpha(colors,0.3), cex=2)
	title(ylab="PCE (%)", mgp=mgp.y, cex.lab=1.5)
	title(xlab=paste(x,y), mgp=mgp.x, cex.lab=1.5)
	graphics.off()

	png(paste(name,  "-", x,"-", y, "-Jsc.png", sep=""), width=1200, height=600);
	par(mar=mar)
	boxplot(Jsc~droplevels(interaction(rev,get(x),get(y),sep="  ")), d, cex.axis=1.5, border=colors, las=2, outline=FALSE, ylim=c(min(d$Jsc),max(d$Jsc)))
	stripchart(Jsc ~ droplevels(interaction(rev,get(x),get(y))), vertical = TRUE, data = d, method = "jitter", add = TRUE, pch = 20, col = add.alpha(colors,0.3), cex=2)
	title(ylab=bquote("J"["sc"]~"(mA/cm"^"2"*")"), mgp=mgp.y, cex.lab=1.5)
	title(xlab=paste(x,y), mgp=mgp.x, cex.lab=1.5)
	graphics.off()

	png(paste(name,  "-", x,"-", y, "-Voc.png", sep=""), width=1200, height=600);
	par(mar=mar)
	boxplot(Voc~droplevels(interaction(rev,get(x),get(y),sep="  ")), d, cex.axis=1.5, border=colors, las=2, outline=FALSE, ylim=c(min(d$Voc),max(d$Voc)))
	stripchart(Voc ~ droplevels(interaction(rev,get(x),get(y))), vertical = TRUE, data = d, method = "jitter", add = TRUE, pch = 20, col = add.alpha(colors,0.3), cex=2)
	title(ylab=bquote("V"["oc"]~"(V)"), mgp=mgp.y, cex.lab=1.5)
	title(xlab=paste(x,y), mgp=mgp.x, cex.lab=1.5)
	graphics.off()

	png(paste(name,  "-", x,"-", y, "-FF.png", sep=""), width=1200, height=600);
	par(mar=mar)
	boxplot(FF~droplevels(interaction(rev,get(x),get(y),sep="  ")), d, cex.axis=1.5, border=colors, las=2, outline=FALSE, ylim=c(min(d$FF),max(d$FF)))
	stripchart(FF ~ droplevels(interaction(rev,get(x),get(y))), vertical = TRUE, data = d, method = "jitter", add = TRUE, pch = 20, col = add.alpha(colors,0.3), cex=2)
	title(ylab="FF", mgp=mgp.y, cex.lab=1.5)
	title(xlab=paste(x,y), mgp=mgp.x, cex.lab=1.5)
	graphics.off()
	}

if(exists("var2")){
	double(x=var1, y=var2, mar=c(12,6,1,1), mgp.y=c(4,1,0), mgp.x=c(10,1,0))
	if(exists("var3")){
		double(x=var1, y=var3, mar=c(12,6,1,1), mgp.y=c(4,1,0), mgp.x=c(10,1,0))
		double(x=var2, y=var3, mar=c(12,6,1,1), mgp.y=c(4,1,0), mgp.x=c(10,1,0))
		}
	}
