library(RColorBrewer)

a <- read.table("table", col.names = c("device","diode","reverse","Jsc","Voc","FF","PCE","comments","notes"),fill = TRUE, na.strings=c("","NA"))
levels(a$notes) <- c(levels(a$notes), "NANA")
a[is.na(a)] <- "NANA"

if(!exists("mycolors")){
  mycolors_once<-brewer.pal(8,"Dark2")
  mycolors<-rep(mycolors_once, each=2)
}

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
  
  #png(paste(name, "-", x, "-PCE.png", sep=""), width=1200, height=600);
  pdf(paste(name, "-", x, "-PCE.pdf", sep=""), width=image_bigpdf_width, height=image_bigpdf_height, pointsize=7);
  par(mar=mar)
  boxplot(PCE~droplevels(interaction(rev,get(x),sep="  ")), d, cex.axis=1.4, border=mycolors, las=2, outline=FALSE, ylim=c(min(d$PCE),max(d$PCE)))
  stripchart(PCE ~ droplevels(interaction(rev,get(x))), vertical = TRUE, data = d, method = "jitter", jitter=0.2, add = TRUE, pch =rep(seq(21,25),each=2), col = add.alpha(mycolors,0.3), cex=1.5)
  mtext(text="PCE (%)", side=2, line=3, cex=1.7, las=1);
  #	title(ylab="PCE (%)", mgp=mgp.y, cex.lab=1.5)
  title(xlab=gsub("_", " ", x), mgp=mgp.x, cex.lab=1.7)
  print("PCE")
  print(aggregate(d$PCE, by=list(interaction(d$rev, d[[x]])), FUN=signif.mean, digits=3))
  print("PCE-stdev")
  print(aggregate(d$PCE, by=list(interaction(d$rev, d[[x]])), FUN=signif.stdev, digits=2))
  graphics.off()
  
  #png(paste(name, "-", x, "-Jsc.png", sep=""), width=1200, height=600);
  pdf(paste(name, "-", x, "-Jsc.pdf", sep=""), width=image_bigpdf_width, height=image_bigpdf_height, pointsize=7);
  par(mar=mar)
  boxplot(Jsc~droplevels(interaction(rev,get(x),sep="  ")), d, cex.axis=1.4, border=mycolors, las=2, outline=FALSE, ylim=c(min(d$Jsc),max(d$Jsc)))
  stripchart(Jsc ~ droplevels(interaction(rev,get(x))), vertical = TRUE, data = d, method = "jitter", jitter=0.2, add = TRUE, pch =rep(seq(21,25),each=2), col = add.alpha(mycolors,0.3), cex=1.5)
  mtext(text=bquote("J"["SC"]~"(mA/cm"^"2"*")"), side=2, line=2, cex=1.7, las=1);
  #	title(ylab=bquote("J"["SC"]~"(mA/cm"^"2"*")"), mgp=mgp.y, cex.lab=1.5)
  title(xlab=gsub("_", " ", x), mgp=mgp.x, cex.lab=1.7)
  print("Jsc")
  print(aggregate(d$Jsc, by=list(interaction(d$rev, d[[x]])), FUN=signif.mean, digits=3))
  print("Jsc-stdev")
  print(aggregate(d$Jsc, by=list(interaction(d$rev, d[[x]])), FUN=signif.stdev, digits=2))
  graphics.off()
  
  #png(paste(name, "-", x, "-Voc.png", sep=""), width=1200, height=600);
  pdf(paste(name, "-", x, "-Voc.pdf", sep=""), width=image_bigpdf_width, height=image_bigpdf_height, pointsize=7);
  par(mar=mar)
  boxplot(Voc~droplevels(interaction(rev,get(x),sep="  ")), d, cex.axis=1.4, border=mycolors, las=2, outline=FALSE, ylim=c(min(d$Voc),max(d$Voc)))
  stripchart(Voc ~ droplevels(interaction(rev,get(x))), vertical = TRUE, data = d, method = "jitter", jitter=0.2, add = TRUE, pch =rep(seq(21,25),each=2), col = add.alpha(mycolors,0.3), cex=1.5)
  mtext(text=bquote("V"["OC"]~"(V)"), side=2, line=4, cex=1.7, las=1);
  #	title(ylab=bquote("V"["OC"]~"(V)"), mgp=mgp.y, cex.lab=1.5)
  title(xlab=gsub("_", " ", x), mgp=mgp.x, cex.lab=1.7)
  print("Voc")
  print(aggregate(d$Voc, by=list(interaction(d$rev, d[[x]])), FUN=signif.mean, digits=3))
  print("Voc-stdev")
  print(aggregate(d$Voc, by=list(interaction(d$rev, d[[x]])), FUN=signif.stdev, digits=2))
  graphics.off()
  
  #png(paste(name, "-", x, "-FF.png", sep=""), width=1200, height=600);
  pdf(paste(name, "-", x, "-FF.pdf", sep=""), width=image_bigpdf_width, height=image_bigpdf_height, pointsize=7);
  par(mar=mar)
  boxplot(FF~droplevels(interaction(rev,get(x),sep="  ")), d, cex.axis=1.4, border=mycolors, las=2, outline=FALSE, ylim=c(min(d$FF),max(d$FF)))
  stripchart(FF ~ droplevels(interaction(rev,get(x))), vertical = TRUE, data = d, method = "jitter", jitter=0.2, add = TRUE, pch = rep(seq(21,25),each=2), col = add.alpha(mycolors,0.3), cex=1.5)
  mtext(text="FF", side=2, line=5, cex=1.7, las=1);
  #	title(ylab="FF", mgp=mgp.y, cex.lab=1.5)
  title(xlab=gsub("_", " ", x), mgp=mgp.x, cex.lab=1.7)
  print("FF")
  print(aggregate(d$FF, by=list(interaction(d$rev, d[[x]])), FUN=signif.mean, digits=3))
  print("FF-stdev")
  print(aggregate(d$FF, by=list(interaction(d$rev, d[[x]])), FUN=signif.stdev, digits=2))
  graphics.off()
}

single(x=var1, mar=c(12,11,1,1), mgp.y=c(4,1,0), mgp.x=c(10,1,0))
if(exists("var2")){single(x=var2, mar=c(12,11,1,1), mgp.y=c(4,1,0), mgp.x=c(10,1,0))
  if(exists("var3")){single(x=var3, mar=c(12,11,1,1), mgp.y=c(4,1,0), mgp.x=c(10,1,0))}}

double <- function(x,y,mar,mgp.y,mgp.x){
  #png(paste(name, "-", x, "-", y, ".png", sep=""), width=1200, height=600);
  pdf(paste(name, "-", x, "-", y, "-PCE.pdf", sep=""), width=image_bigpdf_width, height=image_bigpdf_height, pointsize=7);
  par(mar=mar)
  boxplot(PCE~droplevels(interaction(rev,get(x),get(y),sep="  ")), d, cex.axis=1.4, border=mycolors, las=2, outline=FALSE, ylim=c(min(d$PCE),max(d$PCE)))
  stripchart(PCE ~ droplevels(interaction(rev,get(x),get(y))), vertical = TRUE, data = d, method = "jitter", jitter=0.2, add = TRUE, pch = rep(seq(21,25),each=2), col = add.alpha(mycolors,0.3), cex=1.5)
  mtext(text="PCE (%)", side=2, line=3, cex=1.7, las=1);
  #	title(ylab="PCE (%)", mgp=mgp.y, cex.lab=1.5)
  title(xlab=paste(x,y), mgp=mgp.x, cex.lab=1.7)
  graphics.off()
  
  #png(paste(name,  "-", x,"-", y, "-Jsc.png", sep=""), width=1200, height=600);
  pdf(paste(name, "-", x, "-", y, "-Jsc.pdf", sep=""), width=image_bigpdf_width, height=image_bigpdf_height, pointsize=7);
  par(mar=mar)
  boxplot(Jsc~droplevels(interaction(rev,get(x),get(y),sep="  ")), d, cex.axis=1.4, border=mycolors, las=2, outline=FALSE, ylim=c(min(d$Jsc),max(d$Jsc)))
  stripchart(Jsc ~ droplevels(interaction(rev,get(x),get(y))), vertical = TRUE, data = d, method = "jitter", jitter=0.2, add = TRUE, pch = rep(seq(21,25),each=2), col = add.alpha(mycolors,0.3), cex=1.5)
  mtext(text=bquote("J"["SC"]~"(mA/cm"^"2"*")"), side=2, line=2, cex=1.7, las=1);
  #	title(ylab=bquote("J"["SC"]~"(mA/cm"^"2"*")"), mgp=mgp.y, cex.lab=1.5)
  title(xlab=paste(x,y), mgp=mgp.x, cex.lab=1.7)
  graphics.off()
  
  #png(paste(name,  "-", x,"-", y, "-Voc.png", sep=""), width=1200, height=600);
  pdf(paste(name, "-", x, "-", y, "-Voc.pdf", sep=""), width=image_bigpdf_width, height=image_bigpdf_height, pointsize=7);
  par(mar=mar)
  boxplot(Voc~droplevels(interaction(rev,get(x),get(y),sep="  ")), d, cex.axis=1.4, border=mycolors, las=2, outline=FALSE, ylim=c(min(d$Voc),max(d$Voc)))
  stripchart(Voc ~ droplevels(interaction(rev,get(x),get(y))), vertical = TRUE, data = d, method = "jitter", jitter=0.2, add = TRUE, pch =  rep(seq(21,25),each=2), col = add.alpha(mycolors,0.3), cex=1.5)
  mtext(text=bquote("V"["OC"]~"(V)"), side=2, line=4, cex=1.7, las=1);
  #	title(ylab=bquote("V"["OC"]~"(V)"), mgp=mgp.y, cex.lab=1.5)
  title(xlab=paste(x,y), mgp=mgp.x, cex.lab=1.7)
  graphics.off()
  
  #png(paste(name,  "-", x,"-", y, "-FF.png", sep=""), width=1200, height=600);
  pdf(paste(name, "-", x, "-", y, "-FF.pdf", sep=""), width=image_bigpdf_width, height=image_bigpdf_height, pointsize=7);
  par(mar=mar)
  boxplot(FF~droplevels(interaction(rev,get(x),get(y),sep="  ")), d, cex.axis=1.4, border=mycolors, las=2, outline=FALSE, ylim=c(min(d$FF),max(d$FF)))
  stripchart(FF ~ droplevels(interaction(rev,get(x),get(y))), vertical = TRUE, data = d, method = "jitter", jitter=0.2, add = TRUE, pch = rep(seq(21,25),each=2), col = add.alpha(mycolors,0.3), cex=1.5)
  mtext(text="FF", side=2, line=5, cex=1.7, las=1);
  #	title(ylab="FF", mgp=mgp.y, cex.lab=1.5)
  title(xlab=paste(x,y), mgp=mgp.x, cex.lab=1.7)
  graphics.off()
}

if(exists("var2")){
  double(x=var1, y=var2, mar=c(12,11,1,1), mgp.y=c(4,1,0), mgp.x=c(10,1,0))
  if(exists("var3")){
    double(x=var1, y=var3, mar=c(12,11,1,1), mgp.y=c(4,1,0), mgp.x=c(10,1,0))
    double(x=var2, y=var3, mar=c(12,11,1,1), mgp.y=c(4,1,0), mgp.x=c(10,1,0))
  }
}
