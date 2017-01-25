v1g1=c("")
v1g2=c("")
v1g3=c("")
v1all=c(v1g1, v1g2, v1g3)

v2g1=c("")
v2g2=v1all[! v1all %in% v2g1]
v2all=v1all

file.create("output.txt")
mydata <- import.iv.separated(pattern.excl="\\.png$", pattern="ig.*.txt$", list.excl="output.txt")
results <- lapply(names(mydata), function(x){
z <- factor(strsplit(x,"-")[[1]][2], v1all); 
levels(z) <- list(`` = v1g1, `` = v1g2, `` = v1g3);
#y <- factor(strsplit(x,"-")[[1]][2], v2all);
#levels(y) <- list(`` = v2g1, `` = v2g2);
y <- ""

write(extract.iv(mydata[[x]]$Voltage_V, mydata[[x]]$Current_mA, formatted.output=TRUE, directory=getwd(), sample=x, reverse=as.integer(grepl("reverse", x)), comment=paste(toString(y), toString(z))), file="output.txt", append=TRUE);
extract.iv(mydata[[x]]$Voltage_V, mydata[[x]]$Current_mA, comment=paste(toString(y), toString(z)))});
names(results) <- names(mydata)


