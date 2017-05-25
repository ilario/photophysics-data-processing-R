#name="ig"

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

file.create("output.txt")
mydata <- import.iv.separated(pattern.excl="\\.png$", pattern="^ig.*.txt$", list.excl="output.txt")
results <- lapply(names(mydata), function(x){
z <- factor(strsplit(x,"-")[[1]][2], v1all); 
levels(z) <- list1;
if(exists("var2")){
	y <- factor(strsplit(x,"-")[[1]][2], v2all);
	levels(y) <- list2;
	comment=paste(toString(y), toString(z))}
else {comment=toString(z)}


write(extract.iv(mydata[[x]]$Voltage_V, mydata[[x]]$Current_mA, formatted.output=TRUE, directory=getwd(), sample=x, reverse=as.integer(grepl("reverse", x)), comment=comment), file="output.txt", append=TRUE);
extract.iv(mydata[[x]]$Voltage_V, mydata[[x]]$Current_mA, comment=comment)});
names(results) <- names(mydata)


