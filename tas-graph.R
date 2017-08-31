library(plot3D)
theta=-110
phi=40

directory <- tail(strsplit(getwd(), "/")[[1]], n=1)

x=c(read.table("time.table",header=F,colClasses=c("numeric")))
minx=min(x$V1)
maxx=max(x$V1)
y=c(read.table("wavelengths.table",header=F,colClasses=c("numeric")))
maxy=max(y$V1)
miny=min(y$V1)
z=as.matrix(read.table("deltaODmatrix.table",header=F))
maxz=max(c(z))
minz=min(c(z))

png(paste(directory, "_", theta, "_", phi, ".png",sep=""),height=1200, width=1200)
logx=log(x$V1, 10)
logminx=log(minx, 10)
logmaxx=log(maxx, 10)
ribbon3D(logx, y$V1, z, colkey = F, colvar=(exp(z*5000)), ticktype="detailed", xlab="Time (s)", ylab="Probe wavelength (nm)", zlab="delta OD", axes=FALSE, theta=theta, phi=phi, lighting=T, shade=0.2)

prettyx=pretty(logx)
#text3D(x=prettyx, y=rep(miny*0.95, length(prettyx)), z = rep(minz-0.05*(maxz-minz), length(prettyx)), labels=10^prettyx, add=T, cex=1.6)
prettyy=pretty(y$V1)
#text3D(x=rep(logmaxx,length(prettyy)), y=prettyy, z = rep(minz-0.05*(maxz-minz), length(prettyy)), labels=prettyy, add=T, cex=1.6)
prettyz=pretty(z)
#text3D(x=rep(logminx*1.1,length(prettyz)), y=rep(maxy*1.05, length(prettyz)), z = prettyz, labels=prettyz, add=T, cex=1.6)
graphics.off()

#sub1us = length(x$V1[x$V1 < 1e-6])
#x1us=x$V1[1:sub1us]
#z1us=z[1:sub1us,]
#minx1us=min(x1us)
#maxx1us=max(x1us)
#maxz1us=max(c(z1us))
#minz1us=min(c(z1us))
#png(paste(directory, "_", theta, "_", phi, "-1us.png",sep=""),height=1200, width=1200)
#logx1us=log(x1us, 10)
#logminx1us=log(minx1us, 10)
#logmaxx1us=log(maxx1us, 10)
#ribbon3D(logx1us, y$V1, z1us, colkey = F, ticktype="detailed", xlab="Time (s)", ylab="Probe wavelength (nm)", zlab="delta OD", axes=FALSE, theta=theta, phi=phi, lighting=T, shade=0.2)
#prettyx1us=pretty(logx1us)
#text3D(x=prettyx1us, y=rep(miny*0.95, length(prettyx1us)), z = rep(minz1us-0.05*(maxz1us-minz1us), length(prettyx1us)), labels=10^prettyx1us, add=T, cex=1.6)
#prettyy=pretty(y$V1)
#text3D(x=rep(logmaxx1us,length(prettyy)), y=prettyy, z = rep(minz1us-0.05*(maxz1us-minz1us), length(prettyy)), labels=prettyy, add=T, cex=1.6)
#prettyz1us=pretty(z1us)
#text3D(x=rep(logminx1us*1.1,length(prettyz1us)), y=rep(maxy*1.05, length(prettyz1us)), z = prettyz1us, labels=prettyz1us, add=T, cex=1.6)
#graphics.off()

#source("~/software/photophysics-data-processing-R/wavelength_to_rgb.R")
#col=c()
#lapply(y$V1, function(x){col <<- c(col, wavelength_to_rgb(x/2))})
#png(paste(directory, "-1us.png",sep=""), height=1200, width=1200)
#matplot(z1us, col=col, type="l")
#graphics.off()

startp=3
stopp=39
xpoints=x$V1[startp:stopp]
zpoints=z[startp:stopp,]
minxpoints=min(xpoints)
maxxpoints=max(xpoints)
maxzpoints=max(c(zpoints))
minzpoints=min(c(zpoints))
png(paste(directory, "-10_100ns.png",sep=""), height=600, width=1000)
plot(y$V1, colMeans(z), type="l", xlab="Wavelength (nm)", ylab="delta OD", xaxt="n")
axis(side=1, at=pretty(y$V1, n=15), labels=pretty(y$V1, n=15))
graphics.off()


