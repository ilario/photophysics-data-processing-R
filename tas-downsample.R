logdownsampling <- function(data, s=0.001)
{
  i <- 1
  t <- 0
  means <- c()
  length <- length(data)
  while(i < length)
  {
    by <- floor(10^(t*s))-1
    means <- c(means, mean(data[i:i+by]))
    i <- i+by+1
    t <- t+1
  }
  means=means[!is.na(means)]
  return ( 
      values = means  
  )
}

allfiles = list.files(pattern='*nm.dat')
pngfiles = list.files(pattern='png')
list=allfiles[!allfiles %in% pngfiles] 

print(list[1])
data=read.table(list[1], skip=244, header=F)
time=logdownsampling(data$V1)/10000

wavelength=as.numeric(read.table(list[1], nrows=1, header=F)$V2)
wavelengths=wavelength
negative=read.table(list[1], skip=2, nrows=242, header=F)
bias=mean(negative$V2)
data=read.table(list[1], skip=243, header=F, colClasses=c("numeric","numeric"))
deltaOD=logdownsampling(data$V2)
deltaOD=(deltaOD-bias)/10000
deltaODmatrix=matrix(deltaOD)
list <- list[-1]

invisible(lapply(list, function(x){print(x);
wavelength=as.numeric(read.table(x, nrows=1, header=F)$V2)
wavelengths<<-c(wavelengths, wavelength)

negative=read.table(x, skip=2, nrows=242, header=F)
bias=mean(negative$V2)

data=read.table(x, skip=244, header=F)
deltaOD=logdownsampling(data$V2)
deltaOD=(deltaOD-bias)/10000

deltaODmatrix<<-cbind(deltaODmatrix, deltaOD)
}))

write.table(time, "time.table", row.names=F, col.names=F)
write.table(wavelengths, "wavelengths.table", row.names=F, col.names=F)
write.table(deltaODmatrix, "deltaODmatrix.table", row.names=F, col.names=F)

