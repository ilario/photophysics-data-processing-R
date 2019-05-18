# Copyright (C) 2015-2016 Ilario Gelmetti <iochesonome@gmail.com>
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.

ceFromOutputToGraph <- function(cedir="ce")
{
  print("CE: PLOTTING")
  library(robustbase)
  require(minpack.lm)
  a <- read.table(file.path(cedir, "outputChargeDensityCE.txt"), header=T,stringsAsFactors=F)
  
  #b<-strsplit(a$file, "_")
  #c<-unlist(b)[length(b[[1]])*(1:length(a$file))]
  #d<-as.numeric(gsub("mV", "", c))
  #a$d <- d
  
  expfit <- lm(ChargeDensityCE ~ 0 + Voc, data=a)
  expSuccess = FALSE
  robSuccess = FALSE
  tryCatch({
    expfit <- nlsLM(ChargeDensityCE~ exp(B)*Voc+exp(C)*(exp(exp(D)*Voc)-1), start=list(B=log(max(a$ChargeDensityCE)/max(a$Voc)),C=log(1e-10),D=2), data=a)
    expSuccess = TRUE
  }, error=function(e) print("Failed fit"))
  tryCatch({
    expfit <- nlrob(ChargeDensityCE~ exp(B)*Voc+exp(C)*(exp(exp(D)*Voc)-1), start=list(B=coef(expfit)[[1]],C=coef(expfit)[[2]],D=coef(expfit)[[3]]), data=a)
    robSuccess = TRUE
  }, error=function(e) print("Failed robust fit"))
  if(!robSuccess){
    tryCatch({
      expfit <- nlrob(ChargeDensityCE~ exp(B)*Voc+exp(C)*(exp(exp(D)*Voc)-1), start=list(B=log(max(a$ChargeDensityCE)/max(a$Voc)),C=log(1e-10),D=2), data=a)
      robSuccess = TRUE
    }, error=function(e) print("Failed robust fit"))
  }
  if(!robSuccess){
    tryCatch({
      expfit <- nlrob(ChargeDensityCE~ exp(B)*Voc+exp(C)*(exp(exp(D)*Voc)-1), start=list(B=log(a$ChargeDensityCE[10]/a$Voc[10]),C=log(1e-9),D=3), data=a)
      robSuccess = TRUE
    }, error=function(e) print("Failed robust fit"))
  }
  if(!robSuccess){
    tryCatch({
      expfit <- nlrob(ChargeDensityCE~ exp(B)*Voc+exp(C)*(exp(exp(D)*Voc)-1), start=list(B=log(a$ChargeDensityCE[5]/a$Voc[5]),C=log(1e-11),D=1), data=a)
      robSuccess = TRUE
    }, error=function(e) print("Failed robust fit"))
  }
  
  f <- data.frame(Voc = sort(a$Voc))
  
  directory <- tail(strsplit(getwd(), "/")[[1]], n=2)
  png(file.path(cedir, paste("charge_extraction-", directory[1], ".png", sep="")), width=800, height=800)
  plot(a$Voc, a$ChargeDensityCE, ylab="Charge Density (C/cm2)", xlab="Voltage (V)",cex.lab=1.4, cex.axis=1.4)#, log="y")
  if(exists("expfit")){lines(f$Voc,predict(expfit,f), lwd=2, col="red")}
  graphics.off()
  
  if(exists("expfit")){
    #write.table(t(c("B","Ch0","gamma","GeomCh","ChemCh")), file=file.path(cedir,"outputChargeDensityCE-fit.txt"), append=FALSE, col.names=F, row.names=F);
    write.table(t(c("Cg","neq","m","gamma")), file=file.path(cedir,"outputChargeDensityCE-fit.txt"), append=FALSE, col.names=F, row.names=F);
    Cg=if(expSuccess){exp(coef(expfit)[[1]])} else {coef(expfit)[[1]]}
    neq=if(expSuccess){exp(coef(expfit)[[2]])} else {0}
    gamma=if(expSuccess){exp(coef(expfit)[[3]])} else {Inf}
    # m = q/(gamma kB T)
    m=1.6021766208E-19/(gamma*1.38064852E-23*300)
    # Voc=max(a$Voc)
    # GeomCh=eB*Voc
    # ChemCh=eCh0*(exp(egamma*Voc)-1)
    output <- t(c(Cg, neq, m, gamma))
    write.table(output, file=file.path(cedir,"outputChargeDensityCE-fit.txt"), append=TRUE, col.names=F, row.names=F)
  }
}

