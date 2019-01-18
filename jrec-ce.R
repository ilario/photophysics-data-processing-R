# Copyright (C) 2015-2017 Ilario Gelmetti <iochesonome@gmail.com>
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

jrecCe <- function(cedir="ce", tpvdir="tpv")
{
  print("Jrec-CE: CALCULATING")
  ce <- read.table(file.path(cedir, "outputChargeDensityCE.txt"), header=T,stringsAsFactors=F)
  cefit <- read.table(file.path(cedir,"outputChargeDensityCE-fit.txt"), header=T,stringsAsFactors=F)
  tpv <- read.table(file.path(tpvdir,"output-monoexp.txt"), header=T,stringsAsFactors=F)
  tpvfit <- read.table(file.path(tpvdir,"output-monoexp-fit.txt"), header=T,stringsAsFactors=F)
  
  directory <- tail(strsplit(getwd(), "/")[[1]], n=1)
  
  write.table(t(c("Jrec")), file=paste("jrec-ce-",directory[1],".txt",sep=""), append=FALSE, col.names=F, row.names=F);
  ce1sun <- tail(ce$ChargeDensityCE[order(ce$Voc)], n=1)
  print(paste("CE 1sun",ce1sun))
  gamma <- cefit$gamma
  print(paste("gamma", gamma))
  tau1sun <- tail(tpv$T[order(tpv$Voc)], n=1)
  print(paste("Tau 1sun", tau1sun))
  beta <- -tpvfit$beta
  print(paste("beta", beta))
  output <- ce1sun/((1+beta/gamma)*tau1sun)
  print(paste("Jrec-CE", output))
  write.table(output, file=paste("jrec-ce-",directory[1],".txt",sep=""), append=TRUE, col.names=F, row.names=F)
}
