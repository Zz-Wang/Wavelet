setwd("P:")
rm(list=ls(all=TRUE))
library(spd)
library(MARSS)
library(tseries)
library(rmgarch)
library(copula)
library(QRM)
library(xts)
library(fPortfolio)
library(forecast)
library(VineCopula)
library(mvtnorm)
library(vars)
library(evir)
library(dlm)
library(Rglpk)
library(quadprog)
library(Rsolnp)
library(DEoptim)
library(robustbase)
library(tsDyn)
library(parma)
library(waveslim)
library(hts)
Returns<-read.csv("./Returns (OMXS30 Daily).csv",header=TRUE, sep=";")
DataType<-"OMXS30 Daily"
######################Model#############################
model<-"h.modwt.haar.Wavelet.AR.sGARCH.STD"
########################################################
########################################################
EWL<-100                          # Estimation Window Length
n_level = 2                       #Decomposition Level
h <- 1                            #Forecast Horizon
waveletType = "haar"              #Either haar or la8
waveletBound = "periodic"
waveletMethod = "modwt"
a<-0.1                            #0.1 or 0.05
GarchType<-"sGARCH"               #or gjrGARCH or iGARCH
ARMAOrder<-c(1, 0)                #or c(1,1)
Wave.<-function(x){
  p_modwt <- modwt(x, wf = waveletType, n.level = n_level, boundary = waveletBound)
  p_modwt.m <- matrix(unlist(p_modwt), ncol = n_level+1, byrow = F)
  colnames(p_modwt.m)<-names(p_modwt)
  if (n_level== 1) {
    colnames(p_modwt.m) <- c("D","S")
    y_base <- hts(p_modwt.m, characters = c(1)) # Creating a list containing the base times series
  }
  
  if (n_level== 2) {
    colnames(p_modwt.m) <- c("DD","SD","SS")
    y_base <- hts(p_modwt.m, characters = c(1,1)) # Creating a list containing the base times series
  }
  if (n_level== 3) {
    colnames(p_modwt.m) <- c("DDD","SDD", "SSD","SSS")
    y_base <- hts(p_modwt.m, characters = c(1,1,1)) # Creating a list containing the base times series
  }
  if (n_level== 4) {
    colnames(p_modwt.m) <- c("DDDD","SDDD", "SSDD","SSSD","SSSS")
    y_base <- hts(p_modwt.m, characters = c(1,1,1,1)) 
  }
  if (n_level== 5) {
    colnames(p_modwt.m) <- c("DDDDD","SDDDD", "SSDDD","SSSDD","SSSSD","SSSSS")
    y_base <- hts(p_modwt.m, characters = c(1,1,1,1,1)) 
  }
  if (n_level== 6) {
    colnames(p_modwt.m) <- c("DDDDDD","SDDDDD", "SSDDDD","SSSDDD","SSSSDD","SSSSSD","SSSSSS")
    y_base <- hts(p_modwt.m, characters = c(1,1,1,1,1,1)) 
  }
  if (n_level== 7) {
    colnames(p_modwt.m) <- c("DDDDDDD","SDDDDDD", "SSDDDDD","SSSDDDD","SSSSDDD","SSSSSDD","SSSSSSD","SSSSSSS")
    y_base <- hts(p_modwt.m, characters = c(1,1,1,1,1,1,1)) 
  }
  if (n_level== 8) {
    colnames(p_modwt.m) <- c("DDDDDDDD","SDDDDDDD", "SSDDDDDD","SSSDDDDD","SSSSDDDD","SSSSSDDD","SSSSSSDD","SSSSSSSD","SSSSSSSS")
    y_base <- hts(p_modwt.m, characters = c(1,1,1,1,1,1,1,1)) 
  }
  ally <- allts(y_base)
  list(ally,y_base$nodes)
}

Forecasting.Model<-function(x){
  Wavelet.<-Wave.(x)
  all.series<-Wavelet.[[1]]
  no.<-dim(all.series)[2]
  uspec<-ugarchspec(variance.model = list(model = GarchType, garchOrder = c(1, 
                                                                           1), submodel = NULL, external.regressors = NULL, variance.targeting = FALSE), 
                    mean.model = list(armaOrder = ARMAOrder, include.mean = TRUE, 
                                      archm = FALSE, archpow = 1, arfima = FALSE, external.regressors = NULL, 
                                      archex = FALSE), distribution.model = "std")
  mspec<- multispec( replicate(no.,uspec))
  
  MFit<-multifit(mspec, data = all.series,solver = "hybrid")
  MForecast<-multiforecast(MFit, n.ahead = h, n.roll = 0)
  Mu.hat<- matrix(fitted(MForecast),nr=1)
  S.hat<-matrix(sigma(MForecast),nr=1)
  z<-sapply(1:no.,function(iu) qdist(distribution = "std", a, mu = 0, sigma = 1, lambda = -0.5, 
                                     skew = 1, shape = coef(MFit)['shape',iu]))
  VaR<-Mu.hat+S.hat*z
  Mu.hat.<-aggts(combinef(Mu.hat, Wavelet.[[2]], weights = NULL), levels = 0)
  S.hat.<-aggts(combinef(S.hat, Wavelet.[[2]], weights = NULL), levels = 0)
  VaR.<-aggts(combinef(VaR, Wavelet.[[2]], weights = NULL), levels = 0)
  matrix(cbind(Mu.hat[,1],S.hat[,1],Mu.hat.,S.hat.,VaR[,1],VaR.),nr=1)
}

Forecast<- matrix(NA, ncol = 6, nrow =(3-0))
colnames(Forecast)<-c("Mu","S","WMu","WS","VaR","wVaR")
for(i in 1:nrow(Forecast)){
  cat(i,"",fill = T)
  Forecast[i,]<-Forecasting.Model(Returns[(i+0):(i+0+EWL-1),])
}
#Out.Sample<-Returns[(1+EWL):(dim(Returns )[1]),]
#VaRTest(alpha = a, actual=Out.Sample, VaR=Forecast[,4], conf.level = 0.95)
#plot(Out.Sample,type="l")
#lines(as.matrix(Forecast[,4]),col="red")

#Exceed<-matrix(NA,nr=length(Out.Sample),nc=1)
#for (i in 1:nrow(Exceed)){
#  if(isTRUE(Out.Sample[i]>as.numeric(Forecast[i,4]))==TRUE){Exceed[i,]=0}else{Exceed[i,]=1}
#}

write.csv(Forecast,file=sprintf("Forecast(%s.%d.%d.%d.%s).csv",model,EWL,n_level,a*100,DataType))














