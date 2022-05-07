library(lme4)
library(nlme)
nlsListStart <- function(dat){
  out <- nlsList(logeno ~  log(exp(logCaw) + (Ca- exp(logCaw))*exp(-exp(logDaw)/flow))|id,
                 data=dat, start=list(Ca = 2, logCaw = log(68), logDaw = log(12)), # THESE STARTING VALUES ARE IMPORTANT!
                 control=list(tolerance=0.001))
  out
} 

TS_NLS_S1<-function(dat,X){
  fitnls <- nlsListStart(dat)
  estCoef <- coef(fitnls)
  estCoef$id<-as.numeric(rownames(estCoef))
  estCoef_o<-estCoef[order(estCoef$id,decreasing=FALSE),]
  estCoef_o<-cbind(estCoef_o,X)
  return(estCoef_o)
}

