library(lme4)
library(nlme)
nlsListStart <- function(dat){
  out <- nlsList(logeno ~  log(exp(logCaw) + (Ca- exp(logCaw))*exp(-exp(logDaw)/flow))|yr_id,
                 data=dat, start=list(Ca = 2, logCaw = log(68), logDaw = log(12)), # THESE STARTING VALUES ARE IMPORTANT!
                 control=list(tolerance=0.001))
  out
} 

L_TS_NLS_S1<-function(dat,X){
  fitnls <- nlsListStart(dat)
  estCoef <- coef(fitnls)
  estCoef$id<-as.numeric(unlist(lapply(rownames(estCoef),FUN=function(x) strsplit(x,"_")[[1]][2])))
  estCoef$yr<-as.numeric(unlist(lapply(rownames(estCoef),FUN=function(x) strsplit(x,"_")[[1]][1])))
  estCoef_o<-estCoef[order(estCoef$id,estCoef$yr,decreasing=FALSE),]
  estCoef_o<-cbind(estCoef_o,X)
  return(estCoef_o)
}

