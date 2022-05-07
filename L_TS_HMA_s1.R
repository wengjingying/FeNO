library(MASS)
library(lme4)
library(nlme)
library(reshape2)


source("K:/paper/paper3/HMA_func.R")

L_TS_HMA_S1<-function(dat,X){
  yrid<-unique(dat$yr_id)
  yrs=unique(dat$yr)
  ids=unique(dat$id)
  # Stage I: HMA fits for each person
  estCoef <- as.data.frame(matrix(NA,ncol=7,nrow=length(yrid)))
  colnames(estCoef) <- c("id","yr","Ca","Daw","Jaw","Caw","valid")
  estCoef$id<-as.numeric(unlist(lapply(yrid,FUN=function(x) strsplit(x,"_")[[1]][2])))
  estCoef$yr<-as.numeric(unlist(lapply(yrid,FUN=function(x) strsplit(x,"_")[[1]][1])))
  for(j in 1:length(yrs)){
    for(i in 1:length(ids)){
      dati <- subset(dat, id==ids[i]&yr==yrs[j], select=c(id,flow,eno,pflowlev))
      dati <- subset(dati[order(dati$pflowlev),],pflowlev!=2) # drop flow=50 since HMA uses 3 flow rates (L/M/H)
      avgFlow <- tapply(dati$flow,dati$pflowlev,mean,na.rm=TRUE)
      avgFeNO <- tapply(dati$eno,dati$pflowlev,mean,na.rm=TRUE)
      fL <- avgFlow["1"]
      fM <- avgFlow["3"]
      fH <- avgFlow["4"]  
      eL <- avgFeNO["1"]
      eM <- avgFeNO["3"]
      eH <- avgFeNO["4"]
      # only run Hogman Algorithm if all 3 (low=30,medium=100,high=300) flows are available
      if(any(is.na(c(fL,fM,fH,eL,eM,eH)))){
        HogmanResultsi <- c(Ca=NA,Daw=NA,Jaw=NA,Caw=NA,valid=NA)
      }
      else{
        HogmanResultsi <- Hogman3rdOrderAlg(fL,fM,fH,eL,eM,eH)
        # if algorithm produces any NaN/Inf, replace with all missing (happens 1x in CHS)
        if(any(is.na(HogmanResultsi))){HogmanResultsi <- c(CaNO=NA,DawNO=NA,JawNO=NA,CawNO=NA,valid=NA)}
      }
      estCoef[estCoef$id==ids[i]&estCoef$yr==yrs[j],3:7] <- HogmanResultsi
    }
  }
  
  
  estCoef$logCaw <- log(estCoef$Caw)
  estCoef$logDaw <- log(estCoef$Daw)
  # Stage II: simple linear regressions
  estCoef<-estCoef[order(estCoef$id,estCoef$yr),]
  estCoef<-cbind(estCoef,X)
  return(estCoef)
  }



