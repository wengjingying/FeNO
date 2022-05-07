library(MASS)
library(lme4)
library(nlme)
library(reshape2)

nonLinInitPlogJ <-function(mCall,LHS,data){
  xy <- sortedXyData(mCall[["predictor"]],LHS,data)
  fit <- lm(xy[,"y"]~I(1/xy[,"x"]) + I(1/xy[,"x"]^2))
  est <- fit$coef
  CaNO <- est[1]
  JawNO <- est[2] - 2*est[1]*est[3]/est[2]
  DawNO <- -2*est[3]/est[2]
  logDawNO <- log(max(DawNO,0.0001))
  logJawNO <- log(max(JawNO,0.0001))
  logCawNO <- log(max(JawNO/DawNO,0.0001))
  value <- c(logCawNO,CaNO,logDawNO)
  names(value) <- mCall[c("logCawNO","CaNO","logDawNO")]
  value
}
## quadT self-starter prelim function
nonLinInitTlogJ <-function(mCall,LHS,data){
  xy <- sortedXyData(mCall[["predictor"]],LHS,data)
  fit <- lm(I(xy[,"y"]*xy[,"x"])~I(xy[,"x"]) + I(1/xy[,"x"]))
  est <- fit$coef
  CaNO <- est[2]
  JawNO <- est[1] - 2*est[2]*est[3]/est[1]
  DawNO <- -2*est[3]/est[1]
  logDawNO <- log(max(DawNO,0.0001))
  logJawNO <- log(max(JawNO,0.0001))
  logCawNO <- log(max(JawNO/DawNO,0.0001))
  value <- c(logCawNO,CaNO,logDawNO)
  names(value) <- mCall[c("logCawNO","CaNO","logDawNO")]
  value
}

## model without gradient
nonLinLogModellogJ <- function(predictor,logCawNO,CaNO,logDawNO){
  log(exp(logCawNO) + (CaNO - exp(logCawNO))*exp(-exp(logDawNO)/predictor))
}
# model with gradient
nonLinLogModellogJwithGrad <- deriv( ~log(exp(logCawNO) + (CaNO - exp(logCawNO))*exp(-exp(logDawNO)/predictor)),
                                     c("logCawNO","CaNO","DawNO"),
                                     function(predictor,logCawNO,CaNO,logDawNO){})

#self-starting function with gradient
SSnonLinLogPlogJGrad <- selfStart(nonLinLogModellogJwithGrad, nonLinInitPlogJ, c("logCawNO","CaNO","logDawNO"))
#self-starting function without gradient
SSnonLinLogPlogJwoGrad <- selfStart(nonLinLogModellogJ, nonLinInitPlogJ, c("logCawNO","CaNO","logDawNO"))

L_TS_NLME_S1<-function(dat,tol1,tol2){
  mydat <- subset(dat,select=c(yr_id,logeno,flow,id,yr)) #,id%in%unique(dat$id)[1:50]
  dat2 <- groupedData(logeno~ flow | id/yr, data=mydat)
  # nls separately for each group
  fitnlsList <- nlsList(logeno~SSnonLinLogPlogJwoGrad(flow,logCawNO,CaNO,logDawNO),data=dat2,
                        control=list(tolerance = 0.1))
  # starting values for fixed effects
  cf <- coef(fitnlsList)
  startfix <- unlist(lapply(cf, median, na.rm = TRUE))
  # starting values for id-level RE
  idREstart <- matrix(0,nrow=length(unique(mydat$id)),ncol=3)
  colnames(idREstart) <- names(startfix); rownames(idREstart) <- unique(mydat$id)
  idvisitmat <- matrix(as.integer(unlist(strsplit(rownames(coef(fitnlsList)),"/"))),byrow=T,ncol=2)
  for(i in 1:3){
    fillcol <- tapply(cf[,i],as.factor(idvisitmat[,1]),mean,na.rm=T)
    idREstart[which(!is.nan(fillcol)),i] <- fillcol[which(!is.nan(fillcol))]-startfix[i]
  }
  # starting values for visit(within-id)-level RE
  visitREstart <- matrix(0,nrow=nrow(cf),ncol=3)
  colnames(visitREstart) <- names(startfix); rownames(visitREstart) <- rownames(cf)
  visitREstartIDmeans <- visitREstart
  for(i in unique(mydat$id)){
    myrow=which(idvisitmat[,1]==i)
    visitREstartIDmeans[myrow,] <- matrix(rep(idREstart[which(rownames(idREstart)==i),] + startfix,
                                              length(myrow)),
                                          ncol=3,byrow=TRUE)
  }
  visitREstart <- cf-visitREstartIDmeans
  visitREstart[is.na(visitREstart)] <- 0
  visitREstart <- as.matrix(visitREstart)
  


  fit1.nlme <- nlme(logeno~SSnonLinLogPlogJwoGrad(flow,logCawNO,CaNO,logDawNO),
                    data=mydat,
                    fixed=logCawNO + CaNO +logDawNO ~ 1,
                    random= list(pdLogChol(logCawNO + CaNO +logDawNO ~ 1),
                                 pdDiag(logCawNO + CaNO +logDawNO ~ 1)),
                    groups=~id/yr,
                    start=list( fixed=c(1.82,4.12,2.63),
                                random=list(id=idREstart, # random effects starting values
                                            yr=visitREstart) # only speed up cvg by ~5% overall
                                #  ~3.5 hours with both visits CHS data
                    ),
                    verbose=T,
                    control=list(tolerance = tol1,opt="nlminb",natural=TRUE,niterEM=50))
  fit2.nlme <- update(fit1.nlme,control=list(tolerance =tol2,maxIter=50),verbose=T)
  return(fit2.nlme)
}

L_TS_NLME_S2<-function(dat,fit,X){
  mydat <- subset(dat,select=c(yr_id,logeno,flow,id,yr)) #,id%in%unique(dat$id)[1:50]
  nlmeFit_yrinid_logCawlogDawCa <- fit
  # create a matrix with predicted id-level RE, repeated for each year
  idRE <- ranef(fit,level=1)[as.character(unique(mydat$id)),]
  idREyr <- matrix(0,nrow=nrow(cf),ncol=3)
  colnames(idREyr) <- names(startfix); rownames(idREyr) <- rownames(cf)
  for(i in unique(mydat$id)){
    myrow <- which(idvisitmat[,1]==i)
    idREyr[myrow,] <- matrix(rep(as.numeric(idRE[which(rownames(idRE)==i),]),
                                 length(myrow)),
                             ncol=3,byrow=TRUE)
  }
  nlmeFit_yrinid_logCawlogDawCa_CaNO_est      <- idREyr[rownames(cf),"CaNO"] + ranef(fit,level=2)[rownames(cf),"CaNO"] + fixef(fit)["CaNO"]
  nlmeFit_yrinid_logCawlogDawCa_logCawNO_est  <- idREyr[rownames(cf),"logCawNO"] + ranef(fit,level=2)[rownames(cf),"logCawNO"] + fixef(fit)["logCawNO"]
  nlmeFit_yrinid_logCawlogDawCa_logDawNO_est  <- idREyr[rownames(cf),"logDawNO"] + ranef(fit,level=2)[rownames(cf),"logDawNO"] + fixef(fit)["logDawNO"]
  nlmeFit_est<-cbind(nlmeFit_yrinid_logCawlogDawCa_CaNO_est,cbind(nlmeFit_yrinid_logCawlogDawCa_logCawNO_est,nlmeFit_yrinid_logCawlogDawCa_logDawNO_est))
  colnames(nlmeFit_est)<-c("Ca","logCaw","logDaw")
  nlmeFit_est<-as.data.frame(nlmeFit_est)
  nlmeFit_est$id<-t(as.data.frame(strsplit(rownames(nlmeFit_est),split="/")))[,1]
  nlmeFit_est$yr<-t(as.data.frame(strsplit(rownames(nlmeFit_est),split="/")))[,2]
  nlmeFit_est<-nlmeFit_est[order(nlmeFit_est$id,nlmeFit_est$yr),]
  nlmeFit_est<-cbind(nlmeFit_est,X)
  return(nlmeFit_est)
}

