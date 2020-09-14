library(MASS)
library(reshape2)
library(lme4)
library(nlme)
library(coda)
library(parallel)
library(R2jags)


alpha   <-c(1,3.5,2.5)
# maneuver flow
flow    <- c(rep(30,2),rep(50,2),rep(100,2),rep(300,2))
# measurement error
sderror <-0.1

# Variance-Covariance matrix for NO parameters
sdCa            <- 0.66392040
sdlogCaw        <- 0.78668940
sdlogDaw        <- 0.59812968

corlogCawCa     <- 0.255
corlogCawlogDaw <- -0.311
corlogDawCa     <- -0.357
CovMat_chs_Caw<- matrix(c(sdCa^2,corlogCawCa*sdCa*sdlogCaw,corlogDawCa*sdCa*sdlogDaw,
                          corlogCawCa*sdCa*sdlogCaw,sdlogCaw^2,corlogCawlogDaw*sdlogCaw*sdlogDaw,
                          corlogDawCa*sdCa*sdlogDaw,corlogCawlogDaw*sdlogCaw*sdlogDaw,sdlogDaw^2),3,3)

#function to generate datasets
Obs_FB<-function(alpha,beta,X,Flow,SD,NOcov){
  #set empty matrixs for logeno and NO parameters
  logeno <- matrix(NA,ncol=length(Flow),nrow=length(X))
  colnames(logeno)<-Flow
  NOparam <-matrix(NA,ncol=3,nrow=length(X))
  colnames(NOparam)<-c("Ca","logCaw","logDaw")
  
  for(i in 1:length(X)){
    repeat{
      # NOparam ~ MVN (mu(Ca, logCaw,logDaw),CovMat_chs)
      param <- mvrnorm(1, X[i] *  beta + alpha, NOcov)
      FeNOi  <- exp(param[2])+(param[1]-exp(param[2]))*exp(-1*exp(param[3])/Flow)
      # FeNO shoule be positive
      # log(FeNO)~N(0,SD)
      logFeNOobsi  <- ifelse(FeNOi>0,log(FeNOi) + rnorm(length(Flow),0,SD),NA) 
      # Ca should also be positive, met both constrains will end the loop
      if(!any(is.na(logFeNOobsi)) & param[1]>0){break}#satisfy all
    }
    logeno[i,]  <- logFeNOobsi
    NOparam[i,] <- c(param[1],param[2],param[3])
  }
  return(list(NOparam=NOparam,logeno=logeno,X=X))
}

# environment factor, assumed to be normal scaled
X<-rnorm(Ndat,0,1)
# Scenario 1
beta<-c(0.1,0.1,0.1)
# The populational NO parameter level. For example, Ca_i = alpha_Ca + beta_Ca * X_i




#Jags model
truncationNOX <- function()  {
  for( i in 1 : Ndat) {
    for( j in 1 : Nflows) {
      logeno[i,j] ~ dnorm(logmu[i,j],tau_c)
      mu[i,j] <- exp(logCaw[i])+(Ca[i]-exp(logCaw[i]))*exp(-1*exp(logDaw[i])/flow[j])
      logmu[i,j]<-log(max(mu[i,j] ,0.001))#logFeNO>0.001
    }
    
    Ca[i]       <- param_subject[i, 1]
    logCaw[i]   <- param_subject[i, 2] 
    logDaw[i]   <- param_subject[i, 3] 
    
    param_subject[i,1] ~ dnorm(mu1_22[i],tau1_22);T(0,)
    mu1_22[i]<-Ca_mean[i]+part1[1,1]*(param_subject[i,2]-logCaw_mean[i])+part1[1,2]*(param_subject[i,3]-logDaw_mean[i])
    param_subject[i,2:3] ~ dmnorm(c(logCaw_mean[i],logDaw_mean[i]),tau22[1:2,1:2])
    
    Ca_mean[i]<-alpha_Ca+beta_Ca * X[i]
    logCaw_mean[i]<-alpha_logCaw+beta_logCaw * X[i]
    logDaw_mean[i]<-alpha_logDaw+beta_logDaw * X[i]
  }
  
  tau1_22<-1/cov1_22
  ############# variance of A Conditioned on B
  cov1_22[1,1]<-varCa-(covCalogCaw^2*tau22[1,1]+2*covCalogDaw*tau22[2,1]*covCalogCaw+covCalogDaw^2*tau22[2,2])
  ############# SigmaAB * tauB
  part1[1,1]<-covCalogCaw*tau22[1,1]+covCalogDaw*tau22[2,1]
  part1[1,2]<-covCalogCaw*tau22[1,2]+covCalogDaw*tau22[2,2]
  ############# Sigma B
  tau22[1:2,1:2]<-inverse(cov22[1:2,1:2])
  cov22[1,1] <- varlogCaw
  cov22[1,2] <- covlogCawlogDaw
  cov22[2,1] <- covlogCawlogDaw
  cov22[2,2] <- varlogDaw
  ############# var-cov #SPE: change name from rho to cor for consistency
  covCalogCaw       <-corlogCawCa*sqrt(varCa*varlogCaw)
  covCalogDaw       <-corlogDawCa*sqrt(varCa*varlogDaw)
  covlogCawlogDaw   <-corlogCawlogDaw*sqrt(varlogCaw*varlogDaw)
  
  varCa     <-1/tauCa
  varlogCaw <-1/taulogCaw
  varlogDaw <-1/taulogDaw
  
  ############ SPE: calculate more interpretable values for monitoring
  sdCa     <- sqrt(varCa)
  sdlogCaw <- sqrt(varlogCaw)
  sdlogDaw <- sqrt(varlogDaw)
  sigma_c  <- sqrt(1/tau_c)

  ############# sample correlation
  corlogDawCa   ~ dunif(L_corlogDawCa,U_corlogDawCa)
  L_corlogDawCa <- corlogCawCa*corlogCawlogDaw-sqrt((corlogCawlogDaw^2-1)*(corlogCawCa^2-1))
  U_corlogDawCa <- corlogCawCa*corlogCawlogDaw+sqrt((corlogCawlogDaw^2-1)*(corlogCawCa^2-1))
  
  corlogCawCa        ~ dunif(-1,1)
  corlogCawlogDaw    ~ dunif(-1,1)
  ############# sample tau
  tauCa     ~ dgamma(1.0E-3, 1.0E-3)
  taulogCaw ~ dgamma(1.0E-3, 1.0E-3)
  taulogDaw ~ dgamma(1.0E-3, 1.0E-3)
  tau_c     ~ dgamma(1.0E-3, 1.0E-3) #random error in regression
  ############# sample regression coefficients (mean)
  alpha_Ca     ~ dnorm(alphaC_prior[1],1e-6)
  alpha_logCaw ~ dnorm(alphaC_prior[2],1e-6)
  alpha_logDaw ~ dnorm(alphaC_prior[3],1e-6)
  beta_Ca      ~ dnorm(betaC_prior[1], 1e-6)
  beta_logCaw  ~ dnorm(betaC_prior[2], 1e-6)
  beta_logDaw  ~ dnorm(betaC_prior[3], 1e-6)
} 


###### TS-NLS #####

# the dataset shape is different than we needed for U-HB

# NLS function
nlsListStart <- function(dat){
  out <- nlsList(logeno ~  log(exp(logCaw) + (Ca- exp(logCaw))*exp(-exp(logDaw)/flow))|id,
                 data=dat, start=list(Ca = 2, logCaw = log(68), logDaw = log(12)), # THESE STARTING VALUES ARE IMPORTANT!
                 control=list(tolerance=0.001))
  out
} 

nls_FINAL <- function(mynlsfit,X){
  # Stage I: from nlsList, separate nonlinLog models for each person
  estCoef <- coef(mynlsfit)
  # Stage II: simple linear regressions
  fit_Ca     <- lm(estCoef$Ca~X)
  fit_logCaw <- lm(estCoef$logCaw~X)
  fit_logDaw <- lm(estCoef$logDaw~X)
  # output
  out <- matrix(NA,nrow=6,ncol=5)
  rownames(out) <- c("Ca.(Intercept)","Ca.X","logCaw.(Intercept)","logCaw.X","logDaw.(Intercept)","logDaw.X")
  colnames(out) <- c("est","lb","ub","se","p")
  out <- as.data.frame(out)
  
  out[c("Ca.(Intercept)","Ca.X"),c("est","se","p")] <- summary(fit_Ca)$coef[,c("Estimate","Std. Error","Pr(>|t|)")]
  out[c("logCaw.(Intercept)","logCaw.X"),c("est","se","p")] <- summary(fit_logCaw)$coef[,c("Estimate","Std. Error","Pr(>|t|)")]
  out[c("logDaw.(Intercept)","logDaw.X"),c("est","se","p")] <- summary(fit_logDaw)$coef[,c("Estimate","Std. Error","Pr(>|t|)")]
  
  out[c("Ca.(Intercept)","Ca.X"),c("lb","ub")] <- confint(fit_Ca)[,c("2.5 %","97.5 %")]
  out[c("logCaw.(Intercept)","logCaw.X"),c("lb","ub")] <- confint(fit_logCaw)[,c("2.5 %","97.5 %")]
  out[c("logDaw.(Intercept)","logDaw.X"),c("lb","ub")] <- confint(fit_logDaw)[,c("2.5 %","97.5 %")]
  
  out
}


######## TS-HMA ##########

Hogman3rdOrderAlg <- function(fL,fM,fH,eL,eM,eH){  
  # linT approximation on high flows to get S=CaNO and I
  CaNO_est <- S <- (fH*eH - fM*eM)/(fH-fM)
  # solve for I using linT model at medium flow rate (flow=100)
  I <- (fM*eM - S*fM)
  JawNO_est  <- I
  
  # first order approximation starting value
  Dw0_fo <- I/(eL - S)
  
  # iterative algorithm function, 1st order
  italg <- function(SV,nIterAlg=10){
    # object to store Dw values for each iteration
    outIterAlg <- numeric(nIterAlg)
    # set starting value (SV)
    Dw_n <- SV
    for(i in 1:nIterAlg){
      Dw_nplus1 <- outIterAlg[i] <-  -I*(exp(-Dw_n/fL) - exp(-Dw_n/fM))/(eL-eM) # Inf-NAN-error
      diffDw <- abs(Dw_n - Dw_nplus1)
      #print(i);print(Dw_nplus1)
      # stop if values of Dw_n and Dw_nplus1 are close enough
      if(diffDw <0.001){break}
      Dw_n <- Dw_nplus1  
    }
    out <- list(DawNO=Dw_nplus1,nIter=i,diffDaw=diffDw)
    out
  }
  DawNO_result_fo <- italg(Dw0_fo)
  DawNO_est_fo <- DawNO_result_fo$DawNO
  
  CawNO_est_fo <- (I + S*DawNO_est_fo)/DawNO_est_fo
  
  Ic <- I/(1-(fM+fH)*DawNO_est_fo/(2*(fM*fH)+(fH^3-fM^3)/(6*fM^2*fH^2*(fH-fM)/DawNO_est_fo^2)))
  Sc <- S-DawNO_est_fo*Ic/2/fM/fH+(fM+fH)*DawNO_est_fo^2*Ic/6/fM^2/fH^2
  
  # iterative algorithm function, 3rd order
  italgto <- function(SV,nIterAlg=10){
    # object to store Dw values for each iteration
    outIterAlg <- numeric(nIterAlg)
    # set starting value (SV)
    Dw_n <- SV
    for(i in 1:nIterAlg){
      Dw_nplus1 <- outIterAlg[i] <-  -Ic*(exp(-Dw_n/fL) - exp(-Dw_n/fM))/(eL-eM)
      diffDw <- abs(Dw_n - Dw_nplus1)
      #print(i);print(Dw_nplus1)
      # stop if values of Dw_n and Dw_nplus1 are close enough
      if(diffDw <0.001){break}
      Dw_n <- Dw_nplus1  
    }
    out <- list(DawNO=Dw_nplus1,nIter=i,diffDaw=diffDw)
    out
  }
  DawNO_result_to <- italgto(DawNO_est_fo)
  DawNO_est_to <- DawNO_result_to$DawNO
  
  DawNO_est_final <- DawNO_est_to  
  CawNO_est_final <- Ic/DawNO_est_final + Sc
  CaNO_est_final <- Sc
  JawNO_est_final <- DawNO_est_final*(CawNO_est_final - CaNO_est_final)
  
  outfinal <- c(  CaNO=unname(CaNO_est_final),
                  DawNO=unname(DawNO_est_final),
                  JawNO=unname(JawNO_est_final),
                  CawNO=unname(CawNO_est_final)
  )
  # data checks (both should be true)
  ## test CaNO estimate is positive
  test1 <- CaNO_est_final > 0
  ## test whether measured data is mathematically consistent with model
  LHS <- (eL-eM)/(eM-eH)
  RHS <- fH/fL*((fM-fL)/(fH-fM))
  test2 <- LHS < RHS
  ## To change the function to only give non-missing values for 
  ## valid datasets that produce positive CaNO estimates, remove the 
  ## comments from the following two lines
  #if(!(test2)) outfinal[1:4] <- NA # leave in negative CaNO estimates
  #if(!(test1 & test2)) outfinal[1:4] <- NA
  c(outfinal,valid=unname(test2))
}

HMA_FINAL <- function(dat,X){
  ids <- unique(dat$id)
  # Stage I: HMA fits for each person
  estCoef <- matrix(NA,ncol=5,nrow=length(ids))
  rownames(estCoef) <- ids
  colnames(estCoef) <- c("Ca","Daw","Jaw","Caw","valid")
  for(i in 1:length(ids)){
    dati <- subset(dat, id==ids[i], select=c(id,flow,eno))
    dati <- subset(dati[order(dati$flow),],flow!=50) # drop flow=50 since HMA uses 3 flow rates (L/M/H)
    avgFeNO <- tapply(dati$eno,dati$flow,mean,na.rm=TRUE)
    fL <- 30
    fM <- 100
    fH <- 300
    eL <- avgFeNO["30"]
    eM <- avgFeNO["100"]
    eH <- avgFeNO["300"]
    # only run Hogman Algorithm if all 3 (low=30,medium=100,high=300) flows are available
    if(any(is.na(c(fL,fM,fH,eL,eM,eH)))){
      HogmanResultsi <- c(Ca=NA,Daw=NA,Jaw=NA,Caw=NA,valid=NA)
    }
    else{
      HogmanResultsi <- Hogman3rdOrderAlg(fL,fM,fH,eL,eM,eH)
      # if algorithm produces any NaN/Inf, replace with all missing (happens 1x in CHS)
      if(any(is.na(HogmanResultsi))){HogmanResultsi <- c(CaNO=NA,DawNO=NA,JawNO=NA,CawNO=NA,valid=NA)}
    }
    estCoef[i,] <- HogmanResultsi
  }
  estCoef <- as.data.frame(estCoef)
  estCoef$logCaw <- log(estCoef$Caw)
  estCoef$logDaw <- log(estCoef$Daw)
  # Stage II: simple linear regressions
  fit_Ca     <- lm(estCoef$Ca~X)
  fit_logCaw <- lm(estCoef$logCaw~X)
  fit_logDaw <- lm(estCoef$logDaw~X)
  # output
  out <- matrix(NA,nrow=6,ncol=5)
  rownames(out) <- c("Ca.(Intercept)","Ca.X","logCaw.(Intercept)","logCaw.X","logDaw.(Intercept)","logDaw.X")
  colnames(out) <- c("est","lb","ub","se","p")
  out <- as.data.frame(out)
  
  out[c("Ca.(Intercept)","Ca.X"),c("est","se","p")] <- summary(fit_Ca)$coef[,c("Estimate","Std. Error","Pr(>|t|)")]
  out[c("logCaw.(Intercept)","logCaw.X"),c("est","se","p")] <- summary(fit_logCaw)$coef[,c("Estimate","Std. Error","Pr(>|t|)")]
  out[c("logDaw.(Intercept)","logDaw.X"),c("est","se","p")] <- summary(fit_logDaw)$coef[,c("Estimate","Std. Error","Pr(>|t|)")]
  
  out[c("Ca.(Intercept)","Ca.X"),c("lb","ub")] <- confint(fit_Ca)[,c("2.5 %","97.5 %")]
  out[c("logCaw.(Intercept)","logCaw.X"),c("lb","ub")] <- confint(fit_logCaw)[,c("2.5 %","97.5 %")]
  out[c("logDaw.(Intercept)","logDaw.X"),c("lb","ub")] <- confint(fit_logDaw)[,c("2.5 %","97.5 %")]
  
  out
}


##### NLME ########

# functions for output NLME results
nlme_sep_FINAL <- function(myfit){
  # Extract EB estimates of NO parameters (State I)
  ranefdata         <-ranef(myfit,augFrame=T)
  ranefdata$Ca      <-ranefdata$Ca+fixef(myfit)[1]
  ranefdata$logCaw  <-ranefdata$logCaw+fixef(myfit)[2]
  ranefdata$logDaw  <-ranefdata$logDaw+fixef(myfit)[3]
  # Simple linear regressions (Stage II)
  fit_Ca     <- lm(Ca~X,     data=ranefdata)
  fit_logCaw <- lm(logCaw~X, data=ranefdata)
  fit_logDaw <- lm(logDaw~X, data=ranefdata)
  
  # output
  out <- matrix(NA,nrow=6,ncol=5)
  rownames(out) <- c("Ca.(Intercept)","Ca.X","logCaw.(Intercept)","logCaw.X","logDaw.(Intercept)","logDaw.X")
  colnames(out) <- c("est","lb","ub","se","p")
  out <- as.data.frame(out)
  
  out[c("Ca.(Intercept)","Ca.X"),c("est","se","p")] <- summary(fit_Ca)$coef[,c("Estimate","Std. Error","Pr(>|t|)")]
  out[c("logCaw.(Intercept)","logCaw.X"),c("est","se","p")] <- summary(fit_logCaw)$coef[,c("Estimate","Std. Error","Pr(>|t|)")]
  out[c("logDaw.(Intercept)","logDaw.X"),c("est","se","p")] <- summary(fit_logDaw)$coef[,c("Estimate","Std. Error","Pr(>|t|)")]
  
  out[c("Ca.(Intercept)","Ca.X"),c("lb","ub")] <- confint(fit_Ca)[,c("2.5 %","97.5 %")]
  out[c("logCaw.(Intercept)","logCaw.X"),c("lb","ub")] <- confint(fit_logCaw)[,c("2.5 %","97.5 %")]
  out[c("logDaw.(Intercept)","logDaw.X"),c("lb","ub")] <- confint(fit_logDaw)[,c("2.5 %","97.5 %")]
  
  out
}

nlme_sim_FINAL <- function(myfit){
  startFix <- fixef(myfit)
  # nlmeSim
  fit4.nlme <- update(myfit,
                      fixed = list(Ca ~ X, logCaw ~ X, logDaw ~ X),
                      start = c(startFix[1],0,startFix[2],0,startFix[3],0),
                      control=list(tolerance=0.01)
  )
  Ttable <- summary(fit4.nlme)$tTable
  lb_est_ub <- intervals(fit4.nlme)$fixed
  out <- data.frame(    est=Ttable[,"Value"],
                        lb=lb_est_ub[,"lower"],
                        ub=lb_est_ub[,"upper"],
                        se=Ttable[,"Std.Error"],
                        p=Ttable[,"p-value"]
  )
  out
  # FOR NLME_sim, you likely want to save additional output for more direct comparison with JAGS
  # like the VAR/COV of random effects, SD error
}



