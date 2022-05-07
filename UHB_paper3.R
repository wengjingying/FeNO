rm(list=ls())
set.seed(2021)
library(R2jags)
load("/home/wengjing@PREVMED.USC.EDU/biostat2/spe/data/CHScohortE/eeNOdat_01jun16_noadj_ptype3_inclcov.RData")

Enodata<-dat[,c("id","yr","eno","flow","height","maleENO","age","race2n","aNO","bmipctl","asthcaseENO","alrgystatENO","townname","Analyzer","pflowlev")]

Enodata<-Enodata[complete.cases(Enodata),]
Enodata$logeno<-log(Enodata$eno)
Enodata$invFlow    <-1/Enodata$flow
Enodata$enoflow    <-Enodata$flow*Enodata$logeno
Enodata$yr_id<-paste(Enodata$yr,Enodata$id,sep="_")

#Combine Current allergy and not current allergy as one group
# 2 * 2 allergy * asthma
# 0: never allergy, no asthma as ref, 1: never allery, asthma, 
# 2: current or not current allergy, no asthma; 3: current or not current allergy, asthma

Enodata$AA2<-NA
Enodata$AA2[Enodata$asthcaseENO==1 & Enodata$alrgystatENO == "Current allergy"]="3:asthma&Allergy"
Enodata$AA2[Enodata$asthcaseENO==0 & Enodata$alrgystatENO == "Current allergy"]="2:allergyOnly"
Enodata$AA2[Enodata$asthcaseENO==1 & Enodata$alrgystatENO %in% c("Never allergy","Not current allergy")]="1:asthmaOnly"
Enodata$AA2[Enodata$asthcaseENO==0 & Enodata$alrgystatENO %in% c("Never allergy","Not current allergy")]="0:neither"
Enodata$AA2<-as.factor(Enodata$AA2)

Enodata$G1<-ifelse(Enodata$AA2=="0:neither",1,0)
Enodata$G2<-ifelse(Enodata$AA2=="1:asthmaOnly",1,0)
Enodata$G3<-ifelse(Enodata$AA2=="2:allergyOnly",1,0)
Enodata$G4<-ifelse(Enodata$AA2=="3:asthma&Allergy",1,0)

#center
Enodata$heightC<-Enodata$height-mean(Enodata$height,na.rm=T)
Enodata$ageC<-Enodata$age-mean(Enodata$age,na.rm=T)
Enodata$maleC<-Enodata$maleENO-0.5

save(Enodata,file="/home/wengjing@PREVMED.USC.EDU/CHS/CHS_paper3Data.Rdata")

Enodata8<-Enodata[Enodata$yr==8,]
Enodata10<-Enodata[Enodata$yr==10,]

UBtable8<-table(Enodata8$id)
offset8<-c(1,cumsum(UBtable8)+1)

UBtable10<-table(Enodata10$id)
offset10<-c(1,cumsum(UBtable10)+1)

#customize your Xs
ID_X2_y8<-aggregate(G2~id,data=Enodata8,FUN = mean)
ID_X3_y8<-aggregate(G3~id,data=Enodata8,FUN = mean)
ID_X4_y8<-aggregate(G4~id,data=Enodata8,FUN = mean)

ID_X2_y10<-aggregate(G2~id,data=Enodata10,FUN = mean)
ID_X3_y10<-aggregate(G3~id,data=Enodata10,FUN = mean)
ID_X4_y10<-aggregate(G4~id,data=Enodata10,FUN = mean)
##**************##



Jags_U_CHS <- function()  {
  for( i in 1 : Nid) {
    for( j  in offset[i]:(offset[i+1]-1)) {
      logeno[j] ~ dnorm(logmu[j],tau_c)
      mu[j] <- exp(logCaw[i])+(Ca[i]-exp(logCaw[i]))*exp(-1*exp(logDaw[i])/flow[j])
      logmu[j]<-log(max(mu[j] ,0.001))#logFeNO>0.001
    }
    
    Ca[i]       <- param_subject[i, 1]
    logCaw[i]   <- param_subject[i, 2] 
    logDaw[i]   <- param_subject[i, 3] 
    
    param_subject[i,1] ~ dnorm(mu1_22[i],tau1_22);T(0,)
    mu1_22[i]<-Ca_mean[i]+part1[1,1]*(param_subject[i,2]-logCaw_mean[i])+part1[1,2]*(param_subject[i,3]-logDaw_mean[i])
    param_subject[i,2:3] ~ dmnorm(c(logCaw_mean[i],logDaw_mean[i]),tau22[1:2,1:2])
    
    #customize your Xs
    Ca_mean[i]    <-beta0_Ca  +  beta1_Ca     * X1[i]+beta2_Ca     * X2[i]+beta3_Ca     * X3[i]
    logCaw_mean[i]<-beta0_logCaw+beta1_logCaw * X1[i]+beta2_logCaw * X2[i]+beta3_logCaw * X3[i]
    logDaw_mean[i]<-beta0_logDaw+beta1_logDaw * X1[i]+beta2_logDaw * X2[i]+beta3_logDaw * X3[i]
    ##**************##
  }
  
  
  
  #---sub---#
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
  ############# var-cov 
  covCalogCaw       <-corlogCawCa*sqrt(varCa*varlogCaw)
  covCalogDaw       <-corlogDawCa*sqrt(varCa*varlogDaw)
  covlogCawlogDaw   <-corlogCawlogDaw*sqrt(varlogCaw*varlogDaw)
  
  varCa<-1/tauCa
  varlogCaw<-1/taulogCaw
  varlogDaw<-1/taulogDaw
  sdCa     <- sqrt(varCa)
  sdlogCaw <- sqrt(varlogCaw)
  sdlogDaw <- sqrt(varlogDaw)
  sigma_c<-sqrt(1/tau_c)
  ############# sample correlation
  corlogDawCa   ~ dunif(L_corlogDawCa,U_corlogDawCa)
  L_corlogDawCa <- corlogCawCa*corlogCawlogDaw-sqrt((corlogCawlogDaw^2-1)*(corlogCawCa^2-1))
  U_corlogDawCa <- corlogCawCa*corlogCawlogDaw+sqrt((corlogCawlogDaw^2-1)*(corlogCawCa^2-1))
  
  corlogCawCa        ~ dunif(-1,1)
  corlogCawlogDaw    ~ dunif(-1,1)
  ############# sample tau
  tauCa ~ dgamma(1.0E-3, 1.0E-3)
  taulogCaw ~ dgamma(1.0E-3, 1.0E-3)
  taulogDaw ~ dgamma(1.0E-3, 1.0E-3)
  #---error---#
  tau_c       ~ dgamma(1.0E-3, 1.0E-3)#sigma_C
  ############# sample mean
  beta0_Ca     ~ dnorm(beta0_prior[1],1e-6)
  beta0_logCaw ~ dnorm(beta0_prior[2],1e-6)
  beta0_logDaw ~ dnorm(beta0_prior[3],1e-6)
  
  #customize your Xs
  beta1_Ca      ~ dnorm(betaC_prior[1], 1e-6)
  beta1_logCaw  ~ dnorm(betaC_prior[2], 1e-6)
  beta1_logDaw  ~ dnorm(betaC_prior[3], 1e-6)
  
  beta2_Ca      ~ dnorm(betaC_prior[1], 1e-6)
  beta2_logCaw  ~ dnorm(betaC_prior[2], 1e-6)
  beta2_logDaw  ~ dnorm(betaC_prior[3], 1e-6)
  
  beta3_Ca      ~ dnorm(betaC_prior[1], 1e-6)
  beta3_logCaw  ~ dnorm(betaC_prior[2], 1e-6)
  beta3_logDaw  ~ dnorm(betaC_prior[3], 1e-6)
  ##**************##
} 

rhat=1.1
NO_R_prior=diag(x=1e-3,nrow=3,ncol=3)
beta0_prior=c(2,log(68),log(12))
betaC_prior =c(0,0,0)
#flowarray=Arrayflow
#logenoarray=Arraylogeno
#X=ArrayheightS
#N_yr=N_yr
#idlist=idlist
#yr_perID=yr_perID
data_U_CHS_8 <- list(Nid=length(unique(Enodata8$id)),offset=offset8,
                      logeno=Enodata8$logeno,flow=Enodata8$flow,
                      #customize your Xs
                      X1=ID_X2_y8$G2,X2=ID_X3_y8$G3,X3=ID_X4_y8$G4,
                      ##**************##
                      beta0_prior=beta0_prior,betaC_prior=betaC_prior
)


param.U_CHS_8 <-c("beta0_Ca","beta0_logCaw","beta0_logDaw",
                   "Ca","logCaw","logDaw",
                   #costumize your Xs
                   "beta1_Ca","beta1_logCaw","beta1_logDaw",
                   "beta2_Ca","beta2_logCaw","beta2_logDaw",
                   "beta3_Ca","beta3_logCaw","beta3_logDaw",
                   ##**************##
                   "sdCa","sdlogCaw","sdlogDaw",
                   "corlogDawCa","corlogCawCa","corlogCawlogDaw",
                   "sigma_c")

param.Rhat.U        <-c(
  "beta0_Ca","beta0_logCaw","beta0_logDaw",
  #costumize your Xs
  "beta1_Ca","beta1_logCaw","beta1_logDaw",
  "beta2_Ca","beta2_logCaw","beta2_logDaw",
  "beta3_Ca","beta3_logCaw","beta3_logDaw",
  ##**************##
  "sdCa","sdlogCaw","sdlogDaw",
  "corlogDawCa","corlogCawCa","corlogCawlogDaw",
  "sigma_c")

SubM<-function(){
  sdCa            <- 1/sqrt(runif(1,1,20))
  sdlogCaw        <- 1/sqrt(runif(1,1,20))
  sdlogDaw        <- 1/sqrt(runif(1,1,20))
  
  corlogCawCa     <- runif(1,-1,1)
  corlogCawlogDaw <- runif(1,-1,1)
  L_corlogDawCa        <- corlogCawCa*corlogCawlogDaw-sqrt((corlogCawlogDaw^2-1)*(corlogCawCa^2-1))
  U_corlogDawCa        <- corlogCawCa*corlogCawlogDaw+sqrt((corlogCawlogDaw^2-1)*(corlogCawCa^2-1))
  corlogDawCa     <- runif(1,L_corlogDawCa,U_corlogDawCa)
  Sub_CovMat_chs_Caw<- matrix(c(sdCa^2,corlogCawCa*sdCa*sdlogCaw,corlogDawCa*sdCa*sdlogDaw,
                                corlogCawCa*sdCa*sdlogCaw,sdlogCaw^2,corlogCawlogDaw*sdlogCaw*sdlogDaw,
                                corlogDawCa*sdCa*sdlogDaw,corlogCawlogDaw*sdlogCaw*sdlogDaw,sdlogDaw^2),3,3)
  return(Sub_CovMat_chs_Caw)
}

inits_U <-
  list(
    list(
      tau_c               =runif(1,10,200),
      beta0_Ca            =runif(1,0,5),
      beta0_logCaw        =runif(1,0,5),
      beta0_logDaw        =runif(1,0,5),
      #customize your Xs
      beta1_Ca            =runif(1,-0.5,0.5),
      beta1_logCaw        =runif(1,-0.5,0.5),
      beta1_logDaw        =runif(1,-0.5,0.5),
      beta2_Ca            =runif(1,-0.5,0.5),
      beta2_logCaw        =runif(1,-0.5,0.5),
      beta2_logDaw        =runif(1,-0.5,0.5),
      beta3_Ca            =runif(1,-0.5,0.5),
      beta3_logCaw        =runif(1,-0.5,0.5),
      beta3_logDaw        =runif(1,-0.5,0.5),
      ##**************##
      tau_sub             =solve(SubM())
    ),
    list(
      tau_c               =runif(1,10,200),
      beta0_Ca            =runif(1,0,5),
      beta0_logCaw        =runif(1,0,5),
      beta0_logDaw        =runif(1,0,5),
      #customize your Xs
      beta1_Ca            =runif(1,-0.5,0.5),
      beta1_logCaw        =runif(1,-0.5,0.5),
      beta1_logDaw        =runif(1,-0.5,0.5),
      beta2_Ca            =runif(1,-0.5,0.5),
      beta2_logCaw        =runif(1,-0.5,0.5),
      beta2_logDaw        =runif(1,-0.5,0.5),
      beta3_Ca            =runif(1,-0.5,0.5),
      beta3_logCaw        =runif(1,-0.5,0.5),
      beta3_logDaw        =runif(1,-0.5,0.5),
      ##**************##
      tau_sub             =solve(SubM())
    ),
    list(
      tau_c               =runif(1,10,200),
      beta0_Ca            =runif(1,0,5),
      beta0_logCaw        =runif(1,0,5),
      beta0_logDaw        =runif(1,0,5),
      #customize your Xs
      beta1_Ca            =runif(1,-0.5,0.5),
      beta1_logCaw        =runif(1,-0.5,0.5),
      beta1_logDaw        =runif(1,-0.5,0.5),
      beta2_Ca            =runif(1,-0.5,0.5),
      beta2_logCaw        =runif(1,-0.5,0.5),
      beta2_logDaw        =runif(1,-0.5,0.5),
      beta3_Ca            =runif(1,-0.5,0.5),
      beta3_logCaw        =runif(1,-0.5,0.5),
      beta3_logDaw        =runif(1,-0.5,0.5),
      ##**************##
      tau_sub             =solve(SubM())
    )
  )



#### truncated model
N.iterT=12000
N.burn=8000
N.thinM=40
Max_update=20
n.final=6000
addon.iter=6000

time<-Sys.time()
LongiX_sim_2<- jags(data=data_U_CHS_8, parameters.to.save=param.U_CHS_8, inits=inits_U,
                    model.file=Jags_U_CHS, n.chains=3, n.iter=N.iterT, n.burnin=N.burn,n.thin=N.thinM) 
timespend<-Sys.time()-time
print(timespend)

N_updateT=0
N_final=0
while(N_updateT<Max_update){
  if(max(LongiX_sim_2$BUGSoutput$summary[rownames(LongiX_sim_2$BUGSoutput$summary) %in% param.Rhat.U,"Rhat" ])<rhat){
    LongiX_sim_2<-update(LongiX_sim_2,n.iter=n.final,n.thin=20)
    N_final=N_final+1
    print("final")
    print(N_final)
    print(LongiX_sim_2$BUGSoutput$summary[rownames(LongiX_sim_2$BUGSoutput$summary) %in% param.Rhat.U,"Rhat" ])
    if(max(LongiX_sim_2$BUGSoutput$summary[rownames(LongiX_sim_2$BUGSoutput$summary) %in% param.Rhat.U,"Rhat" ])<rhat){
      break
    }
  }
  LongiX_sim_2<-update(LongiX_sim_2,n.iter=addon.iter,n.thin=10)
  print(LongiX_sim_2$BUGSoutput$summary[rownames(LongiX_sim_2$BUGSoutput$summary) %in% param.Rhat.U,"Rhat" ])
  N_updateT<-N_updateT+1
  print("update")
  print(N_updateT)
  timespend<-Sys.time()-time
  print(timespend)
}
timespend<-Sys.time()-time
print(timespend)
trimmean<-apply(LongiX_sim_2$BUGSoutput$sims.matrix,2,function(x) mean(x, 0.2))
quantile1_99<-t(apply(LongiX_sim_2$BUGSoutput$sims.matrix,2,function(x) quantile(x, c(0.005,0.995))))
result<-list(
  timespend=timespend,
  summary=cbind(LongiX_sim_2$BUGSoutput$summary,"trimmean"=trimmean,quantile1_99),
  Iter_T=N.iterT+N_updateT*addon.iter+n.final*N_final,#total not summary iter
  Iter_track=c(N.iterT=N.iterT,N_updateT=N_updateT,addon.iter=addon.iter,n.final=n.final,N_final=N_final)
)



save(result,LongiX_sim_2,file="new_CHS_JAGS_U_paper3_y8_M1.Rdata")


save(result,file="new_D_CHS_JAGS_U_paper3_y8_M1.Rdata")
