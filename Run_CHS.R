setwd("/Users/wengjingying/Downloads/FeNO/eno/CHS")
source("prepCode_FirstATSCheckIn.R")
load("/Volumes/projects_eh/biostat2/spe/data/CHScohortE/eeNOdat2_diffPlatDefnsVolume28apr15noadj_firstATS.RData")

load("/Volumes/projects_eh/biostat2/spe/data/CHScohortE/Adat_diffPlatDefnsVolume28apr15noadj_firstATS_indoorNOpaper.RData")
fit <- lm(CaNO ~ I(aNO/10)#
          +agekidENO
          +male
          +Race2n
          +asthcase
          +asthmedCatMiss
          +AllergyStatus
          +hourcollCat
          #+bmipct
          +SHS
          +ParentEduc
          +Analyzer
          +ft_town
          ,
          dat=subset(Adat)
)
summary(fit)

Enodata<-dat[,c("id","eno","aNO","flow","agekidENO","male","Race2n",
                "asthcase","asthmedCat","AllergyStatus","hourcollCat",
                "SHS","ParentEduc","Analyzer","ft_town")]
Enodata$asthmedCatMiss<-ifelse(is.na(Enodata$asthmedCat),1,0)
Enodata$logeno<-log(Enodata$eno)
Enodata$aNO10<-Enodata$aNO/10
#use age asthcase SHS completed data
avelogaNO10<-mean(sqrt(Enodata$aNO10))
Enodata$Cl_aNO10<-sqrt(Enodata$aNO10)-avelogaNO10
# As in Adat: no NAs for age, 6 NA for race, no NA in asthcase, AlleryStatus,




CHS_Xs <- function()  {
  for( i in 1 : Ndat) {
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
    
    Ca_mean[i]<-alpha_Ca+beta_age_Ca * X[j,"agekidENO"]+beta_male_Ca * X[j,"male"]+beta_race_Ca * X[j,"Race2n"]
    +beta_asthcase_Ca * X[j,"asthcase"]+beta_asthMiss_Ca * X[j,"asthmedCatMiss"]+beta_Status_Ca * X[j,"AllergyStatus"]
    +beta_hour_Ca * X[j, "hourcollCat"]+beta_SHS_Ca * X[j,"SHS"]+beta_Pedu_Ca * X[j,"ParentEduc"]
    +beta_Anal_Ca * X[j,"Analyzer"]+beta_ft_Ca * X[j,"ft_town"]
    
    logCaw_mean[i]<-alpha_logCaw+beta_age_logCaw * X[j,"agekidENO"]+beta_male_logCaw * X[j,"male"]+beta_race_logCaw * X[j,"Race2n"]
    +beta_asthcase_logCaw * X[j,"asthcase"]+beta_asthMiss_logCaw * X[j,"asthmedCatMiss"]+beta_Status_logCaw * X[j,"AllergyStatus"]
    +beta_hour_logCaw * X[j, "hourcollCat"]+beta_SHS_logCaw * X[j,"SHS"]+beta_Pedu_logCaw * X[j,"ParentEduc"]
    +beta_Anal_logCaw * X[j,"Analyzer"]+beta_ft_logCaw * X[j,"ft_town"]
    
    logDaw_mean[i]<-alpha_logDaw+beta_age_logDaw * X[j,"agekidENO"]+beta_male_logDaw * X[j,"male"]+beta_race_logDaw * X[j,"Race2n"]
    +beta_asthcase_logDaw * X[j,"asthcase"]+beta_asthMiss_logDaw * X[j,"asthmedCatMiss"]+beta_Status_logDaw * X[j,"AllergyStatus"]
    +beta_hour_logDaw * X[j, "hourcollCat"]+beta_SHS_logDaw * X[j,"SHS"]+beta_Pedu_logDaw * X[j,"ParentEduc"]
    +beta_Anal_logDaw * X[j,"Analyzer"]+beta_ft_logDaw * X[j,"ft_town"]
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
  tau_c       ~ dgamma(1.0E-3, 1.0E-3)#random error in regression
  ############# sample mean
  alpha_Ca     ~ dnorm(alphaC_prior[1],1e-6)
  alpha_logCaw ~ dnorm(alphaC_prior[2],1e-6)
  alpha_logDaw ~ dnorm(alphaC_prior[3],1e-6)
  
  beta_age_Ca      ~ dnorm(betaC_prior[1], 1e-6)
  beta_age_logCaw  ~ dnorm(betaC_prior[2], 1e-6)
  beta_age_logDaw  ~ dnorm(betaC_prior[3], 1e-6)
  
  beta_male_Ca      ~ dnorm(betaC_prior[1], 1e-6)
  beta_male_logCaw  ~ dnorm(betaC_prior[2], 1e-6)
  beta_male_logDaw  ~ dnorm(betaC_prior[3], 1e-6)
  
  beta_raceA_Ca      ~ dnorm(betaC_prior[1], 1e-6)
  beta_raceA_logCaw  ~ dnorm(betaC_prior[2], 1e-6)
  beta_raceA_logDaw  ~ dnorm(betaC_prior[3], 1e-6)
  
  beta_raceW_Ca      ~ dnorm(betaC_prior[1], 1e-6)
  beta_raceW_logCaw  ~ dnorm(betaC_prior[2], 1e-6)
  beta_raceW_logDaw  ~ dnorm(betaC_prior[3], 1e-6)
  
  beta_raceH_Ca      ~ dnorm(betaC_prior[1], 1e-6)
  beta_raceH_logCaw  ~ dnorm(betaC_prior[2], 1e-6)
  beta_raceH_logDaw  ~ dnorm(betaC_prior[3], 1e-6)
  
  beta_raceB_Ca      ~ dnorm(betaC_prior[1], 1e-6)
  beta_raceB_logCaw  ~ dnorm(betaC_prior[2], 1e-6)
  beta_raceB_logDaw  ~ dnorm(betaC_prior[3], 1e-6)
  
  beta_raceO_Ca      ~ dnorm(betaC_prior[1], 1e-6)
  beta_raceO_logCaw  ~ dnorm(betaC_prior[2], 1e-6)
  beta_raceO_logDaw  ~ dnorm(betaC_prior[3], 1e-6)
  
  beta_raceM_Ca      ~ dnorm(betaC_prior[1], 1e-6)
  beta_raceM_logCaw  ~ dnorm(betaC_prior[2], 1e-6)
  beta_raceM_logDaw  ~ dnorm(betaC_prior[3], 1e-6)
  
  beta_asthcase_Ca      ~ dnorm(betaC_prior[1], 1e-6)
  beta_asthcase_logCaw  ~ dnorm(betaC_prior[2], 1e-6)
  beta_asthcase_logDaw  ~ dnorm(betaC_prior[3], 1e-6)
  
  beta_asthMiss_Ca      ~ dnorm(betaC_prior[1], 1e-6)
  beta_asthMiss_logCaw  ~ dnorm(betaC_prior[2], 1e-6)
  beta_asthMiss_logDaw  ~ dnorm(betaC_prior[3], 1e-6)
  
  beta_StatusNc_Ca      ~ dnorm(betaC_prior[1], 1e-6)
  beta_StatusNc_logCaw  ~ dnorm(betaC_prior[2], 1e-6)
  beta_StatusNc_logDaw  ~ dnorm(betaC_prior[3], 1e-6)
  
  beta_StatusC_Ca      ~ dnorm(betaC_prior[1], 1e-6)
  beta_StatusC_logCaw  ~ dnorm(betaC_prior[2], 1e-6)
  beta_StatusC_logDaw  ~ dnorm(betaC_prior[3], 1e-6)
  
  beta_StatusNe_Ca      ~ dnorm(betaC_prior[1], 1e-6)
  beta_StatusNe_logCaw  ~ dnorm(betaC_prior[2], 1e-6)
  beta_StatusNe_logDaw  ~ dnorm(betaC_prior[3], 1e-6)
  
  beta_StatusM_Ca      ~ dnorm(betaC_prior[1], 1e-6)
  beta_StatusM_logCaw  ~ dnorm(betaC_prior[2], 1e-6)
  beta_StatusM_logDaw  ~ dnorm(betaC_prior[3], 1e-6)
  
  beta_hour1_Ca      ~ dnorm(betaC_prior[1], 1e-6)
  beta_hour1_logCaw  ~ dnorm(betaC_prior[2], 1e-6)
  beta_hour1_logDaw  ~ dnorm(betaC_prior[3], 1e-6)
  
  beta_hour2_Ca      ~ dnorm(betaC_prior[1], 1e-6)
  beta_hour2_logCaw  ~ dnorm(betaC_prior[2], 1e-6)
  beta_hour2_logDaw  ~ dnorm(betaC_prior[3], 1e-6)
  
  beta_hour3_Ca      ~ dnorm(betaC_prior[1], 1e-6)
  beta_hour3_logCaw  ~ dnorm(betaC_prior[2], 1e-6)
  beta_hour3_logDaw  ~ dnorm(betaC_prior[3], 1e-6)
  
  beta_SHS_Ca      ~ dnorm(betaC_prior[1], 1e-6)
  beta_SHS_logCaw  ~ dnorm(betaC_prior[2], 1e-6)
  beta_SHS_logDaw  ~ dnorm(betaC_prior[3], 1e-6)
  
  beta_Pedu1_Ca      ~ dnorm(betaC_prior[1], 1e-6)
  beta_Pedu1_logCaw  ~ dnorm(betaC_prior[2], 1e-6)
  beta_Pedu1_logDaw  ~ dnorm(betaC_prior[3], 1e-6)
  
  beta_Pedu2_Ca      ~ dnorm(betaC_prior[1], 1e-6)
  beta_Pedu2_logCaw  ~ dnorm(betaC_prior[2], 1e-6)
  beta_Pedu2_logDaw  ~ dnorm(betaC_prior[3], 1e-6)
  
  beta_Pedu3_Ca      ~ dnorm(betaC_prior[1], 1e-6)
  beta_Pedu3_logCaw  ~ dnorm(betaC_prior[2], 1e-6)
  beta_Pedu3_logDaw  ~ dnorm(betaC_prior[3], 1e-6)
  
  beta_PeduM_Ca      ~ dnorm(betaC_prior[1], 1e-6)
  beta_PeduM_logCaw  ~ dnorm(betaC_prior[2], 1e-6)
  beta_PeduM_logDaw  ~ dnorm(betaC_prior[3], 1e-6)
  
  beta_AnalA_Ca      ~ dnorm(betaC_prior[1], 1e-6)
  beta_AnalA_logCaw  ~ dnorm(betaC_prior[2], 1e-6)
  beta_AnalA_logDaw  ~ dnorm(betaC_prior[3], 1e-6)
  
  beta_AnalK_Ca      ~ dnorm(betaC_prior[1], 1e-6)
  beta_AnalK_logCaw  ~ dnorm(betaC_prior[2], 1e-6)
  beta_AnalK_logDaw  ~ dnorm(betaC_prior[3], 1e-6)
  
  beta_AnalN_Ca      ~ dnorm(betaC_prior[1], 1e-6)
  beta_AnalN_logCaw  ~ dnorm(betaC_prior[2], 1e-6)
  beta_AnalN_logDaw  ~ dnorm(betaC_prior[3], 1e-6)
  
  beta_ft_Ca      ~ dnorm(betaC_prior[1], 1e-6)
  beta_ft_logCaw  ~ dnorm(betaC_prior[2], 1e-6)
  beta_ft_logDaw  ~ dnorm(betaC_prior[3], 1e-6)
  
  beta_ft_Ca      ~ dnorm(betaC_prior[1], 1e-6)
  beta_ft_logCaw  ~ dnorm(betaC_prior[2], 1e-6)
  beta_ft_logDaw  ~ dnorm(betaC_prior[3], 1e-6)
  
  beta_ft_Ca      ~ dnorm(betaC_prior[1], 1e-6)
  beta_ft_logCaw  ~ dnorm(betaC_prior[2], 1e-6)
  beta_ft_logDaw  ~ dnorm(betaC_prior[3], 1e-6)
  
  beta_ft_Ca      ~ dnorm(betaC_prior[1], 1e-6)
  beta_ft_logCaw  ~ dnorm(betaC_prior[2], 1e-6)
  beta_ft_logDaw  ~ dnorm(betaC_prior[3], 1e-6)
  
  beta_ft_Ca      ~ dnorm(betaC_prior[1], 1e-6)
  beta_ft_logCaw  ~ dnorm(betaC_prior[2], 1e-6)
  beta_ft_logDaw  ~ dnorm(betaC_prior[3], 1e-6)
  
  beta_ft_Ca      ~ dnorm(betaC_prior[1], 1e-6)
  beta_ft_logCaw  ~ dnorm(betaC_prior[2], 1e-6)
  beta_ft_logDaw  ~ dnorm(betaC_prior[3], 1e-6)
} 


param.CHS_Xs <-c("alpha_Ca","alpha_logCaw","alpha_logDaw",
             "beta_age_Ca","beta_male_Ca","beta_race_Ca","beta_asthcase_Ca","beta_asthMiss_Ca",
             "beta_Status_Ca","beta_hour_Ca","beta_SHS_Ca","beta_Pedu_Ca","beta_Anal_Ca","beta_ft_Ca",
             "beta_age_logCaw","beta_male_logCaw","beta_race_logCaw","beta_asthcase_logCaw","beta_asthMiss_logCaw",
             "beta_Status_logCaw","beta_hour_logCaw","beta_SHS_logCaw","beta_Pedu_logCaw","beta_Anal_logCaw","beta_ft_logCaw",
             "beta_age_logDaw","beta_male_logDaw","beta_race_logDaw","beta_asthcase_logDaw","beta_asthMiss_logDaw",
             "beta_Status_logDaw","beta_hour_logDaw","beta_SHS_logDaw","beta_Pedu_logDaw","beta_Anal_logDaw","beta_ft_logDaw",
             "Ca","logCaw","logDaw",
             "sdCa","sdlogCaw","sdlogDaw","corlogCawCa","corlogCawlogDaw","corlogDawCa",
             "sigma_c")

inits_CHS <- 
  list(
    list(
      tau_c=runif(1,0.01,20), #could be unif too
      alpha_Ca=runif(1,0,5),
      alpha_logCaw=runif(1,0,5),
      alpha_logDaw=runif(1,0,5),
      beta_age_Ca         =runif(1,-1,1),
      beta_age_logCaw     =runif(1,-1,1),
      beta_age_logDaw     =runif(1,-1,1),
      beta_male_Ca        =runif(1,-1,1),
      beta_male_logCaw    =runif(1,-1,1),
      beta_male_logDaw    =runif(1,-1,1),
      beta_race_Ca        =runif(1,-1,1),
      beta_race_logCaw    =runif(1,-1,1),
      beta_race_logDaw    =runif(1,-1,1),
      beta_asthcase_Ca    =runif(1,-1,1),
      beta_asthcase_logCaw=runif(1,-1,1),
      beta_asthcase_logDaw=runif(1,-1,1),
      beta_asthMiss_Ca    =runif(1,-1,1),
      beta_asthMiss_logCaw=runif(1,-1,1),
      beta_asthMiss_logDaw=runif(1,-1,1),
      beta_Status_Ca      =runif(1,-1,1),
      beta_Status_logCaw  =runif(1,-1,1),
      beta_Status_logDaw  =runif(1,-1,1),
      beta_hour_Ca        =runif(1,-1,1),
      beta_hour_logCaw    =runif(1,-1,1),
      beta_hour_logDaw    =runif(1,-1,1),
      beta_SHS_Ca         =runif(1,-1,1),
      beta_SHS_logCaw     =runif(1,-1,1),
      beta_SHS_logDaw     =runif(1,-1,1),
      beta_Pedu_Ca        =runif(1,-1,1),
      beta_Pedu_logCaw    =runif(1,-1,1),
      beta_Pedu_logDaw    =runif(1,-1,1),
      beta_Anal_Ca        =runif(1,-1,1),
      beta_Anal_logCaw    =runif(1,-1,1),
      beta_Anal_logDaw    =runif(1,-1,1),
      beta_ft_Ca          =runif(1,-1,1),
      beta_ft_logCaw      =runif(1,-1,1),
      beta_ft_logDaw      =runif(1,-1,1)
    ),
    list(
      tau_c=runif(1,0.01,20),
      alpha_Ca=runif(1,0,5),
      alpha_logCaw=runif(1,0,5),
      alpha_logDaw=runif(1,0,5),
      beta_age_Ca         =runif(1,-1,1),
      beta_age_logCaw     =runif(1,-1,1),
      beta_age_logDaw     =runif(1,-1,1),
      beta_male_Ca        =runif(1,-1,1),
      beta_male_logCaw    =runif(1,-1,1),
      beta_male_logDaw    =runif(1,-1,1),
      beta_race_Ca        =runif(1,-1,1),
      beta_race_logCaw    =runif(1,-1,1),
      beta_race_logDaw    =runif(1,-1,1),
      beta_asthcase_Ca    =runif(1,-1,1),
      beta_asthcase_logCaw=runif(1,-1,1),
      beta_asthcase_logDaw=runif(1,-1,1),
      beta_asthMiss_Ca    =runif(1,-1,1),
      beta_asthMiss_logCaw=runif(1,-1,1),
      beta_asthMiss_logDaw=runif(1,-1,1),
      beta_Status_Ca      =runif(1,-1,1),
      beta_Status_logCaw  =runif(1,-1,1),
      beta_Status_logDaw  =runif(1,-1,1),
      beta_hour_Ca        =runif(1,-1,1),
      beta_hour_logCaw    =runif(1,-1,1),
      beta_hour_logDaw    =runif(1,-1,1),
      beta_SHS_Ca         =runif(1,-1,1),
      beta_SHS_logCaw     =runif(1,-1,1),
      beta_SHS_logDaw     =runif(1,-1,1),
      beta_Pedu_Ca        =runif(1,-1,1),
      beta_Pedu_logCaw    =runif(1,-1,1),
      beta_Pedu_logDaw    =runif(1,-1,1),
      beta_Anal_Ca        =runif(1,-1,1),
      beta_Anal_logCaw    =runif(1,-1,1),
      beta_Anal_logDaw    =runif(1,-1,1),
      beta_ft_Ca          =runif(1,-1,1),
      beta_ft_logCaw      =runif(1,-1,1),
      beta_ft_logDaw      =runif(1,-1,1)
    ),
    list(
      tau_c=runif(1,0.01,20),
      alpha_Ca=runif(1,0,5),
      alpha_logCaw=runif(1,0,5),
      alpha_logDaw=runif(1,0,5),
      beta_age_Ca         =runif(1,-1,1),
      beta_age_logCaw     =runif(1,-1,1),
      beta_age_logDaw     =runif(1,-1,1),
      beta_male_Ca        =runif(1,-1,1),
      beta_male_logCaw    =runif(1,-1,1),
      beta_male_logDaw    =runif(1,-1,1),
      beta_race_Ca        =runif(1,-1,1),
      beta_race_logCaw    =runif(1,-1,1),
      beta_race_logDaw    =runif(1,-1,1),
      beta_asthcase_Ca    =runif(1,-1,1),
      beta_asthcase_logCaw=runif(1,-1,1),
      beta_asthcase_logDaw=runif(1,-1,1),
      beta_asthMiss_Ca    =runif(1,-1,1),
      beta_asthMiss_logCaw=runif(1,-1,1),
      beta_asthMiss_logDaw=runif(1,-1,1),
      beta_Status_Ca      =runif(1,-1,1),
      beta_Status_logCaw  =runif(1,-1,1),
      beta_Status_logDaw  =runif(1,-1,1),
      beta_hour_Ca        =runif(1,-1,1),
      beta_hour_logCaw    =runif(1,-1,1),
      beta_hour_logDaw    =runif(1,-1,1),
      beta_SHS_Ca         =runif(1,-1,1),
      beta_SHS_logCaw     =runif(1,-1,1),
      beta_SHS_logDaw     =runif(1,-1,1),
      beta_Pedu_Ca        =runif(1,-1,1),
      beta_Pedu_logCaw    =runif(1,-1,1),
      beta_Pedu_logDaw    =runif(1,-1,1),
      beta_Anal_Ca        =runif(1,-1,1),
      beta_Anal_logCaw    =runif(1,-1,1),
      beta_Anal_logDaw    =runif(1,-1,1),
      beta_ft_Ca          =runif(1,-1,1),
      beta_ft_logCaw      =runif(1,-1,1),
      beta_ft_logDaw      =runif(1,-1,1)
    )
  )


rhat=1.1
addon.iter=100
Max_update=2
n.final=100
N.iterT=100
sigma2=1e-3

UBtable<-table(Enodata$id)
offset<-c(1,cumsum(UBtable)+1)
Ndat<-length(unique(Enodata$id))
flow<-Enodata$flow
alphaC_prior=c(2,log(68),log(12))
betaC_prior =c(0,0,0)
data_cov <- list(Ndat=Ndat,offset=offset,flow=flow,logeno=Enodata$logeno,
                 X=Enodata[,c("agekidENO","male","Race2n","asthcase","asthmedCat","AllergyStatus","hourcollCat","SHS","ParentEduc","Analyzer","ft_town")],
                 betaC_prior=betaC_prior,alphaC_prior=alphaC_prior)


time<-proc.time()[3]
#### truncated model

fit.CHS<- jags(data=data_cov, parameters.to.save=param.CHS_Xs, inits=inits_CHS,
                       model.file=CHS_Xs, n.chains=3, n.iter=N.iterT, n.burnin=N.iterT/3) 
N_updateT=0
N_final=0
while(N_updateT<Max_update){
  if(max(fit.CHS$BUGSoutput$summary[rownames(fit.CHS$BUGSoutput$summary) %in% param.CHS_Xs,"Rhat" ])<rhat){
    fit.CHS<-update(fit.CHS,n.iter=n.final)
    N_final=N_final+1
    print("final")
    print(N_final)
    prit(fit.CHS$BUGSoutput$summary[rownames(fit.CHS$BUGSoutput$summary) %in% param.CHS_Xs,"Rhat" ])
    if(max(fit.CHS$BUGSoutput$summary[rownames(fit.CHS$BUGSoutput$summary) %in% param.CHS_Xs,"Rhat" ])<rhat){
      break
    }
  }
  fit.CHS<-update(fit.CHS,n.iter=addon.iter)
  N_updateT<-N_updateT+1
  print("update")
  print(N_updateT)
}

trimmean<-apply(fit.CHS$BUGSoutput$sims.matrix,2,function(x) mean(x, 0.2))
quantile1_99<-t(apply(fit.CHS$BUGSoutput$sims.matrix,2,function(x) quantile(x, c(0.005,0.995))))
result<-list(
  timespend=timespend,
  summary=cbind(fit.CHS$BUGSoutput$summary,"trimmean"=trimmean,quantile1_99),
  Iter_T=N.iterT+N_updateT*addon.iter+n.final*N_final,#total not summary iter
  Iter_track=c(N.iterT=N.iterT,N_updateT=N_updateT,addon.iter=addon.iter,n.final=n.final,N_final=N_final),
)
timespend<-proc.time()[3]-time
print(timespend)

save(result, file="/home/rcf-proj/spe1/wengjing/FeNO/simulation/Rdata/CHS/MultiX.Rdata")

