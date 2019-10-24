load("/home/rcf-proj/spe1/wengjing/FeNO/simulation/Rdata/Summary/todolist.Rdata")
load("/Users/wengjingying/Downloads/FeNO/eno/Seedform.Rdata")


slurm_arrayid <-as.integer(Sys.getenv('SLURM_ARRAY_TASK_ID'))
repeatid<-unlist(faillist)[slurm_arrayid]
type<-as.character(stack(faillist)$ind[slurm_arrayid])

seed<-Seedmatrix[repeatid,type]
set.seed(seed)

setwd(paste("/home/rcf-proj/spe1/wengjing/FeNO/simulation/Rdata/freq-",type,sep=""))



load(paste("data_",repeatid,".Rdata",sep=""))
.libPaths("/home/rcf-proj/spe1/wengjing/R_packages")
source("/home/rcf-proj/spe1/wengjing/FeNO/simulation/Script/Func.R")
rhat=1.1
addon.iter=6000
Max_update=40
n.final=6000
Ndat=1000
N.iterT=12000
sigma2=1e-3

alphaC_prior=c(2,log(68),log(12))
betaC_prior =c(0,0,0)

flow <- c(rep(30,2),rep(50,2),rep(100,2),rep(300,2)) #specify flow rates

data_cov <- list(Ndat=length(enodata$X),Nflows=length(flow),flow=flow,logeno=enodata$logeno,X=enodata$X,
                 betaC_prior=betaC_prior,alphaC_prior=alphaC_prior)

inits_FcovC <- 
  list(
    list(
      tau_c=runif(1,0.05,400), # to have smaller SD init (0.05,5)
      alpha_Ca=runif(1,0,5),
      alpha_logCaw=runif(1,0,5),
      alpha_logDaw=runif(1,0,5),
      beta_Ca=runif(1,-1,1),
      beta_logCaw=runif(1,-1,1),
      beta_logDaw=runif(1,-1,1)
    ),
    list(
      tau_c=runif(1,0.05,400), # to have smaller SD init (0.05,5)
      alpha_Ca=runif(1,0,5),
      alpha_logCaw=runif(1,0,5),
      alpha_logDaw=runif(1,0,5),
      beta_Ca=runif(1,-1,1),
      beta_logCaw=runif(1,-1,1),
      beta_logDaw=runif(1,-1,1)
    ),
    list(
      tau_c=runif(1,0.05,400), # to have smaller SD init (0.05,5)
      alpha_Ca=runif(1,0,5),
      alpha_logCaw=runif(1,0,5),
      alpha_logDaw=runif(1,0,5),
      beta_Ca=runif(1,-1,1),
      beta_logCaw=runif(1,-1,1),
      beta_logDaw=runif(1,-1,1)
    )
  )


param.R2jags_NOTX <-c("alpha_Ca","alpha_logCaw","alpha_logDaw",
                      "beta_Ca","beta_logCaw","beta_logDaw",
                      "Ca","logCaw","logDaw",
                      "sdCa","sdlogCaw","sdlogDaw","corlogCawCa","corlogCawlogDaw","corlogDawCa",
                      "sigma_c")
time<-proc.time()[3]
#### truncated model
fit.R2jags_NOXT<- jags(data=data_cov, parameters.to.save=param.R2jags_NOTX, inits=inits_FcovC,
                       model.file=truncationNOX, n.chains=3, n.iter=N.iterT, n.burnin=N.iterT/3) 
N_updateT=0
N_final=0
while(N_updateT<Max_update){
  if(max(fit.R2jags_NOXT$BUGSoutput$summary[rownames(fit.R2jags_NOXT$BUGSoutput$summary) %in% param.R2jags_NOTX,"Rhat" ])<rhat){
    fit.R2jags_NOXT<-update(fit.R2jags_NOXT,n.iter=n.final)
    N_final=N_final+1
    if(max(fit.R2jags_NOXT$BUGSoutput$summary[rownames(fit.R2jags_NOXT$BUGSoutput$summary) %in% param.R2jags_NOTX,"Rhat" ])<rhat){
      break
    }
  }
  fit.R2jags_NOXT<-update(fit.R2jags_NOXT,n.iter=addon.iter)
  N_updateT<-N_updateT+1
}


timespend<-proc.time()[3]-time
trimmean<-apply(fit.R2jags_NOXT$BUGSoutput$sims.matrix,2,function(x) mean(x, 0.2))
quantile1_99<-t(apply(fit.R2jags_NOXT$BUGSoutput$sims.matrix,2,function(x) quantile(x, c(0.005,0.995))))
result<-list(
  timespend=timespend,
  summary=cbind(fit.R2jags_NOXT$BUGSoutput$summary,"trimmean"=trimmean,quantile1_99),
  Iter_T=N.iterT+N_updateT*addon.iter+n.final*N_final,#total not summary iter
  Iter_track=c(N.iterT=N.iterT,N_updateT=N_updateT,addon.iter=addon.iter,n.final=n.final,N_final=N_final),
  inits=inits_FcovC
  )


save(result, file=paste("/home/rcf-proj/spe1/wengjing/FeNO/simulation/Rdata/freq-",type,"/true_JagsT","result_",repeatid,".Rdata",sep=""))

print(fit.R2jags_NOXT$BUGSoutput$summary[rownames(fit.R2jags_NOXT$BUGSoutput$summary)%in% param.R2jags_NOTX,"Rhat" ])
      



