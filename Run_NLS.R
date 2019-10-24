rm(list=ls())
load("/home/rcf-proj/spe1/wengjing/FeNO/simulation/Rdata/Seedform.Rdata")
source("/home/rcf-proj/spe1/wengjing/FeNO/simulation/Script/NEWfreqFunc.R")
for(type in colnames(Seedmatrix)){
  starttime0=proc.time()[3]
  setwd(paste("/home/rcf-proj/spe1/wengjing/FeNO/simulation/Rdata/freq-",type,sep=""))
  slurm_arrayid <-as.integer(Sys.getenv('SLURM_ARRAY_TASK_ID'))
  
  filenameresult<-paste("NLS",type,"result_",sep="")
  filelist<-list.files()
  freqNEWresult<-filelist[grepl(filenameresult,filelist)]
  resultlist<-as.numeric(unlist(lapply(freqNEWresult,FUN=function(x) strsplit(strsplit(x,split="result_")[[1]][2],split="\\.")[[1]][1])))
  print(slurm_arrayid %in% resultlist)

  seed<-Seedmatrix[slurm_arrayid,type]
  set.seed(seed)
  
  load(paste("data_",slurm_arrayid,".Rdata",sep=""))
  Failrecord<-c("HMA"=0,"NLS"=0,"NLSw"=0,"nlmeStart"=0,"nlmeSep"=0,"nlmeSim"=0)
  flow <- c(  rep(30, 2),rep(50, 2),rep(100,2),rep(300,2))
  Ndat<-length(enodata$X)
  logenodata<-melt(data.frame(id=c(1:Ndat),enodata$logeno),id="id")
  logenodata<-logenodata[order(logenodata$id,decreasing = FALSE),]
  datOut <- data.frame(id=rep(c(1:Ndat),each=length(flow)),eno=NA,logeno=NA,flow=rep(flow,Ndat),X=rep(enodata$X,each=length(flow)))
  datOut$logeno <- logenodata$value
  datOut$eno    <-exp(datOut$logeno)
  datOut$invFlow<-1/datOut$flow
  datOut$enoflow<-datOut$flow*datOut$eno
  
  dat <- datOut
  X <- enodata$X
  
  starttime=proc.time()[3]
  
  fitnls <- nlsListStart(dat)
  nls_NA<-sum(is.na(coef(fitnls)[,1]))
  
  fitC.nls<-lapply(c(1:1000), FUN=function(sub){
    tryCatch(nlsConstrained(dat,sub),error = function(x) {
      print(paste(sub,"Error with nls_C",sep=":"));
      warnings();
      c(    Ca=NA,
                     logCaw=NA,
                     logDaw=NA
      )
    })
  })
  
  fitCm.nls<-data.frame(t(as.data.frame(fitC.nls,row.names = c("Ca","logCaw","logDaw"),col.names=c(1:1000))))
  nlsC_NA<-sum(is.na(fitCm.nls[,1]))
  
  NLSout <- tryCatch(nls_FINAL(fitnls,X), error = function(x) {
    print("Error with nls");
    warnings();
    data.frame(    est=rep(NA,6),
                   lb=rep(NA,6),
                   ub=rep(NA,6),
                   se=rep(NA,6),
                   p=rep(NA,6)
    )
  }
  )
  
  NLSout_w <- tryCatch(nlsW_FINAL(fitnls,X), error = function(x) {
    print("Error with nls_w");
    warnings();
    data.frame(    est=rep(NA,6),
                   lb=rep(NA,6),
                   ub=rep(NA,6),
                   se=rep(NA,6),
                   p=rep(NA,6)
    )
  }
  )
  
  NLSout_C <- tryCatch(nlsC_FINAL(fitCm.nls,X), error = function(x) {
    print("Error with nls");
    warnings();
    data.frame(    est=rep(NA,6),
                   lb=rep(NA,6),
                   ub=rep(NA,6),
                   se=rep(NA,6),
                   p=rep(NA,6)
    )
  }
  )
  
  elapsed=proc.time()[3]-starttime
  Freqout<-cbind(NLSout,NLSout_w,NLSout_C)
  colnames(Freqout)<-paste(rep(c("NLS","NLS_w","NLS_C"),each=5),rep(c("est","lb","ub","se","p"),3),sep="_")
  save(elapsed,Freqout,nls_NA,nlsC_NA, file=paste("NLS",type,"result_",slurm_arrayid,".Rdata",sep=""))

  print(type)
  proc.time()[3]-starttime0
}


