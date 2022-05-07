library(MASS)
library(lme4)
library(nlme)
library(reshape2)
#library(gplots)
library(plotrix)
library(latex2exp)

jpeg("K:/paper/paper3/P3_plot1.jpeg", width = 900,height = 720)
par(mfrow=c(1,1),mgp=c(1,0.5,0),oma=c(3.5,1,2,0),mar=c(3,6,1,1),xaxs="i",yaxs="i")
xvals=c(1,4,7,10,14,17,20,23,27,30,33,36)

#NLS M1/M2
plotCI(x=xvals,err="y",
       y=c(0,NLSY8$Ca_est[2:4],0,NLSY8$logCaw_est[2:4],0,NLSY8$logDaw_est[2:4])+10.5,
       li=c(0,NLSY8$Ca_lb[2:4],0,NLSY8$logCaw_lb[2:4],0,NLSY8$logDaw_lb[2:4])+10.5,
       ui=c(0,NLSY8$Ca_ub[2:4],0,NLSY8$logCaw_ub[2:4],0,NLSY8$logDaw_ub[2:4])+10.5,
       yaxt="n",
       xaxt="n",xlim=c(0,38.5),ylim=c(0,12),xlab="",ylab="",pch=19,lty=1,cex=1.5)
plotCI(x=xvals+1,
       y=c(0,NLSY8_M2$Ca_est[2:4],0,NLSY8_M2$logCaw_est[2:4],0,NLSY8_M2$logDaw_est[2:4])+10.5,
       li=c(0,NLSY8_M2$Ca_lb[2:4],0,NLSY8_M2$logCaw_lb[2:4],0,NLSY8_M2$logDaw_lb[2:4])+10.5,
       ui=c(0,NLSY8_M2$Ca_ub[2:4],0,NLSY8_M2$logCaw_ub[2:4],0,NLSY8_M2$logDaw_ub[2:4])+10.5,
       pch=15,lty=2,add=TRUE,cex=1.5)

#HMA M1/M2
plotCI(x=xvals,err="y",
       y=c(0,HMAY8$Ca_est[2:4],0,HMAY8$logCaw_est[2:4],0,HMAY8$logDaw_est[2:4])+8.5,
       li=c(0,HMAY8$Ca_lb[2:4],0,HMAY8$logCaw_lb[2:4],0,HMAY8$logDaw_lb[2:4])+8.5,
       ui=c(0,HMAY8$Ca_ub[2:4],0,HMAY8$logCaw_ub[2:4],0,HMAY8$logDaw_ub[2:4])+8.5,pch=19,lty=1,add=TRUE,cex=1.5)
plotCI(x=xvals+1,
       y=c(0,HMAY8_M2$Ca_est[2:4],0,HMAY8_M2$logCaw_est[2:4],0,HMAY8_M2$logDaw_est[2:4])+8.5,
       li=c(0,HMAY8_M2$Ca_lb[2:4],0,HMAY8_M2$logCaw_lb[2:4],0,HMAY8_M2$logDaw_lb[2:4])+8.5,
       ui=c(0,HMAY8_M2$Ca_ub[2:4],0,HMAY8_M2$logCaw_ub[2:4],0,HMAY8_M2$logDaw_ub[2:4])+8.5,
       pch=15,slty=2,add=TRUE,cex=1.5)

#TSNLME M1/M2
plotCI(x=xvals,err="y",
       y=c(0,TSNLMEY8_M1$Ca_est[2:4],0,TSNLMEY8_M1$logCaw_est[2:4],0,TSNLMEY8_M1$logDaw_est[2:4])+6.5,
       li=c(0,TSNLMEY8_M1$Ca_lb[2:4],0,TSNLMEY8_M1$logCaw_lb[2:4],0,TSNLMEY8_M1$logDaw_lb[2:4])+6.5,
       ui=c(0,TSNLMEY8_M1$Ca_ub[2:4],0,TSNLMEY8_M1$logCaw_ub[2:4],0,TSNLMEY8_M1$logDaw_ub[2:4])+6.5,pch=19,lty=1,add=TRUE,cex=1.5)

plotCI(x=xvals+1,
       y=c(0,TSNLMEY8_M2$Ca_est[2:4],0,TSNLMEY8_M2$logCaw_est[2:4],0,TSNLMEY8_M2$logDaw_est[2:4])+6.5,
       li=c(0,TSNLMEY8_M2$Ca_lb[2:4],0,TSNLMEY8_M2$logCaw_lb[2:4],0,TSNLMEY8_M2$logDaw_lb[2:4])+6.5,
       ui=c(0,TSNLMEY8_M2$Ca_ub[2:4],0,TSNLMEY8_M2$logCaw_ub[2:4],0,TSNLMEY8_M2$logDaw_ub[2:4])+6.5,
       pch=15,slty=2,add=TRUE,cex=1.5)

#UNLME M1/M2
plotCI(x=xvals,err="y",
       y=c(0,UNLME_Y8$Ca_est[2:4],0,UNLME_Y8$logCaw_est[2:4],0,UNLME_Y8$logDaw_est[2:4])+4.5,
       li=c(0,UNLME_Y8$Ca_lb[2:4],0,UNLME_Y8$logCaw_lb[2:4],0,UNLME_Y8$logDaw_lb[2:4])+4.5,
       ui=c(0,UNLME_Y8$Ca_ub[2:4],0,UNLME_Y8$logCaw_ub[2:4],0,UNLME_Y8$logDaw_ub[2:4])+4.5,pch=19,lty=1,add=TRUE,cex=1.5)
plotCI(x=xvals+1,
       y=c(0,UNLME_Y8_M2$Ca_est[2:4],0,UNLME_Y8_M2$logCaw_est[2:4],0,UNLME_Y8_M2$logDaw_est[2:4])+4.5,
       li=c(0,UNLME_Y8_M2$Ca_lb[2:4],0,UNLME_Y8_M2$logCaw_lb[2:4],0,UNLME_Y8_M2$logDaw_lb[2:4])+4.5,
       ui=c(0,UNLME_Y8_M2$Ca_ub[2:4],0,UNLME_Y8_M2$logCaw_ub[2:4],0,UNLME_Y8_M2$logDaw_ub[2:4])+4.5,
       pch=15,slty=2,add=TRUE,cex=1.5)

#TSHB 
plotCI(x=xvals,err="y",
       y=c(0,TSHBY8_M1$Ca_est[2:4],0,TSHBY8_M1$logCaw_est[2:4],0,TSHBY8_M1$logDaw_est[2:4])+2.5,
       li=c(0,TSHBY8_M1$Ca_lb[2:4],0,TSHBY8_M1$logCaw_lb[2:4],0,TSHBY8_M1$logDaw_lb[2:4])+2.5,
       ui=c(0,TSHBY8_M1$Ca_ub[2:4],0,TSHBY8_M1$logCaw_ub[2:4],0,TSHBY8_M1$logDaw_ub[2:4])+2.5,pch=19,lty=1,add=TRUE,cex=1.5)
plotCI(x=xvals+1,
       y=c(0,TSHBY8_M2$Ca_est[2:4],0,TSHBY8_M2$logCaw_est[2:4],0,TSHBY8_M2$logDaw_est[2:4])+2.5,
       li=c(0,TSHBY8_M2$Ca_lb[2:4],0,TSHBY8_M2$logCaw_lb[2:4],0,TSHBY8_M2$logDaw_lb[2:4])+2.5,
       ui=c(0,TSHBY8_M2$Ca_ub[2:4],0,TSHBY8_M2$logCaw_ub[2:4],0,TSHBY8_M2$logDaw_ub[2:4])+2.5,
       pch=15,slty=2,add=TRUE,cex=1.5)
#UHB
plotCI(x=xvals,err="y",
       y=c(0,UHBY8_M1$Ca_est[2:4],0,UHBY8_M1$logCaw_est[2:4],0,UHBY8_M1$logDaw_est[2:4])+0.5,
       li=c(0,UHBY8_M1$Ca_lb[2:4],0,UHBY8_M1$logCaw_lb[2:4],0,UHBY8_M1$logDaw_lb[2:4])+0.5,
       ui=c(0,UHBY8_M1$Ca_ub[2:4],0,UHBY8_M1$logCaw_ub[2:4],0,UHBY8_M1$logDaw_ub[2:4])+0.5,pch=19,lty=1,add=TRUE,cex=1.5)



abline(v=c(12.5,25.5))
abline(h=c(2,4,6,8,10))
abline(h=c(0,2,4,6,8,10)+0.5,lty=2)
#title(ylab="Estimation",cex.lab=1,line=1.5)
axis(3,at=c(6,19,32),cex.lab=1.2,labels=c(TeX('$\\beta_{C_{A}}$'),TeX('$\\beta_{logCaw}$'),TeX('$\\beta_{logDaw}$')),tick=F,line=-0.5,cex.axis=1.2)
#axis(2,at=c(0,2,4,6,8,10)+0.5,cex.lab=1.2,labels=rep(0,6),las=2)
axis(1,at=xvals,labels=rep(c("None","Asthma","Allergy","Both"),3),cex.axis=0.8)
axis(2,at=c(1,3,5,7,9,11),labels=c("L_U_HB","L_TS_HB","L_U_NLME","L_TS_NLME","L_TS_HMA","L_TS_NLS"),las=1,cex.axis=1.2)

par(fig=c(0,1,0,1),oma=c(0,0,0,0),mar=c(0,0,0,0),new=TRUE)
plot(0,0,type="l",bty="n",xaxt="n",yaxt="n")
legend("bottom",
       legend=c("Unadjusted","Adjusted"),
       pch=c(19,15),
       cex=1.5,bty="n",xpd=TRUE,horiz=TRUE)
dev.off()