## R program for HIV data analysis in the <<Statistica Sinica>> paper titled
## "Sieve estimation of a class of partially linear transformation models 
## with interval-censored competing risks data";
## Author: Xuewen Lu; Edited by Dipankar Bandyopadhyay
## Date: July 16, 2021


#After running "Gridsearch.R", or this program with nboot=0 (it runs faster),
#with grid search, we obtained optimal alpha=(2,0), with maximum log-likelihood value. 
#Then, we run bootstrap with nboot=500 to get more accurate SE estimates;
#See the R output results at the end;


rm(list=ls())
sink("BIC-optimal-boot.txt", append=T)
cat("The program started:", date(), "\n")
options(width=60,length=200)
set.seed(20835)
library(basefun)  
library(foreach)
library(doParallel)
library(intccr)


#Function to generate (m+1) Bernstein polynomial basis functions;

m=2
nbt=500 #number of bootstrap samples;
Lambda.basis<- function(w,m){ 	
  minTX=min(w)
  maxTX=max(w)
  #bb.TX <- Bernstein_basis(numeric_var("x", support = c(minTX, maxTX)), order = 4, ui = "increasing")       #For m=4, 5 basis functions;
  bb.TX <- Bernstein_basis(numeric_var("x", support = c(minTX, maxTX)), order = m, ui = "increasing")       #For m=2, then 3 basis functions;
  rep(minTX,length(w))->MIN
  Lam.basis<-as.data.frame(bb.TX(w)-bb.TX(MIN))
  return(list(Lam.basis))
}



# Call your dataset here 
repro<-read.csv("iedea.csv",header=TRUE)

repro[,5] -> age
repro[,6] -> male
repro[,7] -> art_cd4

ns=length(age)
Age.mat<-matrix(0, ns, (m+1))
Bernstein.polynomial <- Lambda.basis(age,m)[[1]]	

#for (j in (1:(m+1))){
#	Age.mat[,j]=Bernstein.polynomial[,j]
#}
#Age=as.data.frame(Age.mat)	
	
Age1 <- Bernstein.polynomial[,1]
Age2 <- Bernstein.polynomial[,2]
Age3 <- Bernstein.polynomial[,3]
#Age4 <- Bernstein.polynomial[,4]
#Age5 <- Bernstein.polynomial[,5]

#cbind(repro,Age1,Age2,Age3,Age4,Age5) -> repro

cbind(repro,Age1,Age2, Age3) -> repro 
#need this to use Bootstrap for resampling;

#alpha1 <- seq(12.5,50,0.5)
#alpha2 <- seq(12.5,50,0.5)

#alpha1 <- seq(0,1,0.5)
#alpha2 <- seq(0,1,0.5)

#alpha1 <- seq(0,4,0.5)
#alpha2 <- seq(0,4,0.5)

##only run at optimal alpha=(2,0) selected from grid search;
alpha1 <- c(2)
alpha2 <- c(0)

alpha.init <- merge(alpha1,alpha2)
llh <- rep(0,dim(alpha.init)[1])
nsim=dim(alpha.init)[1]

#cl <- makeCluster(5)
#registerDoParallel(cl)
#Q <- foreach(i=1:nsim,.combine='rbind',.packages=c("intccr","basefun","doParallel")) %dopar%
#  {
   # fit <- ciregicBerns(formula=Surv2(v, u, c)~male+art_cd4+Age1+Age2+Age3+Age4+Age5, #data=repro, alpha=c(alpha[i,1],alpha[i,2]), do.par= TRUE, nboot=0)
#   fit <- ciregic(formula=Surv2(v, u, c)~male+art_cd4+Age1+Age2+Age3+Age4+Age5, data=repro, alpha=c(alpha[i,1],alpha[i,2]), do.par= TRUE, nboot=0)
#    cat("This is the", i, "th iteration, loglik=", fit$loglikelihood,  "\n") 
#    c(i,fit$loglikelihood)
#  }
  
#stopCluster(cl)

Q=matrix(0, nsim,2)
BIC=Q
#est.coef=matrix(0,nsim,14)      #14=7x2;
#Since only Age1 and Age4 are significant, we remove Age2, Age3 and Age5 and refit;
#Use m=2, and get Age1, Age2 and Age3, only;

est.coef=matrix(0,nsim,(2+(m+1))*2)      # 2= # of beta parameters, m+1=# of basis functions;

Alist=vector("list", nsim)

for (i in 1:nsim){
#fit <- ciregic(formula=Surv2(v, u, event=c)~male+art_cd4+Age1+Age2+Age3+Age4+Age5, data=repro, #alpha=c(alpha.init[i,1],alpha.init[i,2]), do.par= TRUE, nboot=0)

#fit <- ciregic(formula=Surv2(v, u, event=c)~male+art_cd4+Age1+Age2+Age3, data=repro, alpha=c(alpha.init[i,1],alpha.init[i,2]), k=1, do.par= TRUE, nboot=0)

fit <- ciregic(formula=Surv2(v, u, event=c)~male+art_cd4+Age1+Age2+Age3, data=repro, alpha=c(alpha.init[i,1],alpha.init[i,2]), do.par= TRUE, nboot=nbt) #default k=1;

#fit <- ciregic(formula=Surv2(v, u, c)~male+art_cd4+Age1, data=repro, alpha=c(alpha.init[i,1],alpha.init[i,2]), do.par= TRUE, nboot=0)
## If nboot = 0, the function ciregic provides the variance estimator of the regression 
## Parameter estimates using the least-squares method and does not perform the bootstrap method.
Alist[[i]]=fit
est.coef[i,]=fit$coef
nb=length(fit$gamma)     #total numbeer of B-spline bases functions =22;
cat("Number of B-spline basis functions used in each alpha choice=", nb, "\n")
Q[i,]<-c(i,fit$loglikelihood)
BIC.val=-2*fit$loglikelihood+log(ns)*(2*2+nb+(m+1)*2)
BIC[i,]<-c(i,BIC.val)
cat("This is the", i, "th iteration", "alpha.init=", c(alpha.init[i,1],alpha.init[i,2]),"loglik=", fit$loglikelihood, "BIC=", BIC.val, "\n \n") 
}

aloglik=cbind(alpha.init,Q[,2],BIC[,2])
loglik.maxind=which.max(Q[,2])
BIC.minind=which.min(BIC[,2])

alpha.optm.loglik=alpha.init[loglik.maxind,]
loglik.max=Q[,2][loglik.maxind]
alpha.optm.BIC=alpha.init[BIC.minind,]
BIC.min=BIC[,2][BIC.minind]

cat("number of Bernstein polynomial basis functions (m+1)=", m+1, "\n \n")
cat("number of B-spline basis functions nb               =", nb, "\n \n")
cat("number of bootstrap samples       nbt               =", nbt, "\n \n")

cat("optimal alpha selected by loglik=", as.numeric(alpha.optm.loglik), "\n")
cat("the corresponding loglik value  =", loglik.max, "\n \n")

cat("optimal alpha selected by BIC =", as.numeric(alpha.optm.BIC), "\n")
cat("the corresponding BIC value   =", BIC.min, "\n \n")


BIC.fit=Alist[[BIC.minind]]
print(summary(BIC.fit))

#agew1=cbind(Age1,Age2,Age3,Age4,Age5)%*%fit$coef[3:7]
#agew2=cbind(Age1,Age2,Age3,Age4,Age5)%*%fit$coef[10:14]
#agew1=cbind(Age1,Age2,Age3)%*%fit$coef[3:(2+m+1)]
#agew2=cbind(Age1,Age2,Age3)%*%fit$coef[(2+m+1+3):((2+m+1)*2)]

BIC.est.coef=est.coef[BIC.minind,]

agew1=cbind(Age1,Age2,Age3)%*%BIC.est.coef[3:(2+m+1)]
agew2=cbind(Age1,Age2,Age3)%*%BIC.est.coef[(2+m+1+3):((2+m+1)*2)]

age.curve=data.frame(age,agew1,agew2)
age.order=order(age)
age.curve=age.curve[age.order,]

dev.new()

pdf(file="age-curve-nbp3-3.pdf")
par(mfrow=c(2,2))
par(mfrow=c(1,1))

plot(age.curve[,1], age.curve[,2],"l", lty=1,lwd=2, col="black",ylim=c(-2,2),ylab=expression(paste(phi[1],"(age)"," and ",phi[2],"(age)")),main=expression(paste("Estimation of risk functions ", phi[1], " and ", phi[2])),xlab="age")
lines(age.curve[,1], age.curve[,3],"l",lty=2,lwd=2)
legend(16,2,c(expression(paste(phi[1],"(age) for Loss to care")), expression(paste(phi[2],"(age) for Death"))), lty=c(1,2), cex=0.8) 

dev.off()

#llh <- rbind(alpha.init,Q[,2])
#cat(llh[which.max(llh),])

cat("The program ended: ", date(), "\n \n")
sink()



###############Output of analysis (running the motivating HIV "iedea" data ###########

#> source("Gridsearch.R")
#The program started: Tue May 12 08:07:30 2020 
#This is the 1 th iteration alpha.init= 0 0 loglik= -3143.516 
#This is the 2 th iteration alpha.init= 0.5 0 loglik= -3142.029 
#This is the 3 th iteration alpha.init= 1 0 loglik= -3147.804 
#This is the 4 th iteration alpha.init= 0 0.5 loglik= -3143.278 
#This is the 5 th iteration alpha.init= 0.5 0.5 loglik= -3142.253 
#This is the 6 th iteration alpha.init= 1 0.5 loglik= -3141.517 
#This is the 7 th iteration alpha.init= 0 1 loglik= -3143.147 
#This is the 8 th iteration alpha.init= 0.5 1 loglik= -3142.493 
#This is the 9 th iteration alpha.init= 1 1 loglik= -3141.647 
#optimal alpha selected=         1 0.5 
#the corresponding loglik value= -3141.517 
#Call:
#ciregic.default(formula = Surv2(v, u, event = c) ~ male + art_cd4 + 
#    Age1 + Age2 + Age3, data = repro, alpha = c(alpha.init[i, 
#    1], alpha.init[i, 2]), do.par = TRUE, nboot = 0)

#Event type 1
#       Estimate Std. Error z value Pr(>|z|)    
#male      0.2300     0.1041   2.210   0.0271 *  
#art_cd4   0.0000     0.0003  -0.016   0.9876    
#Age1      1.0337     0.1576   6.559   <2e-16 ***
#Age2     -1.0135     0.4360  -2.324   0.0201 *  
#Age3      0.0607     0.4441   0.137   0.8913    
#---
#Signif. codes:  
#0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#Event type 2
#        Estimate Std. Error z value Pr(>|z|)    
#male      0.5447     0.1915   2.845   0.0044 ** 
#art_cd4  -0.0035     0.0009  -3.739   0.0002 ***
#Age1     -1.2150     0.3515  -3.457   0.0005 ***
#Age2      1.2845     0.9390   1.368   0.1713    
#Age3      0.0115     0.8046   0.014   0.9886    
#---
#Signif. codes:  
#0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#The program ended:  Tue May 12 08:13:01 2020 



##Run "Gridsearch.R" with nboot=0 to get alpha values, 
##then run this program "Optimal.R" with nboot=500 to get SE;

#The program started: Sat Aug 22 07:21:59 2020 
#Number of B-spline basis functions used in each alpha choice= 38 
#This is the 1 th iteration alpha.init= 2 0 loglik= -3140.884 BIC= 6666.914 
 
#number of Bernstein polynomial basis functions (m+1)= 3 
 
#number of B-spline basis functions nb               = 38 
 
#number of bootstrap samples       nbt               = 500 
 
#optimal alpha selected by loglik= 2 0 
#the corresponding loglik value  = -3140.884 
 
#optimal alpha selected by BIC = 2 0 
#the corresponding BIC value   = 6666.914 
 
#Call:
#ciregic.default(formula = Surv2(v, u, event = c) ~ male + art_cd4 + 
#   Age1 + Age2 + Age3, data = repro, alpha = c(alpha.init[i, 
#   1], alpha.init[i, 2]), do.par = TRUE, nboot = nbt)

#Event type 1
#       Estimate Std. Error z value Pr(>|z|)    
#male      0.2822     0.1228   2.298   0.0216 *  
#art_cd4   0.0000     0.0004  -0.076   0.9398    
#Age1      1.3532     0.2541   5.325   <2e-16 ***
#Age2     -0.9919     0.6325  -1.568   0.1168    
#Age3      0.1847     0.4530   0.408   0.6834    
#---
#Signif. codes:  
#0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#Event type 2
#       Estimate Std. Error z value Pr(>|z|)   
#male      0.5173     0.1785   2.898   0.0038 **
#art_cd4  -0.0034     0.0015  -2.233   0.0255 * 
#Age1     -1.0384     0.3653  -2.842   0.0045 **
#Age2      1.4166     0.8611   1.645   0.0999 . 
#Age3      0.1677     0.6887   0.244   0.8076   
#---
#Signif. codes:  
#0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1



