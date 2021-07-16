## R program for HIV data analysis in the <<Statistica Sinica>> paper titled
## "Sieve estimation of a class of partially linear transformation models 
## with interval-censored competing risks data";
## Author: Xuewen Lu; Edited by Dipankar Bandyopadhyay
## Date: July 16, 2021


#Run this program "Gridsearch.R", to get an idea for optimal alpha #value; set this program with nboot=0 (runs faster).
#With grid search, we obtained optimal alpha=(2,0), with the maximum log-likelihood value.
#Then setting alph=(2,0), run "Optimal.R" with nboot=500 to get more accurate SE estimates.
#See the R output results at the end;
#We found m=2 as the best choice based on BIC (for our dataset) 
#Our HIV dataset is not uploaded in GitHub for security reasons


rm(list=ls())
cat("The program started:", date(), "\n")
options(width=60,length=200)
set.seed(20835)
library(basefun)  
library(foreach)
library(doParallel)
library(intccr)

#Function to generate (m+1) Bernstein polynomial basis functions;

m=2
Lambda.basis<- function(w,m){ 	
  minTX=min(w)
  maxTX=max(w)
  #bb.TX <- Bernstein_basis(numeric_var("x", support = c(minTX, maxTX)), order = 4, ui = "increasing")       #For m=4, 5 basis functions;
  bb.TX <- Bernstein_basis(numeric_var("x", support = c(minTX, maxTX)), order = m, ui = "increasing")       #For m=1, 2 basis functions;
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

#cbind(repro,Age1,Age2, Age3) -> repro

#alpha1 <- seq(12.5,50,0.5)
#alpha2 <- seq(12.5,50,0.5)

#alpha1 <- seq(0,1,0.5)
#alpha2 <- seq(0,1,0.5)


alpha1 <- seq(0,4,0.5)
alpha2 <- seq(0,4,0.5)


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
fit <- ciregic(formula=Surv2(v, u, event=c)~male+art_cd4+Age1+Age2+Age3, data=repro, alpha=c(alpha.init[i,1],alpha.init[i,2]), k=1, do.par= TRUE, nboot=0)
#fit <- ciregic(formula=Surv2(v, u, c)~male+art_cd4+Age1, data=repro, alpha=c(alpha.init[i,1],alpha.init[i,2]), do.par= TRUE, nboot=0)
## If nboot = 0, the function ciregic provides the variance estimator of the regression 
## parameter estimates using the least-squares method and does not perform the bootstrap method.
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



###############Output of data analysis###########
##> source("Gridsearch.R")
#Number of B-spline basis functions used in each alpha choice= 38 
#This is the 81 th iteration alpha.init= 4 4 loglik= -3142.56 BIC= 6670.267 
 
#number of Bernstein polynomial basis functions (m+1)= 3 
 
#number of B-spline basis functions nb               = 38 
 
#optimal alpha selected by loglik= 2 0 
#the corresponding loglik value  = -3140.884 
 
#optimal alpha selected by BIC = 2 0 
# Call:
#ciregic.default(formula = Surv2(v, u, event = c) ~ male + art_cd4 + 
#    Age1 + Age2 + Age3, data = repro, alpha = c(alpha.init[i,1], 
#    alpha.init[i, 2]), k = 1, do.par = TRUE, nboot = 0)

#Event type 1
#        Estimate Std. Error z value Pr(>|z|)    
#male      0.2822     0.1206   2.339   0.0193 *  
#art_cd4   0.0000     0.0004  -0.083   0.9342    
#Age1      1.3532     0.1841   7.350   <2e-16 ***
#Age2     -0.9919     0.5016  -1.978   0.0480 *  
#Age3      0.1847     0.5039   0.367   0.7139    
#---
#Signif. codes:  
#0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#Event type 2
#        Estimate Std. Error z value Pr(>|z|)    
#male      0.5173     0.1849   2.797   0.0052 ** 
#art_cd4  -0.0034     0.0009  -3.606   0.0003 ***
#Age1     -1.0384     0.3427  -3.030   0.0024 ** 
#Age2      1.4166     0.9100   1.557   0.1195    
#Age3      0.1677     0.7747   0.216   0.8286    
#---
#Signif. codes:  
#0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#The program ended:  Tue May 12 17:52:15 2020 
