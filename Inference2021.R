args= commandArgs(trailingOnly=TRUE)

set.seed(c(as.numeric(args[2])))

library(MASS)
library(purrr)
library(glmnet)
library(flare)
library(ggplot2)
library(HDclassif)
library(xLLiM)


p=c(800,850,900,950,1000)

sp=15 


omega0 = 0.3
eta = 0.3 #try eta=0.5 if doesn't work
T0=30
kappa=0.3
C=1.8



n=300
result = matrix(ncol=6,nrow=5)
for(o in 1:5){
  
  #generate data
  b1 = rep(0,p[o])
  b1[1:sp[sp.index]]=.45
  b2 = rep(0,p[o])
  b2[(p[o]/2+1):(p[o]/2+sp[sp.index])]=-.45
  sig1 = toeplitz(seq(0.4,0,length.out = p[o]/10))
  Sig = bdiag(rep(list(sig1),10))+diag(rep(0.6,p[o]))
  
  
  X = mvrnorm(n, mu=rep(0,p[o]), Sigma=Sig)
  z=rbinom(n,1,prob=omega0)
  y=c() 
  for(i in 1:n){
    if(z[i]==1) y[i] = X[i,] %*% b1 + rnorm(1)
    if(z[i]==0) y[i] = X[i,] %*% b2 + rnorm(1)
  }
  
  ####################
  ######MIREM1
  ####################
  
  #initialization
  my.fit=glmnet(x = X, y = y, family = "gaussian", alpha = 1,  intercept=F, 
                lambda = 5*sqrt(log(p[o])/n))$beta
  support = which(my.fit!=0)
  
  init.class=emgm(t(cbind(y,X[,support])),init=2)
  while(sum(init.class$R[,1]!=0) < n/4 | sum(init.class$R[,1]!=0) > n*3/4){
    init.class=emgm(t(cbind(y,X[,support])),init=2)
  }
  init.class= init.class$label-1
  my.fit=glmnet(x = X[init.class==0,], y = y[init.class==0], family = "gaussian", alpha = 0.5,  intercept=F, 
                lambda = 2*sqrt(log(p[o])/n))
  b1.01 = my.fit$beta
  my.fit=glmnet(x = X[init.class==1,], y = y[init.class==1], family = "gaussian", alpha = 0.5,  intercept=F, 
                lambda = 2*sqrt(log(p[o])/n))
  b2.01 = my.fit$beta
  
  omega = 1/2
  gamma = omega*exp(-(y-X%*%b1.01)^2/2)/(omega*exp(-(y-X%*%b1.01)^2/2)+(1-omega)*exp(-(y-X%*%b2.01)^2/2))
  
  
  #start EM
  lambda=sqrt(log(p[o])/n)
  for(t in 1:(T0-1)){
    b1.temp1 = slim(verbose=FALSE,X=diag(as.vector(sqrt(gamma)),nrow=length(gamma)) %*% X, Y=diag(as.vector(sqrt(gamma)),nrow=length(gamma)) %*% y,
                    lambda = lambda, method = "lasso")$beta
    b2.temp1 = slim(verbose=FALSE,X=diag(as.vector(sqrt(1-gamma)),nrow=length(gamma)) %*% X, Y=diag(as.vector(sqrt(1-gamma)),nrow=length(gamma)) %*% y,
                    lambda = lambda, method = "lasso")$beta
    sigma1.hat =mean( (diag(as.vector(sqrt(gamma)),nrow=length(gamma)) %*% y-diag(as.vector(sqrt(gamma)),nrow=length(gamma)) %*% X %*% b1.temp1)^2)
    sigma2.hat =mean((diag(as.vector(sqrt(1-gamma)),nrow=length(gamma)) %*% y-diag(as.vector(sqrt(1-gamma)),nrow=length(gamma)) %*% X %*% b2.temp1)^2)
    sigma.hat = sqrt((sigma1.hat+sigma2.hat)/2)
    lambda = kappa*lambda+C*sqrt(log(p[o])/n)
    gamma = omega*exp(-(y-X%*%b1.temp1)^2/2/sigma.hat^2)/(omega*exp(-(y-X%*%b1.temp1)^2/2/sigma.hat^2)+(1-omega)*exp(-(y-X%*%b2.temp1)^2/2/sigma.hat^2))
    gamma[which(is.na(gamma))] = 0
    omega = mean(gamma)
  }
  
  #inference & bias correction
  #M = sugm(verbose =FALSE,data= X,method = "clime",lambda=eta*sqrt(log(p[o])/n ))
M=M.safe
  Cov1=t(diag(as.vector(sqrt(gamma)),
              nrow=length(gamma)) %*% X) %*% diag(as.vector(sqrt(gamma)),
                                                  nrow=length(gamma)) %*% X
  Cov1=Cov1/n
  SIGxy1 = t(diag(as.vector(gamma),nrow=length(gamma)) %*% X) %*% y
  SIGxy1=SIGxy1/n
  b1.hat = b1.temp1+M$icov[[1]] %*% (SIGxy1 - Cov1 %*% b1.temp1)
  Cov2=t(diag(as.vector(sqrt(1-gamma)),nrow=length(gamma)) %*% X) %*% diag(as.vector(sqrt(1-gamma)),nrow=length(gamma)) %*% X
  Cov2=Cov2/n
  SIGxy2 = t(diag(as.vector(1-gamma),nrow=length(gamma)) %*% X) %*% y
  SIGxy2=SIGxy2/n
  b2.hat = b2.temp1+M$icov[[1]] %*% (SIGxy2 - Cov2 %*% b2.temp1)
  
  
  #variance estimation
  T11.1 = -Cov1
  T22.1 = -Cov2
  
  T11.w= 2*(y- X%*% b1.temp1)^2/sigma.hat^2/(omega+(1-omega)*exp((2*y-(X%*%(b1.temp1+b2.temp1)))*(X%*% (b2.temp1-b1.temp1))/2/sigma.hat^2))/(omega+(1-omega)*exp(-(2*y-(X%*%(b1.temp1+b2.temp1)))*(X%*% (b2.temp1-b1.temp1))/2/sigma.hat^2))
  T11.x = diag(as.vector(T11.w)) %*% X
  T11 = (t(X) %*% T11.x)/n
  
  T22.w= 2*(y- X%*% b2.temp1)^2/sigma.hat^2/(omega+(1-omega)*exp((2*y-(X%*%(b1.temp1+b2.temp1)))*(X%*% (b2.temp1-b1.temp1))/2/sigma.hat^2))/(omega+(1-omega)*exp(-(2*y-(X%*%(b1.temp1+b2.temp1)))*(X%*% (b2.temp1-b1.temp1))/2/sigma.hat^2))
  T22.x = diag(as.vector(T22.w)) %*% X
  T22 = (t(T22.x) %*% X)/n
  
  V1 = T11-T11.1
  v1 = diag(M$icov[[1]]%*% V1 %*% M$icov[[1]])
  V2 = T22-T22.1
  v2 = diag(M$icov[[1]]%*% V2 %*% M$icov[[1]])
  
  #obtain z scores
  T1 = sqrt(n)*(b1.hat)/sqrt(v1)
  T2 = sqrt(n)*(b2.hat)/sqrt(v2)
  
  
  #FDR procedure
  cutoff = function(x,y) p[o]*(1-pchisq(x^2, df=1))/max(1,sum(abs(T1) > x | abs(T2) > x))-y
  alpha=0.1
  t = sqrt(2*log(p[o]))
  x= seq(0, sqrt(2*log(p[o])-2*log(log(p[o]))), 0.001)
  for(k in 1:length(x)){
    if(cutoff(x[k], alpha/2)<0) t = min(x[k],t)
  }
  test.m1=rep(1,p[o])
  test.m1[abs(T1)<t & abs(T2)<t] = 0
  
  
  
  ########################
  ######### BY procedure
  ########################
  
  BY = (p.adjust(2-2*pnorm(abs(T1)),method="BY")<alpha/2)|(p.adjust(2-2*pnorm(abs(T2)),method="BY")<alpha/2)
  
  
  ###################
  ###########dLasso
  #####################
  b.lasso = as.numeric(glmnet(x = X, y = y, family = "gaussian", alpha = 1,  intercept=F, 
                              lambda = 5*sqrt(log(p[o])/n))$beta)
  b.dlasso = b.lasso+M$icov[[1]]%*%((t(X)%*% y)/n - ((t(X)%*%X)/n) %*% b.lasso)
  v.lasso = diag(M$icov[[1]]%*% ((t(X)%*%X)/n) %*% M$icov[[1]])*sigma.hat^2
  T.lasso = sqrt(n)*(b.dlasso)/sqrt(v.lasso)
  
  cutoff = function(x) p[o]*(1-pchisq(x^2, df=1))/max(1,sum(abs(T.lasso) > x))
  t = sqrt(2*log(p[o]))
  x= seq(0, sqrt(2*log(p[o])-2*log(log(p[o]))), 0.001)
  for(k in 1:length(x)){
    if(cutoff(x[k])<alpha) t = min(x[k],t)
  }
  
  test.la=rep(1,p[o])
  test.la[abs(T.lasso)<t] = 0
  
  
  
  ###########################
  ####evaluate power and fdr
  ############################
  
  power.temp1 = sum(test.m1[c(1:sp[sp.index],(p[o]/2+1):(p[o]/2+sp[sp.index]))])/2/sp[sp.index]
  fdr.temp1 = sum(test.m1[-c(1:sp[sp.index],(p[o]/2+1):(p[o]/2+sp[sp.index]))])/(sum(test.m1==1))
  
  power.temp.by1 = sum(BY[c(1:sp[sp.index],(p[o]/2+1):(p[o]/2+sp[sp.index]))])/2/sp[sp.index]
  fdr.temp.by1 = sum(BY[-c(1:sp[sp.index],(p[o]/2+1):(p[o]/2+sp[sp.index]))])/(length(BY))

  power.la.temp= sum(test.la[c(1:sp[sp.index],(p[o]/2+1):(p[o]/2+sp[sp.index]))])/2/sp[sp.index]
  fdr.la.temp = sum(test.la[-c(1:sp[sp.index],(p[o]/2+1):(p[o]/2+sp[sp.index]))])/(sum(test.la==1))
  
  
  result[o,] = c(power.temp1,power.la.temp, power.temp.by1,
                 fdr.temp1,fdr.la.temp,fdr.temp.by1)
  
}

save(result, file= "Simu_inf.RData")