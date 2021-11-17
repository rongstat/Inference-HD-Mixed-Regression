
library(MASS)
library(purrr)
library(glmnet)
library(flare)
library(HDclassif)
library(xLLiM)


p=600 #600-1000

sp=c(10,15,20,25,30)


omega0 = 0.3
T0=30
kappa=0.3
C=1.8


n=400

result=matrix(ncol=5,nrow=5)
for(sp.index in 1:5){

#generate data
b1 = rep(0,p)
b1[1:sp[sp.index]]=.85
b2 = rep(0,p)
b2[(p/2+1):(p/2+sp[sp.index])]=-.85
sig1 = toeplitz(seq(0.4,0,length.out = p/10))
Sig = bdiag(rep(list(sig1),10))+diag(rep(0.6,p))


X = mvrnorm(n, mu=rep(0,p), Sigma=Sig)
z=rbinom(n,1,prob=omega0)
y=c() 
for(i in 1:n){
  if(z[i]==1) y[i] = X[i,] %*% b1 + rnorm(1)
  if(z[i]==0) y[i] = X[i,] %*% b2 + rnorm(1)
}

#################
#####MIREM 1
#################

#initialization
my.fit=glmnet(x = X, y = y, family = "gaussian", alpha = 1,  intercept=F, 
              lambda = 5*sqrt(log(p)/n))$beta
support = which(my.fit!=0)

init.class=emgm(t(cbind(y,X[,support])),init=2)
while(sum(init.class$R[,1]!=0) < n/4 | sum(init.class$R[,1]!=0) > n*3/4){
  init.class=emgm(t(cbind(y,X[,support])),init=2)
}
init.class= init.class$label-1
my.fit=glmnet(x = X[init.class==0,], y = y[init.class==0], family = "gaussian", alpha = 0.5,  intercept=F, 
              lambda = 2*sqrt(log(p)/n))
b1.01 = my.fit$beta
my.fit=glmnet(x = X[init.class==1,], y = y[init.class==1], family = "gaussian", alpha = 0.5,  intercept=F, 
              lambda = 2*sqrt(log(p)/n))
b2.01 = my.fit$beta

omega = 1/2
gamma = omega*exp(-(y-X%*%b1.01)^2/2)/(omega*exp(-(y-X%*%b1.01)^2/2)+(1-omega)*exp(-(y-X%*%b2.01)^2/2))

#start EM
lambda=sqrt(log(p)/n)
for(t in 1:(T0-1)){
  b1.temp1 = slim(verbose=FALSE,X=diag(as.vector(sqrt(gamma)),nrow=length(gamma)) %*% X, Y=diag(as.vector(sqrt(gamma)),nrow=length(gamma)) %*% y,
                  lambda = lambda, method = "lasso")$beta
  b2.temp1 = slim(verbose=FALSE,X=diag(as.vector(sqrt(1-gamma)),nrow=length(gamma)) %*% X, Y=diag(as.vector(sqrt(1-gamma)),nrow=length(gamma)) %*% y,
                  lambda = lambda, method = "lasso")$beta
  sigma1.hat =mean( (diag(as.vector(sqrt(gamma)),nrow=length(gamma)) %*% y-diag(as.vector(sqrt(gamma)),nrow=length(gamma)) %*% X %*% b1.temp1)^2)
  sigma2.hat =mean((diag(as.vector(sqrt(1-gamma)),nrow=length(gamma)) %*% y-diag(as.vector(sqrt(1-gamma)),nrow=length(gamma)) %*% X %*% b2.temp1)^2)
  sigma.hat = sqrt((sigma1.hat+sigma2.hat)/2)
  lambda = kappa*lambda+C*sqrt(log(p)/n)
  gamma = omega*exp(-(y-X%*%b1.temp1)^2/2/sigma.hat^2)/(omega*exp(-(y-X%*%b1.temp1)^2/2/sigma.hat^2)+(1-omega)*exp(-(y-X%*%b2.temp1)^2/2/sigma.hat^2))
  gamma[which(is.na(gamma))] = 0
  omega = mean(gamma)
}


#################
#####MIREM 2
#################

#initialization
init.class2=hddc((cbind(y,X[,support])),K=2,show=FALSE)$class-1
my.fit=glmnet(x = X[init.class2==0,], y = y[init.class2==0], family = "gaussian", alpha = 0.5,  intercept=F, 
              lambda = 2*sqrt(log(p)/n))
b1.02 = my.fit$beta
my.fit=glmnet(x = X[init.class2==1,], y = y[init.class2==1], family = "gaussian", alpha = 0.5,  intercept=F, 
              lambda = 2*sqrt(log(p)/n))
b2.02 = my.fit$beta

omega = 1/2
gamma = omega*exp(-(y-X%*%b1.02)^2/2)/(omega*exp(-(y-X%*%b1.02)^2/2)+(1-omega)*exp(-(y-X%*%b2.02)^2/2))

#start EM
lambda=sqrt(log(p)/n)
for(t in 1:(T0-1)){
  b1.temp2 = slim(verbose=FALSE,X=diag(as.vector(sqrt(gamma)),nrow=length(gamma)) %*% X, Y=diag(as.vector(sqrt(gamma)),nrow=length(gamma)) %*% y,
                  lambda = lambda, method = "lasso")$beta
  b2.temp2 = slim(verbose=FALSE,X=diag(as.vector(sqrt(1-gamma)),nrow=length(gamma)) %*% X, Y=diag(as.vector(sqrt(1-gamma)),nrow=length(gamma)) %*% y,
                  lambda = lambda, method = "lasso")$beta
  lambda = kappa*lambda+C*sqrt(log(p)/n)
  
  sigma1.hat =mean( (diag(as.vector(sqrt(gamma)),nrow=length(gamma)) %*% y-diag(as.vector(sqrt(gamma)),nrow=length(gamma)) %*% X %*% b1.temp1)^2)
  sigma2.hat =mean((diag(as.vector(sqrt(1-gamma)),nrow=length(gamma)) %*% y-diag(as.vector(sqrt(1-gamma)),nrow=length(gamma)) %*% X %*% b2.temp1)^2)
  sigma.hat = sqrt(sigma1.hat+sigma2.hat)
  gamma = omega*exp(-(y-X%*%b1.temp2)^2/2/sigma.hat^2)/(omega*exp(-(y-X%*%b1.temp2)^2/2/sigma.hat^2)+(1-omega)*exp(-(y-X%*%b2.temp2)^2/2/sigma.hat^2))
  gamma[which(is.na(gamma))] = 0
  omega = mean(gamma)
}


###################
####GLLiM
####################

b.gllim=gllim(tapp=t(y),t(X),in_K=2)$A
b.gllim1=b.gllim[1:p]
b.gllim2=b.gllim[(p+1):(2*p)]



error.est1 = min(c(sqrt(sum((b.gllim1-b1)^2))+ sqrt(sum((b.gllim2-b2)^2)),sqrt(sum((b.gllim1-b2)^2))+ sqrt(sum((b.gllim2-b1)^2))))
error.est3 = min(c(sqrt(sum((b1.02-b1)^2))+ sqrt(sum((b2.02-b2)^2)),sqrt(sum((b1.02-b2)^2))+ sqrt(sum((b2.02-b1)^2))))
error.est2 = min(c(sqrt(sum((b1.01-b1)^2))+ sqrt(sum((b2.01-b2)^2)),sqrt(sum((b1.01-b2)^2))+ sqrt(sum((b2.01-b1)^2))))
error.est5 = min(c(sqrt(sum((b1.temp2-b1)^2))+ sqrt(sum((b2.temp2-b2)^2)),sqrt(sum((b1.temp2-b2)^2))+ sqrt(sum((b2.temp2-b1)^2))))
error.est4 = min(c(sqrt(sum((b1.temp1-b1)^2))+ sqrt(sum((b2.temp1-b2)^2)),sqrt(sum((b1.temp1-b2)^2))+ sqrt(sum((b2.temp1-b1)^2))))


result.temp = c(error.est1,error.est2,error.est3,error.est4,error.est5)
result[sp.index,] = round(result.temp,3)

}

save(result, file= "Simu_est.RData")