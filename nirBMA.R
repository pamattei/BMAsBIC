library(chemometrics)
library(sBIC)
library(FactoMineR)
source("GSPPCAfunctions.R")
set.seed(123)
M=30
N=20
p=50

data(NIR)
x=NIR$xNIR
n=nrow(x)
## Preprocessing with GSPPCA
dv=estim_ncp(scale(x),method = "Smooth")$ncp
u=VEM(scale(x),dv,alpha=100,epsi=1e-4)$u

## Load preprocessed data
x=scale(x)
dat=list(X=t(x[,order(u,decreasing=TRUE)[(N+1):p]]), Y=t(x[,order(u,decreasing=TRUE)[1:N]]))
X=t(dat$X)[,]
Y=t(dat$Y)[,]


## Do the experiment
nrep=1000
BIC=numeric(nrep)
BMA_BIC=numeric(nrep)
sBIC=numeric(nrep)
BMA_sBIC=numeric(nrep)
OLS=numeric(nrep)
rreg = ReducedRankRegressions(numResponses=N, numCovariates=M, maxRank=min(M, N))

for (k0 in 1:nrep){
  train=sample(1:n,floor(n*0.5))
  test=(1:n)[-train]
  Xt=X[train,]
  Yt=Y[train,]
  Bols=chol2inv(chol(crossprod(Xt)))%*%t(Xt)%*%Yt
  
  
  dat<-list(X=t(Xt),Y=t(Yt))
  res = sBIC(dat, rreg)
  
  v=svd(Xt%*%Bols)$v
  vbic=v[,1:(which.max(res$BIC) - 1)]
  vsbic=v[,1:(which.max(res$sBIC) - 1)]
  bicpost=BICpost(res$BIC)
  sbicpost=BICpost(res$sBIC)
  B1=Bols%*%tcrossprod(vbic)
  B2=Bols%*%tcrossprod(vsbic)
  Bav1=matrix(0,ncol=N,nrow=M)
  for (i in which(bicpost>0)){
    if (i>1) Bav1=Bav1+bicpost[i]*Bols%*%tcrossprod(v[,(1:(i-1))])
  }
  Bav2=matrix(0,ncol=N,nrow=M)
  for (i in which(sbicpost>0)){
    if (i>1) Bav2=Bav2+sbicpost[i]*Bols%*%tcrossprod(v[,(1:(i-1))])
  }
  
  
  BIC[k0]=sum((Y[test,]-X[test,]%*%B1)^2)/length(test)
  sBIC[k0]=sum((Y[test,]-X[test,]%*%B2)^2)/length(test)
  BMA_BIC[k0]=sum((Y[test,]-X[test,]%*%Bav1)^2)/length(test)
  BMA_sBIC[k0]=sum((Y[test,]-X[test,]%*%Bav2)^2)/length(test)
  OLS[k0]=sum((Y[test,]-X[test,]%*%Bols)^2)/length(test)
}

boxplot(BIC,sBIC,BMA_BIC,BMA_sBIC,OLS)
boxplot(BMA_BIC,BMA_sBIC)
mean(BMA_BIC)

t(lapply(list(BIC,sBIC,BMA_BIC,BMA_sBIC,OLS),mean))
t(lapply(list(BIC,sBIC,BMA_BIC,BMA_sBIC,OLS),sd))





