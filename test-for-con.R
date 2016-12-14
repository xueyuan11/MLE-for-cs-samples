setwd("C:/Users/o0/Desktop/ordinal data/simulation/MLE-for-cs-samples")
source("function.R")

library(maxLik)
p=0.4
N=100000

theta=c(3.89,4.59,5.51)
beta=log(1.6)

n=1000



J=4
#likelihood.mod constraints
ui.mod=matrix(0,2*J-3,2*J-1)
for (i in 1:(J-2)){
  ui.mod[i,i]=-1
  ui.mod[i,i+1]=1
}
ui.mod[(J-1):(2*J-3),(J+1):(2*J-1)]=diag(J-1)
ci.mod=rep(0,2*J-3)
#likelihood.pro constraints
ui.pro=matrix(0,J-2,J)
for (i in 1:(J-2)){
  ui.pro[i,i]=-1
  ui.pro[i,i+1]=1
}

ci.pro=rep(0,J-2)


sta.pro=c(3,4,5,1)
sta.mod=c(3,4,5,1,5,6,7)

T=1000
pro0=matrix(0,T,J)
mod0=matrix(0,T,2*J-1)
mod1=matrix(0,T,2*J-1)

ptm=proc.time()
  for (k in 1:T){
    Popdata=gen_popdata(theta,beta,p,N)
    ccdata=cc_data(n,Popdata)
    obs=table(ccdata$g,ccdata$y)
    pro0[k,]=constrOptim(sta.pro,likelihood.pro, grad=NULL, ui.pro, ci.pro, control=list(fnscale=-1), outer.iterations=100000, obs=obs)$par
    mod0[k,]=constrOptim(sta.mod,likelihood.mod, grad=NULL, ui.mod, ci.mod, control=list(fnscale=-1), outer.iterations=100000, obs=obs)$par
}

proc.time()-ptm

boxplot(pro0[,J],mod0[,J], at = 1:2,col=c("2","3"),boxwex = 0.5,xlim=c(0,4),xaxt = "n")

boxplot(pro0[,J],pro1[,J],pro2[,J], at = 1:3,col=c("2","3","5"),boxwex = 0.5,xlim=c(0,20),xaxt = "n")
boxplot(mod0[,J],mod1[,J],mod2[,J], at = 1:3+4,col=c("2","3","5"),add=T,boxwex = 0.5,xaxt = "n")
boxplot(mod0[,J],mod1[,J],mod2[,J], at = 1:3+8,col=c("2","3","5"),add=T,boxwex = 0.5,xaxt = "n")
boxplot(mod[,1],mod1[,1],mod2[,1], at = 1:3+12,col=c("2","3","5"),add=T,boxwex = 0.5,xaxt = "n")
legend("topright", c("constrOptim","constrOptim.grad","maxLik"), fill = c("2", "3","5"))
abline(h=beta,lty =2)
abline(h=0.5481,lty = 2)
library(rms)

lrm(y~g,data=ccdata)
#the true lamda
py=table(Popdata$y)/N
cpy=colSums(obs)/sum(obs)
la=py[1]*cpy/py