
setwd("C:/Users/o0/Desktop/ordinal data/simulation/MLE-for-cs-samples")
source("function.R")

##--------------------------------------some test--------------------------------------------------------
eta.pro=c(1,2,3,1)
likelihood.pro(eta.pro,obs)
likelihood.pro.grad(eta.pro,obs)
cpy=colSums(obs)/sum(obs)
eta.mod=c(eta.pro,rep(cpy[1],3))
likelihood.mod(eta.mod,obs)
likelihood.mod.grad(eta.mod,obs)
#gennerate ccdata
lrm(y~g,data=Popdata)$coefficients
ccdata=cc_data(n,Popdata)
obs=table(ccdata$g,ccdata$y)
sta.pro=c(1,2,3,1)
system.time({
  p11=constrOptim(sta.pro,likelihood.pro, grad=NULL, ui.pro, ci.pro, control=list(fnscale=-1), outer.iterations=100000, obs=obs)
})
system.time({
  p12=constrOptim(sta.pro,likelihood.pro, grad=likelihood.pro.grad, ui.pro, ci.pro, control=list(fnscale=-1), outer.iterations=100000, obs=obs)
})
system.time({
  p21=maxLik(likelihood.pro,grad=NULL,start=sta.pro,constraints = list(ineqA=ui.pro,ineqB=ci.pro),obs=obs)
})
system.time({
  p22=maxLik(likelihood.pro,grad=likelihood.pro.grad,start=sta.pro,constraints = list(ineqA=ui.pro,ineqB=ci.pro),obs=obs)
})

pr=cbind(c(p11$par,p11$value),c(p12$par,p12$value),c(p21$estimate,p21$maximum),c(p22$estimate,p22$maximum))
##init value
sta.mod=c(1,2,3,1,5,6,7)
system.time({
  re11=constrOptim(sta.mod,likelihood.mod, grad=NULL, ui.mod, ci.mod, control=list(fnscale=-1), outer.iterations=100000, obs=obs)
})
#variance is less than without grad
system.time({
  re12=constrOptim(sta.mod,likelihood.mod, grad=likelihood.mod.grad, ui.mod, ci.mod, control=list(fnscale=-1), outer.iterations=100000, obs=obs)
})
likelihood.mod.grad(re11$par,obs)
likelihood.mod.grad(re12$par,obs)
##not correct as above
system.time({
  re21=maxLik(likelihood.mod,grad=NULL,start=sta.mod,constraints = list(ineqA=ui.mod,ineqB=ci.mod),obs=obs)
})
##sometime it is not converage
system.time({
  re22=maxLik(likelihood.mod,grad=likelihood.mod.grad,start=sta.mod,constraints = list(ineqA=ui.mod,ineqB=ci.mod),obs=obs)
})
likelihood.mod.grad(re21$estimate,obs)
likelihood.mod.grad(re22$estimate,obs)
result=cbind(c(sta.mod,1),c(re11$par,re11$value),c(re12$par,re12$value),c(re21$estimate,re21$maximum),c(re22$estimate,re22$maximum))

#true lamda
py=table(Popdata$y)/N
cpy=colSums(obs)/sum(obs)
la=py[1]*cpy/py


##--------------------------------------stimulation--------------------------------------------------------
library(rms)
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
sta.pro=c(3,4,5,1)
ui.pro=matrix(0,J-2,J)
for (i in 1:(J-2)){
  ui.pro[i,i]=-1
  ui.pro[i,i+1]=1
}
ci.pro=rep(0,J-2)
# sta.mod=c(1,2,3,1,1,2,3)
# sta.mod=c(1,2,3,2,1,2,3)
sta.mod=c(1,4,5,1,5,6,7)
# sta.mod=c(1,4,5,2,1,4,5)
T=1000
# pro=matrix(0,T,J)
pro0=matrix(0,T,J)
pro1=matrix(0,T,J)
mod0=matrix(0,T,2*J-1)
mod1=matrix(0,T,2*J-1)
ptm=proc.time()
for (k in 1:T){
  Popdata=gen_popdata(theta,beta,p,N)
  ccdata=cc_data(n,Popdata)
  obs=table(ccdata$g,ccdata$y)
  pro0[k,]=constrOptim(sta.pro,likelihood.pro, grad=NULL, ui.pro, ci.pro, control=list(fnscale=-1), outer.iterations=100000, obs=obs)$par
  mod0[k,]=constrOptim(sta.mod,likelihood.mod, grad=NULL, ui.mod, ci.mod, control=list(fnscale=-1), outer.iterations=100000, obs=obs)$par
  mod1[k,]=constrOptim(sta.mod,likelihood.mod, grad=likelihood.mod.grad, ui.mod, ci.mod, control=list(fnscale=-1), outer.iterations=100000, obs=obs)$par
}
proc.time()-ptm


r=cbind(mod0[,J],mod1[,J])
boxplot(pro1[,J],mod0[,J],mod1[,J], at = 1:3,col=c("2","3","5"),boxwex = 0.5,xlim=c(0,20),xaxt = "n")
abline(h=beta,lty = 2)
boxplot(pro0[,J],mod0[,J],mod1[,J], at = 1:3+4 ,col=c("2","3","5"),add=T,boxwex = 0.5,xaxt = "n")
boxplot(pro0[,J],mod0[,J],mod1[,J], at = 1:3+8,col=c("2","3","5"),add=T,boxwex = 0.5,xaxt = "n")
boxplot(pro[,J],mod0[,J],mod1[,J], at = 1:3+12,col=c("2","3","5"),add=T,boxwex = 0.5,xaxt = "n")
# axis(1, at = seq(2,10,4), labels = c("1","2","3"), tick = TRUE)
# legend("topright", c("constrOptim","constrOptim.grad","maxLik"), fill = c("2", "3","5"))

pbias=mean((pro[,1]-beta))
mbias=mean((mod[,1]-beta))
nbias=mean((newb[,1]-beta))
pmse=mean((pro[,1]-beta)^2)
mmse=mean((mod[,1]-beta)^2)
nmse=mean((newb[,1]-beta)^2)
var(pro[,J])
var(mod[,J])
var(newb[,1])
#the true lamda
py=table(Popdata$y)/N
cpy=colSums(ccdata)/sum(ccdata)
la=py[J]*cpy/py
