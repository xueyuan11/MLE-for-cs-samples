generate.data = function( alphay, betay, N) {
  z = rnorm(N,0,1)
  y = numeric(N)
  py = (1 + exp(- outer(alphay, betay*z ,"+"))) ^ (-1)
  aa = runif(N)
  for(i in 1:N)
    y[i] = sum(aa[i] > py[,i])
    y = as.numeric(as.factor(y))
    data.frame(y=y,z=z)
}

choicedata<-function(n,data,p){
  Nv<-table(data$y)
  ny=length(Nv)
  d<-list()
  f<-NULL
  for(i in 1:ny){ 
    if (Nv[i]<n*p[i]){
      cat("Error:","The population don not have enough individuals to choice in ",i,"\n")
      cat("Please increase more population or change p")
      s=0
      break
    }
    else
      s=1
  } 
  if (s==1){
    for (i in 1:ny){
      d[[i]]=subset(data,y==i)
      c<-sample(length(d[[i]]$y),n*p[i])
      d[[i]]<-d[[i]][c,]
      f<-rbind(f,d[[i]])
    }
  }
  return(f)
}

library(rms)
N = 1000000
NREPL = 1000
n=500
# 1. X_5. Y_4
alphay = c(3.89, 4.59, 5.51)
betay = 0
p=c(0.5,0.25,0.15,0.1)
beta=rep(0,NREPL)
system.time({data = generate.data(alphay, betay, N)})

ptm<-proc.time()
for (i in 1:NREPL){
 ccdata=choicedata(n,data,p)
 beta[i]=lrm(y~z,ccdata)$coefficients[4]
}
proc.time()-ptm
boxplot(beta)
d<-lrm(y~z,ccdata)
