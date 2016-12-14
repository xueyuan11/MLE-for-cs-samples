
setwd("C:/Users/o0/Desktop/ordinal data/simulation/MLE-for-cs-samples")
source("mle.R")
#---------------------------------generate data---------------------------------

##the population
gen_data = function( alpha, beta,p, N) {
  pG=c(p^2,2*p*(1-p),(1-p)^2)
  g = sample(0:2,N,replace = TRUE,prob = pG)
  y = numeric(N)
  py = (1 + exp(- outer(alpha, -beta*g ,"+"))) ^ (-1)
  aa = runif(N)
  for(i in 1:N)
    y[i] = sum(aa[i] > py[,i])
  y = as.numeric(as.factor(y))
  data.frame(y=y,g=g)
}
###the case-control samples from population
choicedata<-function(n,data){
  Nc<-table(data$y)
  Nv<-c(Nc[1],N-Nc[1])
  ncy=length(Nc)
  csc<-c("case","control")
  for(i in 1:2){ 
    if (Nv[i]<n){
      cat("Error:","The population don not have enough individuals to choice in ",csc[i],"/n")
      cat("Please increase more population ")
      s=0
      break
    }
    else
      s=1
  } 
  if (s==1){
      cs<-list()
      control=subset(data,y==1)
      cs[[1]]<-control[sample(Nv[1],n),2]
      case=subset(data,y>1)
      ca=case[sample(Nv[2],n),]
      for (j in 2:ncy){
        cs[[j]]=subset(ca,y==j)[,2]
      }
   
  }
  names(cs)<-c("y=1","y=2","y=3")
  return(cs)
}

p=0.4
N=100000

theta=c(3.48,4.6)
beta=log(1.4)

n=500
# py=c(0.5,0.33,0.17)
# ny=n*py


Popdata=gen_data(theta,beta,p,N)
ccdata=choicedata(n,Popdata)


#--------------------------------the infromation about case-control samples------------

#proporation of G in case-control sample
pg=table(unlist(ccdata))/(2*n)#which is same as the population
#the test use lrm

#proporation of Y in case-control samples
py=unlist(lapply(ccdata,length))/(2*n)

# number of categories
J=length(ccdata)

#---------------------------------test the function in mle---------------------------------
para=c(theta,beta)
proln(para=para)
gproln(para=c(2,5,3))
lamda=rep(py[J],2)
mpara=c(para,lamda)
rtrln(mpara=mpara)
grtrln(mpara = c(2,5,3,lamda))

##------------------------------MLE for prospecitive likelihhod-----------------
library(maxLik)
re=maxLik(proln,gproln,start = c(1,2,3))
summary(re)
gproln(re$estimate)


##-----------------MLE for modifiled likelihood-----------------

#maximum all parameter
A=matrix(0,3,5)
for (i in 1:(J-2)){
  A[i,i]=-1
  A[i,i+1]=1
}
A[2:3,4:5]=diag(2)
rem=maxLik(rtrln,start=c(1,2,3,0.5,2),constraints = list(ineqA=A,ineqB=rep(0,3)))
rem=maxLik(rtrln,grtrln,start=c(1,2,3,0.4,2))
#two-stage
oprla=function(lamda){rtrln(c(1,3,log(1.2),lamda))}
##lamda>0
lamda=maxLik(oprla,start = c(0.5,2),constraints = list(ineqA=diag(J-1),ineqB=rep(0,J-1)))
opr=function(para){rtrln(c(para,lamda$estimate[1:2]))}
#theta1<theta2
Amat<-matrix(0,J-2,J)
for (i in 1:(J-2)){
  Amat[i,i]=-1
  Amat[i,i+1]=1
}
Bmat=rep(0,J-2)
opra=maxLik(opr,start=c(1,2,1),constraints = list(ineqA=Amat,ineqB=Bmat))
## Iteration
Amat<-matrix(0,J-2,J)
for (i in 1:(J-2)){
  Amat[i,i]=-1
  Amat[i,i+1]=1
}
Bmat=rep(0,J-2)
e=0.01
t=1
theta_new=c(1,2,2)
theta_old=c(1,3,2)
lamda_new=c(0.05,2)
lamda_old=c(1,2)
while((sum(abs(theta_new-theta_old))>e)&&(sum(abs(lamda_new-lamda_old))>e)){
  theta_old=theta_new
  lamda_old=lamda_new
  opla=function(lamda){rtrln(c(theta_old,lamda))}
  lamda_new=maxLik(opla,start = lamda_old,constraints = list(ineqA=diag(J-1),ineqB=rep(0,J-1)))$estimate
  optheta=function(para){rtrln(c(para,lamda_new[1:2]))}
  theta_new=maxLik(optheta,start=theta_old,constraints = list(ineqA=Amat,ineqB=Bmat))$estimate
  t=t+1
}
#give the population message

r=table(Popdata$y)
la=r[J]/N*py/(r/N)
popr=function(para){rtrln(c(para,la[1:2]))}
gpopr=function(para){grtrln(c(para,la[1:2]))[1:3]}
rep1=maxLik(popr,start=c(10,20,3))
rep2=maxLik(popr,gpopr,start=c(1,2,3))
rtrln(c(2,3,log(1.8),la[1:2]))
grtrln(rem$estimate)



#---------------------------------the test use lrm in rms---------------------------------
require(rms)
x=NULL
y=NULL
ny=unlist(lapply(ccdata,length))
cny<-c(0,cumsum(ny))
J=length(ccdata)
for (j in 1:J){
  for (i in 1:ny[j]){
    x[i+cny[j]]=ccdata[[j]][i]
    y[i+cny[j]]=j
  }
}
y<-as.factor(y)
dc<-lrm(y~x)
#case-control samples
dc$coefficients
dt=lrm(y~g,data)
#the population
lrm(y~g,ccdata)
dt$coefficients#will have consistent estimate

#---------------------------------simulation---------------------------------
library(maxLik)

p=0.4
N=100000

theta=c(3.58,4.6)
beta=log(1.6)

n=500
# py=c(0.5,0.33,0.17)
# ny=n*py

Popdata=gen_data(theta,beta,p,N)
J=3
#the true lamda
r=table(Popdata$y)
Amat<-matrix(0,J-2,J)
for (i in 1:(J-2)){
  Amat[i,i]=-1
  Amat[i,i+1]=1
}
Bmat=rep(0,J-2)
A=matrix(0,3,5)
for (i in 1:(J-2)){
  A[i,i]=-1
  A[i,i+1]=1
}
A[2:3,4:5]=diag(2)
require(rms)
T=1000
proMLE=matrix(,T,2)
proMLE1=matrix(,T,2)
modMLE1=matrix(,T,2)
modMLE2=matrix(,T,2)
modMLE3=matrix(,T,2)
ptm=proc.time()
for (k in (1:T)){
  ccdata=choicedata(n,Popdata)
##proMLE  
  # proMLE[k,1]=d$coefficients[J]
  # proMLE[k,2]=d$var[J,J]
  
  re=maxLik(proln,gproln,start = c(1,2,3))
  proMLE[k,1]=re$estimate[J]
  proMLE[k,2]=solve(-re$hessian)[J,J]
  #
  
  x=NULL
  y=NULL
  ny=unlist(lapply(ccdata,length))
  cny<-c(0,cumsum(ny))
  J=length(ccdata)
  for (j in 1:J){
    for (i in 1:ny[j]){
      x[i+cny[j]]=ccdata[[j]][i]
      y[i+cny[j]]=j
    }
  }
  dc<-lrm(y~x)
  proMLE1[k,1]=dc$coefficients[J]
  proMLE1[k,2]=dc$var[J,J]
  # proMLE[k,2]=d$var[J,J]
#modMLE
  #two-stage
  # oprla=function(lamda){rtrln(c(1,3,log(1.2),lamda))}
  # ##lamda>0
  # lamda=maxLik(oprla,start =c(1:(J-1)),constraints = list(ineqA=diag(J-1),ineqB=rep(0,J-1)))
  # opr1=function(para){rtrln(c(para,lamda$estimate[1:2]))}
  # #theta1<theta2
  # Amat<-matrix(0,J-2,J)
  # for (i in 1:(J-2)){
  #   Amat[i,i]=-1
  #   Amat[i,i+1]=1
  # }
  # Bmat=rep(0,J-2)
  # op1=maxLik(opr1,start=c(1:(J-1),1),constraints = list(ineqA=Amat,ineqB=Bmat))
  # modMLE1[k,1]=op1$estimate[J]
  # modMLE1[k,2]=solve(-op1$hessian)[J,J]
 #Know the population
  #the true lamda
  # py=unlist(lapply(ccdata,length))/(2*n)
  # la=r[J]/N*py/(r/N)
  # opr2=function(para){rtrln(c(para,la[1:2]))}
  # op2=maxLik(opr2,start=c(1:(J-1),1),constraints = list(ineqA=Amat,ineqB=Bmat))
  # modMLE2[k,1]=op2$estimate[J]
  # modMLE2[k,2]=solve(-op2$hessian)[J,J]
  # ###

  op3=maxLik(rtrln,start=c(1,2,3,1,2),constraints = list(ineqA=A,ineqB=rep(0,3)))
  modMLE3[k,1]=op3$estimate[J]
  modMLE3[k,2]=solve(-op3$hessian)[J,J]
  ##
  # e=0.01
  # t=1
  # theta_new=c(1,2,2)
  # theta_old=c(1,3,2)
  # lamda_new=c(0.05,2)
  # lamda_old=c(1,2)
  # while((sum(abs(theta_new-theta_old))>e)&&(sum(abs(lamda_new-lamda_old))>e)){
  #   theta_old=theta_new
  #   lamda_old=lamda_new
  #   opla=function(lamda){rtrln(c(theta_old,lamda))}
  #   lamda_new=maxLik(opla,start = lamda_old,constraints = list(ineqA=diag(J-1),ineqB=rep(0,J-1)))$estimate
  #   optheta=function(para){rtrln(c(para,lamda_new[1:2]))}
  #   theta_new=maxLik(optheta,start=theta_old,constraints = list(ineqA=Amat,ineqB=Bmat))$estimate
  #   t=t+1
  # }
  # modMLE3[k,1]=theta_new[J]
}
proc.time()-ptm

boxplot(proMLE[,1],proMLE1[,1],modMLE3[,1], at = 1:3,col=c("2","3","5"),boxwex = 0.5,xlim=c(0,13),xaxt = "n")
boxplot(proMLE[,1],modMLE3[,1],modMLE2[,1], at = 1:3+4 ,col=c("2","3","5"),add=T,boxwex = 0.5,xaxt = "n")
boxplot(proMLE[,1],modMLE3[,1],modMLE2[,1], at = 1:3+8,col=c("2","3","5"),add=T,boxwex = 0.5,xaxt = "n")
axis(1, at = seq(2,10,4), labels = c("1","2","3"), tick = TRUE)
legend("topright", c("pro","mo1","mo2"), fill = c("2", "3","5"))
abline(h=beta,lty = 2)
pmse=mean((proMLE[,1]-beta)^2)
mmse=mean((modMLE3[,1]-beta)^2)
var(proMLE[,1])
var(modMLE3[,1])
boxplot(proMLE[,2],modMLE3[,2],ylim=c(-0.2,0.2))

