phi<-function(x){return((1+exp(-x))^(-1))}

gen_csdata<-function(p,para,ny,G=c(0,1,2)){
  J=length(ny)
  beta=para[-(1:(J-1))]
  theta=c(-Inf,para[1:(J-1)],Inf)
  pG=c(p^2,2*p*(1-p),(1-p)^2)
  P<-matrix(,J,3)
  for (i in 1:3){
    for (j in 1:J){
      P[i,j]=pG[i]*(phi(theta[j+1]-beta*G[i])-phi(theta[j]-beta*G[i]))
    }
  }
  csum<-colSums(P)
  for (j in 1:J){
    P[,j]=P[,j]/csum[j]
  }
  data<-list()
  # gen_csdata()
  for (j in 1:J){
    data[[j]]=sample(0:2,ny[j],replace = TRUE,prob = P[,j])
  }
  names(data)<-c("g1","g2","g3")
  return(data)
}

proln<-function(data,para){
  J=length(data)
  beta=para[-(1:(J-1))]
  th=c(-Inf,para[1:(J-1)],Inf)
  ny<-unlist(lapply(data,length))
  ln=0
  for (j in 1:J){
    for (i in 1:ny[j]){
      ln=ln+log(phi(th[j+1]-crossprod(beta,data[[j]][i]))-phi(th[j]-crossprod(beta,data[[j]][i])))
    }
  }
  return(ln)
}

gproln<-function(data,para){
  J=length(data)
  beta=para[-(1:(J-1))]
  th=c(-Inf,para[1:(J-1)],Inf)
  ny<-unlist(lapply(data,length))
  galn=rep(0,length(beta))
  gbln=rep(0,J-1)
  for (j in 1:J){
    for (i in 1:ny[j]){
      galn[1]=galn[1]-(1-phi(th[j+1]-crossprod(beta,data[[j]][i]))-phi(th[j]-crossprod(beta,data[[j]][i])))*data[[j]][i]
    }
  }
  for (j in 1:(J-1)){
    f=b=0
    for (i in 1:ny[j]){
      p<-phi(th[j+1]-crossprod(beta,data[[j]][i]))*(1-phi(th[j+1]-crossprod(beta,data[[j]][i])))
      q<-phi(th[j+1]-crossprod(beta,data[[j]][i]))-phi(th[j]-crossprod(beta,data[[j]][i]))
      f<-f+p/q
    }
    for (i in 1:ny[j+1]){
      p<-phi(th[j+1]-crossprod(beta,data[[j+1]][i]))*(1-phi(th[j+1]-crossprod(beta,data[[j+1]][i])))
      q<-phi(th[j+2]-crossprod(beta,data[[j+1]][i]))-phi(th[j+1]-crossprod(beta,data[[j+1]][i]))
      b<-b+p/q
    }
    gbln[j]=f-b
  }
  return(c(gbln,galn))
}

proItheta<-function(data,para){
  J=length(data)
  beta=para[-(1:(J-1))]
  th=c(-Inf,para[1:(J-1)],Inf)
  ny<-unlist(lapply(data,length))
  npar<-length(para)
  l<-matrix(0,npar,npar)
  for (j in 1:J){
    for (i in 1:ny[j]){
      a<-rep(0,length(beta))
      b<-rep(0,J-1)
      a[1]=(1-phi(th[j+1]-crossprod(beta,data[[j]][i]))-phi(th[j]-crossprod(beta,data[[j]][i])))*data[[j]][i]
      if (j==1){
         b[j]= 1-phi(th[j+1]-crossprod(beta,data[[j]][i]))
      }
      else{ 
        if(j<J){
        b[j-1]= -phi(th[j]-crossprod(beta,data[[j]][i]))*(1-phi(th[j]-crossprod(beta,data[[j]][i]))) / 
          (phi(th[j+1]-crossprod(beta,data[[j]][i]))-phi(th[j]-crossprod(beta,data[[j]][i])))
      
        b[j]= phi(th[j+1]-crossprod(beta,data[[j]][i]))*(1-phi(th[j+1]-crossprod(beta,data[[j]][i]))) / 
          (phi(th[j+1]-crossprod(beta,data[[j]][i]))-phi(th[j]-crossprod(beta,data[[j]][i])))
      }
      else{
         b[j-1]=- phi(th[j]-crossprod(beta,data[[j]][i]))
      }
      }
      c<-c(b,a)
      l<-l+c%*%t(c)
    }
  }
  return(l/sum(ny))
}
rtrln=function(data,mpara){
  J=length(data)
  beta=mpara[J]
  th=c(-Inf,mpara[1:(J-1)],Inf)
  ny<-unlist(lapply(data,length))
  n=sum(ny)
  la=c(mpara[-(1:J)],ny[J]/n)
  m=table(unlist(data))
  ln1=sum(ny*log(la))
  ln2=0
  for (j in 1:J){
    for (i in 1:ny[j]){
      ln2=ln2+log(phi(th[j+1]-crossprod(beta,data[[j]][i]))-phi(th[j]-crossprod(beta,data[[j]][i])))
    }
  }
  ln3<-rep(0,3)

  for (j in 1:J){
    ln3[1]=ln3[1]+la[j]*(phi(th[j+1])-phi(th[j]))
    ln3[2]=ln3[2]+la[j]*(phi(th[j+1]-beta)-phi(th[j]-beta))
    ln3[3]=ln3[3]+la[j]*(phi(th[j+1]-2*beta)-phi(th[j]-2*beta))
}
  ln=ln1+ln2-sum(m*log(ln3))
  return(ln)
}


grtrln=function(data,mpara){
  J=length(data)
  beta=mpara[J]
  th=c(-Inf,mpara[1:(J-1)],Inf)
  ny<-unlist(lapply(data,length))
  n=sum(ny)
  la=c(mpara[-(1:J)],ny[J]/n)
  m=table(unlist(data))
  galn=0
  gbln=rep(0,J-1)
  gcln=rep(0,J-1)
  #beta
  galn1=0
  v1=v2=de0=de1=de2=nu1=nu2=0
  for (j in 1:J){
    for (i in 1:ny[j]){
      galn1=galn1-(1-phi(th[j+1]-crossprod(beta,data[[j]][i]))-phi(th[j]-crossprod(beta,data[[j]][i])))*data[[j]][i]
    }
    v0=la[j]*(phi(th[j+1])-phi(th[j]))
    de0=de0+v0
    v1=la[j]*(phi(th[j+1]-beta)-phi(th[j]-beta))
    de1=de1+v1
    nu1=nu1+m[2]*v1*(1-phi(th[j+1]-beta)-phi(th[j]-beta))
    v2=la[j]*(phi(th[j+1]-2*beta)-phi(th[j]-2*beta))
    de2=de2+v2
    nu2=nu2+2*m[3]*v2*(1-phi(th[j+1]-2*beta)-phi(th[j]-2*beta))
  }
  galn2=nu1/de1
  galn3=nu2/de2
  galn=galn1+galn2+galn3
  
  #theta

  for (j in 1:(J-1)){
    f= b=0
    for (i in 1:ny[j]){
      p<-phi(th[j+1]-crossprod(beta,data[[j]][i]))*(1-phi(th[j+1]-crossprod(beta,data[[j]][i])))
      q<-phi(th[j+1]-crossprod(beta,data[[j]][i]))-phi(th[j]-crossprod(beta,data[[j]][i]))
      f<-f+p/q
    }
    for (i in 1:ny[j+1]){
      p<-phi(th[j+1]-crossprod(beta,data[[j+1]][i]))*(1-phi(th[j+1]-crossprod(beta,data[[j+1]][i])))
      q<-phi(th[j+2]-crossprod(beta,data[[j+1]][i]))-phi(th[j+1]-crossprod(beta,data[[j+1]][i]))
      b<-b+p/q
    }
   
   gbln[j]=f-b
   -m[1]*(la[j]-la[j+1])*phi(th[j+1])*(1-phi(th[j+1]))/de0
   -m[2]*(la[j]-la[j+1])*phi(th[j+1]-beta)*(1-phi(th[j+1]-beta))/de1
   -m[3]*(la[j]-la[j+1])*phi(th[j+1]-2*beta)*(1-phi(th[j+1]-2*beta))/de2
   
   gcln[j]=ny[j]/la[j]
   -m[1]*(phi(th[j+1])-phi(th[j]))/de0
   -m[2]*(phi(th[j+1]-beta)-phi(th[j]-beta))/de1
   -m[3]*(phi(th[j+1]-2*beta)-phi(th[j]-2*beta))/de2
  }
  ##lamda
 
  return(c(gbln,galn,gcln))   
}
opproln=function(para){-proln(data,para)}
oprtrln=function(mpara){-rtrln(data,mpara)}
opgproln=function(para){-gproln(data,para)}
opgrtrln=function(mpara){-grtrln(data,mpara)}


p=0.3
# G=c(0,1,2)
# pG=c(p^2,2*p*(1-p),(1-p)^2)
para=c(3.48,4.6,1)
mpara=c(3.48,4.6,1,5,1)
n=500
py=c(0.5,0.3,0.2)
ny=n*py
data<-gen_csdata(p,para,ny)
proln(data,para)
rtrln(data,mpara)
gproln(data,para)
grtrln(data,mpara)
I=proItheta(data,para)
T=solve(I)
Amat<-matrix(0,1,3)
for (i in 1:1){
  Amat[i,i]=-1
  Amat[i,i+1]=1
}
constrOptim(c(10,15,3),opproln,opgproln,ui=Amat,ci=rep(0,1))
constrOptim(c(1,1.5,2),opproln,NULL,ui=Amat,ci=rep(0,1))$par
Bmat=cbind(Amat,0,0)
re=constrOptim(c(2,4.6,2,5,1),oprtrln,opgrtrln,ui=Bmat,ci=rep(0,1))$par
rtrln(data,re)
# mpara<-c(1,1.5,2,0.5,0.5)
require(rms)
x=NULL
y=NULL
cny<-c(0,cumsum(ny))
J=length(data)
for (j in 1:J){
  for (i in 1:ny[j]){
    x[i+cny[j]]=data[[j]][i]
    y[i+cny[j]]=j
  }
}
y<-as.factor(y)
d<-lrm(y~x)
d$var
d$coefficients
library(MASS)
fit.polr
fit.polr <- polr(y ~ x)
library(gnlm)
