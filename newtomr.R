

indi<-function(y){ny=length(table(y))
c<-matrix(,n,ny)
for (i in 1:n){
  for (j in 1:ny){
  c[i,j]=ifelse(y[i]==j,1,0)
  }}
return(c)}

phi<-function(x){
  return((1+exp(-x))^(-1))
}
rtrln<-function(x,y,para,lamda){
  ln=0
  ny=length(table(y))
  c=indi(y)
  th<-c(-Inf,para[-1],Inf)
  for (i in 1:n){
    to=0
    for ( k in 1:ny){
      to<-to+lamda[k]*(phi(th[k+1]-crossprod(para[1],x[i]))-phi(th[k]-crossprod(para[1],x[i])))
    }
    for (j in 1:ny){
      ln=ln+c[i,j]*(lamda[j]*(phi(th[j+1]-crossprod(para[1],x[i]))-phi(th[j]-crossprod(para[1],x[i])))/to)
    }
  } 
  return(ln)
}
proln<-function(x,y,para){
  ln=0
  ny=length(table(y))
  c=indi(y)
  th<-c(-Inf,para[-1],Inf)
  for (i in 1:n){
    for (j in 1:ny){
      ln=ln+c[i,j]*(phi(th[j+1]-crossprod(para[1],x[i]))-phi(th[j]-crossprod(para[1],x[i])))
    }
  } 
  return(ln)
}
gproln<-function(x,y,para){
  gln=rep(0,length(para))
  ny=length(table(y))
  c=indi(y)
  th<-c(-Inf,para[-1],Inf)
  q=NULL
   p=NULL
  
  for (i in 1:n){
    for (j in 1:ny){
      gln[1]=gln[1]+c[i,j]*(-1)*(1-phi(th[j+1]-crossprod(para[1],x[i]))-phi(th[j]-crossprod(para[1],x[i])))*x[i]
    }
    q[1]<-c[i,1]*(phi(th[2]-crossprod(para[1],x[i]))-phi(th[1]-crossprod(para[1],x[i])))^(-1)
    for (j in 1:(ny-1)){
      p[j]<-phi(phi(th[j+1]-crossprod(para[1],x[i])))*(1-phi(th[j+1]-crossprod(para[1],x[i])))
      q[j+1]<-c[i,j+1]*(phi(th[j+2]-crossprod(para[1],x[i]))-phi(th[j+1]-crossprod(para[1],x[i])))^(-1)
      gln[j+1]<-gln[j+1]+p[j]*(q[j]-q[j+1])
    }
  } 
  gln
}



x=ccdata$z
y=ccdata$y
opproln<-function(par){
  -proln(x,y,par)
}
opgproln<-function(par){
  -gproln(x,y,par)
}

Amat<-matrix(0,2,4)
for (i in 1:2){
    Amat[i,i+1]=-1
    Amat[i,i+2]=1
}
constrOptim(c(0,1,1.5,3),opproln,opgproln,ui=Amat,ci=rep(0,2))
proln(x,y,para)
rtrln(x,y,para,lamda)
gproln(x,y,para)
