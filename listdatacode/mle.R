phi<-function(x){return((1+exp(-x))^(-1))}


###prospecitive likelihood

proln<-function(para){
  J=length(ccdata)
  beta=para[-(1:(J-1))]
  th=c(-Inf,para[1:(J-1)],Inf)
  ny<-unlist(lapply(ccdata,length))
  ln=0
  a=ifelse(order(th)==1:(J+1),1,0)
  if (sum(a)==J+1){
    for (j in 1:J){
      for (i in 1:ny[j]){
        ln=ln+log(phi(th[j+1]-crossprod(beta,ccdata[[j]][i]))-phi(th[j]-crossprod(beta,ccdata[[j]][i])))
      }
    }
  }
  else{
    ln=-Inf
  }

  return(ln)
}


##first gradient of prospecitive likelihood 
gproln<-function(para){
  J=length(ccdata)
  beta=para[-(1:(J-1))]
  th=c(-Inf,para[1:(J-1)],Inf)
  ny<-unlist(lapply(ccdata,length))
  galn=rep(0,length(beta))
  gbln=rep(0,J-1)
  for (j in 1:J){
    for (i in 1:ny[j]){
      galn[1]=galn[1]-(1-phi(th[j+1]-crossprod(beta,ccdata[[j]][i]))-phi(th[j]-crossprod(beta,ccdata[[j]][i])))*ccdata[[j]][i]
    }
  }
  for (j in 1:(J-1)){
    f=b=0
    for (i in 1:ny[j]){
      p<-phi(th[j+1]-crossprod(beta,ccdata[[j]][i]))*(1-phi(th[j+1]-crossprod(beta,ccdata[[j]][i])))
      q<-phi(th[j+1]-crossprod(beta,ccdata[[j]][i]))-phi(th[j]-crossprod(beta,ccdata[[j]][i]))
      f<-f+p/q
    }
    for (i in 1:ny[j+1]){
      p<-phi(th[j+1]-crossprod(beta,ccdata[[j+1]][i]))*(1-phi(th[j+1]-crossprod(beta,ccdata[[j+1]][i])))
      q<-phi(th[j+2]-crossprod(beta,ccdata[[j+1]][i]))-phi(th[j+1]-crossprod(beta,ccdata[[j+1]][i]))
      b<-b+p/q
    }
    gbln[j]=f-b
  }
  return(c(gbln,galn))
}

##Fisher information of pro likelihood
proItheta<-function(para){
  J=length(ccdata)
  beta=para[-(1:(J-1))]
  th=c(-Inf,para[1:(J-1)],Inf)
  ny<-unlist(lapply(ccdata,length))
  npar<-J+length(beta)-1
  l<-matrix(0,npar,npar)
  for (j in 1:J){
    for (i in 1:ny[j]){
      a<-rep(0,length(beta))
      b<-rep(0,J-1)
      a[1]=(1-phi(th[j+1]-crossprod(beta,ccdata[[j]][i]))-phi(th[j]-crossprod(beta,ccdata[[j]][i])))*ccdata[[j]][i]
      if (j==1){
         b[j]= 1-phi(th[j+1]-crossprod(beta,ccdata[[j]][i]))
      }
      else{ 
        if(j<J){
        b[j-1]= -phi(th[j]-crossprod(beta,ccdata[[j]][i]))*(1-phi(th[j]-crossprod(beta,ccdata[[j]][i]))) / 
          (phi(th[j+1]-crossprod(beta,ccdata[[j]][i]))-phi(th[j]-crossprod(beta,ccdata[[j]][i])))
      
        b[j]= phi(th[j+1]-crossprod(beta,ccdata[[j]][i]))*(1-phi(th[j+1]-crossprod(beta,ccdata[[j]][i]))) / 
          (phi(th[j+1]-crossprod(beta,ccdata[[j]][i]))-phi(th[j]-crossprod(beta,ccdata[[j]][i])))
      }
      else{
         b[j-1]=- phi(th[j]-crossprod(beta,ccdata[[j]][i]))
      }
      }
      c<-c(b,a)
      l<-l+c%*%t(c)
    }
  }
  return(l/sum(ny))
}

#modified likelihood
# rtrln1=function(ccdata,mpara,la){
#   J=length(ccdata)
#   beta=mpara[J]
#   th=c(-Inf,mpara[1:(J-1)],Inf)
#   ny<-unlist(lapply(ccdata,length))
#   n=sum(ny)
#   la=c(la,ny[J]/n)
#   # la=c(mpara[-(1:J)],ny[J]/n)
#   ln=0
#   for (j in 1:J){
#     for (i in 1:ny[j]){
#       d=0
#       for (k in 1:J){
#         d=d+la[k]*(phi(th[k+1]-crossprod(beta,ccdata[[j]][i]))-phi(th[k]-crossprod(beta,ccdata[[j]][i])))
#       }
#       ln=ln+log(la[j]*(phi(th[j+1]-crossprod(beta,ccdata[[j]][i]))-phi(th[j]-crossprod(beta,ccdata[[j]][i])))/d)
#     }
#   }
#   return(ln)
# }
rtrln=function(mpara){
  J=length(ccdata)
  beta=mpara[J]
  th=c(-Inf,mpara[1:(J-1)],Inf)
  lamda=mpara[-(1:J)]
  ny=unlist(lapply(ccdata,length))
  py<-ny/sum(ny)
  la=c(lamda,py[J])
  # la=c(mpara[-(1:J)],ny[J]/n)
  m=table(unlist(ccdata))
  ln=0
    ln1=sum(ny*log(la))
    ln2=0
    a=ifelse(order(th)==1:(J+1),1,0)
    if (sum(a)==J+1){
    for (j in 1:J){
      for (i in 1:ny[j]){
        ln2=ln2+log(phi(th[j+1]-crossprod(beta,ccdata[[j]][i]))-phi(th[j]-crossprod(beta,ccdata[[j]][i])))
      }
    }
    ln3<-rep(0,3)
    for (j in 1:J){
      ln3[1]=ln3[1]+la[j]*(phi(th[j+1])-phi(th[j]))
      ln3[2]=ln3[2]+la[j]*(phi(th[j+1]-beta)-phi(th[j]-beta))
      ln3[3]=ln3[3]+la[j]*(phi(th[j+1]-2*beta)-phi(th[j]-2*beta))
  }
    ln=ln1+ln2-sum(m*log(ln3))
    }
    else{
      ln=-Inf
    }
  return(ln)
}



##the first gradient of modified likelihood
grtrln=function(mpara){
  J=length(ccdata)
  beta=mpara[J]
  th=c(-Inf,mpara[1:(J-1)],Inf)
  lamda=mpara[-(1:J)]
  ny<-unlist(lapply(ccdata,length))
  n=sum(ny)
  la=c(lamda,ny[J]/n)
  # la=c(mpara[-(1:J)],ny[J]/n)
  m=table(unlist(ccdata))
  galn=0
  gbln=rep(0,J-1)
  gcln=rep(0,J-1)
  #beta
  galn1=0
  de0=de1=de2=nu1=nu2=0
  for (k in 1:J){
    for (i in 1:ny[k]){
      galn1=galn1-(1-phi(th[k+1]-crossprod(beta,ccdata[[k]][i]))-phi(th[k]-crossprod(beta,ccdata[[k]][i])))*ccdata[[k]][i]
    }
    v0=la[k]*(phi(th[k+1])-phi(th[k]))
    de0=de0+v0
    v1=la[k]*(phi(th[k+1]-beta)-phi(th[k]-beta))
    de1=de1+v1
    nu1=nu1+v1*(1-phi(th[k+1]-beta)-phi(th[k]-beta))
    v2=la[k]*(phi(th[k+1]-2*beta)-phi(th[k]-2*beta))
    de2=de2+v2
    nu2=nu2+v2*(1-phi(th[k+1]-2*beta)-phi(th[k]-2*beta))
  }

  galn2=m[2]*nu1/de1
  galn3=2*m[3]*nu2/de2
  galn=galn1+galn2+galn3
  
  #theta
  
  for (j in 1:(J-1)){
    f= b=0
    for (i in 1:ny[j]){
      p<-phi(th[j+1]-crossprod(beta,ccdata[[j]][i]))*(1-phi(th[j+1]-crossprod(beta,ccdata[[j]][i])))
      q<-phi(th[j+1]-crossprod(beta,ccdata[[j]][i]))-phi(th[j]-crossprod(beta,ccdata[[j]][i]))
      f<-f+p/q
    }
    for (i in 1:ny[j+1]){
      p<-phi(th[j+1]-crossprod(beta,ccdata[[j+1]][i]))*(1-phi(th[j+1]-crossprod(beta,ccdata[[j+1]][i])))
      q<-phi(th[j+2]-crossprod(beta,ccdata[[j+1]][i]))-phi(th[j+1]-crossprod(beta,ccdata[[j+1]][i]))
      b<-b+p/q
    }
    
    gbln[j]=f-b-m[1]*(la[j]-la[j+1])*phi(th[j+1])*(1-phi(th[j+1]))/de0-m[2]*(la[j]-la[j+1])*phi(th[j+1]-beta)*(1-phi(th[j+1]-beta))/de1-m[3]*(la[j]-la[j+1])*phi(th[j+1]-2*beta)*(1-phi(th[j+1]-2*beta))/de2
    
    gcln[j]=ny[j]/la[j]-m[1]*(phi(th[j+1])-phi(th[j]))/de0-m[2]*(phi(th[j+1]-beta)-phi(th[j]-beta))/de1-m[3]*(phi(th[j+1]-2*beta)-phi(th[j]-2*beta))/de2
  }
  ##lamda
  
  # list(a=c(gbln,galn),b=gcln)   
  return(c(gbln,galn,gcln))   
}






# wln=function(ccdata,theta,beta,w){
#   J=length(ccdata)
#   th=c(-Inf,theta,Inf)
#   ny<-unlist(lapply(ccdata,length))
#   n=sum(ny)
#   ln=0
#   de=rep(0,J)
#   w=c(w,1-w[1]-w[2])
#   for (j in 1:J){
#     for (k in 1:J){
#       for (i in 1:ny[k]){
#         de[j]= de[j]+w[ccdata[[k]][i]+1]*(phi(th[j+1]-crossprod(beta,ccdata[[k]][i]))-phi(th[j]-crossprod(beta,ccdata[[k]][i])))
#       }
#     }
#   }
#   de=de/n
#   for (j in 1:J){
#     for (i in 1:ny[j]){
#       ln=ln+log((w[ccdata[[j]][i]+1]*(phi(th[j+1]-crossprod(beta,ccdata[[j]][i]))-phi(th[j]-crossprod(beta,ccdata[[j]][i]))))/de[j])
#     }
#   }
#   list(ln=ln,q=de)
# }
# w=c(1/3,1/3)
# wln(ccccdata,mpara = rertr,w)
# proln(ccccdata,mpara)
# Bmatla<-diag(2)
# mpara=c(2,3,1)
# # w=rep(0.5,500)
# # wln(ccdata,para,w)
# opwln=function(w){-wln(ccccdata,mpara,w)$ln}
# Bmatw=rbind(diag(2),c(-1,-1))
# 
# opw=constrOptim(c(0.2,0.3),opwln,NULL,ui=Bmatw,ci=c(0,0,-1))$par
# lamda=py/wln(ccccdata,mpara,opw)$q
# lamda




