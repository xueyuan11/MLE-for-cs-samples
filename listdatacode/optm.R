

##------------------------------MLE for prospecitive likelihhod-----------------

Amat<-matrix(0,1,3)
for (i in 1:1){
  Amat[i,i]=-1
  Amat[i,i+1]=1
}
initv=c(1,2,3)
repro=constrOptim(c(1,2,3),proln,gproln,ui=Amat,ci=rep(0,1))$par
names(initv)="para"
gproln(ccdata,repro[1:(J-1)],repro[J])
#mle
mleln=function(theta1=1,theta2=2,beta=3){-proln(ccdata,theta=c(unlist(theta1),unlist(theta2)),beta=beta)}
resu=mle(mleln,method ="L-BFGS-B",lower= rep(0,3))

#optim
pproln=function(data,para){-proln(data,theta=para[1:(J-1)],beta=para[-(1:(J-1))])}
gpproln=function(data,para){-gproln(data,theta=para[1:(J-1)],beta=para[-(1:(J-1))])}
optim(c(1,2,3),pproln,ccdata=ccdata)$par
constrOptim(initv,pproln,data=ccdata,NULL,,ui=Amat,ci=rep(0,3))$par
constrOptim(initv,pproln,data=ccdata,gpproln,ui=Amat,ci=rep(0,3))$par







##-----------------MLE for modifiled likelihood-----------------
opla=function(lamda){-rtrln(ccdata,mpara,lamda)}
opgt=function(mpara){-grtrln(ccdata,mpara,lamda)}

mpara=c(1,2,1)
lamda=optimise(opla,c(0,1))$minimum
# la=constrOptim(0.5,opla,NULL,ui=lamat,ci=rep(0,1))$par
# la=py/wln(ccdata,mpara,opw)$q

r=table(data$y)
lamda=r[J]/N*py/(table(data$y)/N)
opr=function(mpara){-rtrln(ccdata,mpara,lamda[1:2])}
opgr=function(mpara){-grtrln(ccdata,mpara,lamda[1:2])[1:3]}
rertr=constrOptim(c(2,3,1),opr,NULL,ui=Amat,ci=rep(0,1))$par
constrOptim(c(2,3,1),opr,opgr,ui=Amat,ci=rep(0,1))$par
grtrln(ccdata,rertr,la)


#mle
mlerln=function(theta1=1,theta2=2,beta=3,lamda1=1,lamda2=2){-rtrln(ccdata,theta=c(unlist(theta1),unlist(theta2)),beta=beta,lamda=c(unlist(lamda1),unlist(lamda2)))}
resu=mle(mlerln,method ="L-BFGS-B",lower= rep(0,5))


#maximum all parameters
opra=function(t){-rtrln(ccdata,theta=t[1:(J-1)],beta=t[J],lamda=t[-(1:J)])}
opgra=function(t){-grtrln(ccdata,mpara=t[1:J],lamda=t[-(1:J)])}
Cmat=rbind(c(-1,1,0,0,0),c(0,0,0,1,0),c(0,0,0,0,1))
initvalue=c(1,3,2,1,2)
optim(initvalue,opra)
re=constrOptim(initvalue,opra,NULL,ui=Cmat,ci=rep(0,3))$par
res=constrOptim(initvalue,opra,opgra,ui=Cmat,ci=rep(0,3))$par