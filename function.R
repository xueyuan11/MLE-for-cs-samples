phi<-function(x){return((1+exp(-x))^(-1))}


gen_popdata = function(theta,beta,p, N) {
  pG=c(p^2,2*p*(1-p),(1-p)^2)
  g = sample(0:2,N,replace = TRUE,prob = pG)
  y = numeric(N)
  py = (1 + exp(- outer(theta, -beta*g ,"+"))) ^ (-1)
  aa = runif(N)
  for(i in 1:N)
    y[i] = sum(aa[i] > py[,i])
  y = as.numeric(as.factor(y))
  data.frame(y=y,g=g)
}

cc_data=function(n,data){
  Nc<-table(data$y)
  Nv<-c(Nc[1],sum(Nc)-Nc[1])
  #choose control samples
  pcontrol=subset(data,y==1)
  control=pcontrol[sample(Nv[1],n/2),]
  #chosse case samples
  pcase=subset(data,y>1)
  case=pcase[sample(Nv[2],n/2),]
  rbind(control,case)
}

likelihood.pro<-function(eta,obs){
  J <- length(obs)/3
  beta <- eta[J]
  theta <- c(-Inf, eta[1:(J-1)], Inf)
  G.num <- t(matrix(obs, nrow=3, ncol=J))
  n <-  rowSums(G.num)
  N <- sum(obs)
  L <- 0
  for (j in 1:J)
  {
    theta.j <- theta[j+1]
    theta.j.minus.1 <- theta[j]
    temp.0 <-  (phi(theta.j) - phi(theta.j.minus.1))
    temp.1 <-  (phi(theta.j-beta) - phi(theta.j.minus.1-beta))
    temp.2 <-  (phi(theta.j-2*beta) - phi(theta.j.minus.1-2*beta))
    L <- L + sum(G.num[j,]*log(c(temp.0,temp.1,temp.2)))
  }
  L
}
likelihood.pro.grad<-function(eta,obs){
  J <- length(obs)/3
  beta <- eta[J]
  theta <- c(-Inf, eta[1:(J-1)], Inf)
  G.num <- t(matrix(obs, nrow=3, ncol=J))
  n <-  rowSums(G.num)
  N <- sum(obs)
  g<-rep(0,J)

  for (j in 1:(J-1))
  {theta.j <- theta[j+1]
   theta.j.minus.1 <- theta[j]
   theta.j.add.1<-theta[j+2]
   f0=phi(theta.j)*(1-phi(theta.j))/(phi(theta.j) - phi(theta.j.minus.1))
   f1=phi(theta.j-beta)*(1-phi(theta.j-beta))/(phi(theta.j-beta) - phi(theta.j.minus.1-beta))
   f2=phi(theta.j-2*beta)*(1-phi(theta.j-2*beta))/(phi(theta.j-2*beta) - phi(theta.j.minus.1-2*beta))
   b0=phi(theta.j)*(1-phi(theta.j))/(phi(theta.j.add.1) - phi(theta.j))
   b1=phi(theta.j-beta)*(1-phi(theta.j-beta))/(phi(theta.j.add.1-beta) - phi(theta.j-beta))
   b2=phi(theta.j-2*beta)*(1-phi(theta.j-2*beta))/(phi(theta.j.add.1-2*beta) - phi(theta.j-2*beta))
   g[j]=sum(G.num[j,]*(c(f0,f1,f2)))-sum(G.num[j+1,]*(c(b0,b1,b2)))
  }
  for (j in 1:J)
  {
    theta.j <- theta[j+1]
    theta.j.minus.1 <- theta[j]
    temp.0 <- 0*(1- phi(theta.j) - phi(theta.j.minus.1))
    temp.1 <- -(1-phi(theta.j-beta) - phi(theta.j.minus.1-beta))
    temp.2 <- -2*(1- phi(theta.j-2*beta) - phi(theta.j.minus.1-2*beta))
    g[J] <- g[J] + sum(G.num[j,]*(c(temp.0,temp.1,temp.2)))
  }
  g
}
# ccdata=cc_data(N,data)
likelihood.mod <- function(eta, obs)   
{
  
  J <- length(obs)/3
  beta <- eta[J]
  theta <- c(-Inf, eta[1:(J-1)], Inf)
  lambda0 <- eta[(J+1):(2*J-1)]
  G.num <- t(matrix(obs, nrow=3, ncol=J))
  n <-  rowSums(G.num)
  N <- sum(obs)
  lambda <- c(n[1]/N,lambda0)
  L1 <- 0
  for (j in 1:J)
  { 
    theta.j <- theta[j+1]
    theta.j.minus.1 <- theta[j]
    temp.0 <- lambda[j] * (phi(theta.j) - phi(theta.j.minus.1))
    temp.1 <- lambda[j] * (phi(theta.j-beta) - phi(theta.j.minus.1-beta))
    temp.2 <- lambda[j] * (phi(theta.j-2*beta) - phi(theta.j.minus.1-2*beta))
    L1 <- L1 + sum(G.num[j,]*log(c(temp.0,temp.1,temp.2)))
  }
  
  temp.0 <- phi(theta)[2:(J+1)] - phi(theta)[1:J]
  temp.1 <- phi(theta-beta)[2:(J+1)] - phi(theta-beta)[1:J]
  temp.2 <- phi(theta-2*beta)[2:(J+1)] - phi(theta-2*beta)[1:J]
  r <- colSums(G.num)
  L2 <- r[1]*log(sum(lambda*temp.0)) + r[2]*log(sum(lambda*temp.1)) + 
    r[3]*log(sum(lambda*temp.2))
  L1-L2    
}

likelihood.mod.grad<-function(eta,obs)
{
  J <- length(obs)/3
  beta <- eta[J]
  theta <- c(-Inf, eta[1:(J-1)], Inf)
  lambda0 <- eta[(J+1):(2*J-1)]
  G.num <- t(matrix(obs, nrow=3, ncol=J))
  r <- colSums(G.num)
  n <-  rowSums(G.num)
  N <- sum(obs)
  lambda <- c(n[1]/N,lambda0)
  g<-rep(0,2*J-1)
  ##beta
  for (j in 1:J)
  {
    theta.j <- theta[j+1]
    theta.j.minus.1 <- theta[j]
    t0 <- 0*(1- phi(theta.j) - phi(theta.j.minus.1))
    t1 <- -(1-phi(theta.j-beta) - phi(theta.j.minus.1-beta))
    t2 <- -2*(1- phi(theta.j-2*beta) - phi(theta.j.minus.1-2*beta))
    g[J] <- g[J] + sum(G.num[j,]*(c(t0,t1,t2)))
    
  }
  temp.0 <- phi(theta)[2:(J+1)] - phi(theta)[1:J]
  temp.1 <- phi(theta-beta)[2:(J+1)] - phi(theta-beta)[1:J]
  temp.2 <- phi(theta-2*beta)[2:(J+1)] - phi(theta-2*beta)[1:J]
  de0=sum(lambda*temp.0)
  de1=sum(lambda*temp.1)
  de2=sum(lambda*temp.2)
  t.nu1<-temp.1*(1-phi(theta-beta)[2:(J+1)] - phi(theta-beta)[1:J])
  t.nu2<-temp.2*(1-phi(theta-2*beta)[2:(J+1)] - phi(theta-2*beta)[1:J])
  g[J]=g[J]+r[2]*sum(lambda*t.nu1)/de1+2*r[3]*sum(lambda*t.nu2)/de2
  #theta
  for (j in 1:(J-1))
    
  { ph.jm.0<-phi(theta[j])
    ph.jm.1<-phi(theta[j]-beta)
    ph.jm.2<-phi(theta[j]-2*beta)
    ph.j.0 <- phi(theta[j+1])
    ph.j.1 <- phi(theta[j+1]-beta)
    ph.j.2 <- phi(theta[j+1]-2*beta)
    ph.ja.0<-phi(theta[j+2])
    ph.ja.1<-phi(theta[j+2]-beta)
    ph.ja.2<-phi(theta[j+2]-2*beta)
   
  
  f0=ph.j.0*(1-ph.j.0)/(ph.j.0 - ph.jm.0)
  f1=ph.j.1*(1-ph.j.1)/(ph.j.1 - ph.jm.1)
  f2=ph.j.2*(1-ph.j.2)/(ph.j.2 - ph.jm.2)
  b0=ph.j.0*(1-ph.j.0)/(ph.ja.0 - ph.j.0)
  b1=ph.j.1*(1-ph.j.1)/(ph.ja.1 - ph.j.1)
  b2=ph.j.2*(1-ph.j.2)/(ph.ja.2 - ph.j.2)
  
  
  g[j]=sum(G.num[j,]*(c(f0,f1,f2)))-sum(G.num[j+1,]*(c(b0,b1,b2)))-r[1]*(lambda[j]-lambda[j+1])* ph.j.0*(1- ph.j.0)/de0-r[2]*(lambda[j]-lambda[j+1])*ph.j.1*(1-ph.j.1)/de1-r[3]*(lambda[j]-lambda[j+1])*ph.j.2*(1-ph.j.2)/de2
 
  g[J+j]=n[j]/lambda[j]-r[1]*(ph.j.0 - ph.jm.0)/de0-r[2]*(ph.j.1 - ph.jm.1)/de1-r[3]*(ph.j.2 - ph.jm.2)/de2
   }
 g
}
