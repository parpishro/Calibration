log.post1<-function(d,odut,tem.dup, y, lambda, gamma){
  rho=exp(lambda)/(1+exp(lambda))
  phi=-4*log(rho);
  p<-1+1/(1+exp(-gamma))
  n<-dim(d)[1]
  dup<-NULL
  for(i in 1:dim(tem.dup)[2]){
    dup<-cbind(dup, (abs((tem.dup[,i])))^(p[i]) )
  }
  d[odut]<-exp(-(dup%*%matrix(phi,ncol=1,byrow=F)))
  d<-d+t(d)
  diags = seq(0,n*(n-1),by=n) + 1:n
  d[diags] = 1
  R<-d
  U<-try(chol(R),silent=T)
  if (class(U)=="try-error")
  {
    logpost=-9e99; mu.hat=logpost; sigma.hat=logpost
  }
  else{
    F1=mat.or.vec(n,1)+1
    k=1
    S=backsolve(U,F1,transpose=T)
    S1=crossprod(S)
    G=backsolve(U,y,transpose=T)
    A=crossprod(S,G)
    ##mu.hat
    mu.hat=solve(S1,A)
    B=y-mu.hat*F1
    B1=backsolve(U,B,transpose=T)
    B2=crossprod(B1)
    sigma.hat=(1/(n-k))*B2
    logdet = sum(log(diag(U)))
    S2=log(det(crossprod(S)))
    tem1=0.5*(exp(-lambda)/(1+exp(-lambda)))^(-0.5)*(1+exp(lambda))^(-2)*exp(lambda)
    tem11<-sum(log(tem1))
    tem2<-rep((1/40), length(tem1))
    tem22<-sum(log(tem2))
    logprior<-tem11+tem22
    logpost=logprior-(0+(n-k)/2)*log(0+(n-k)*sigma.hat/2)-0.5*S2-logdet
  }

  para=list(logpost=logpost,R=R,mu=mu.hat,sigma=sigma.hat)

  return(para)

}
