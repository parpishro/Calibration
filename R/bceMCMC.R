bceMCMC<-function(nmcmc,burn,thin,x,y,xtest1, lambda.ini, lambda.w.ini, gamma.ini, gamma.w.ini){
  dp<-ncol(x)
  j=0
  mcmc.ma.lambda<-matrix(nrow=nmcmc,ncol=dp,byrow=T)
  accept.lambda<-matrix(0,nrow=1,ncol=dp,byrow=T)
  reasonable.lambda<-matrix(0,nrow=1,ncol=dp,byrow=T)
  mcmc.ma.gamma<-matrix(nrow=nmcmc,ncol=dp,byrow=T)
  accept.gamma<-matrix(0,nrow=1,ncol=dp,byrow=T)

  lambda<-lambda.ini
  gamma<-gamma.ini

  res<-matrix(nrow=(nmcmc-burn)/thin,ncol=nrow(xtest1))
  v.term2<-matrix(nrow=(nmcmc-burn)/thin,ncol=nrow(xtest1))

  ff<-dist.R(x)
  cov.dd<-cov.r1.dis(x,xtest1)


  for(i in 1:1){

    for(k in 1:dp){
      lambda.cond<-NULL
      lambda.cond<-lambda
      temlambda<-0.95*rnorm(1, mean=lambda[k], sd=2.38*lambda.w.ini[k])+
        0.05*rnorm(1, mean=lambda[k], sd=0.1)
      lambda.cond[k]<-temlambda
      com.phi<-log.post1(d=ff$d, odut=ff$odut,tem.dup=ff$tem.dup, y, lambda.cond, gamma)$logpost-
        log.post1(d=ff$d, odut=ff$odut, tem.dup=ff$tem.dup, y, lambda, gamma)$logpost
      u<-runif(1)
      if(log(u)<com.phi){
        lambda[k]<-lambda.cond[k]
        accept.lambda[1,k]=accept.lambda[1,k]+1
      }
      reasonable.lambda[1, k]=reasonable.lambda[1,k]+1
      ###########gamma
      gamma.cond<-NULL
      gamma.cond<-gamma
      temgamma<-0.95*rnorm(1, mean=gamma[k], sd=2.38*gamma.w.ini[k])+
        0.05*rnorm(1, mean=gamma[k], sd=0.1)
      temg<-(temgamma)
      if((temg >=-20) & (temg <=20)){
        gamma.cond[k]<-temg
        com.gamma<-log.post1(d=ff$d, odut=ff$odut,tem.dup=ff$tem.dup, y,lambda, gamma.cond)$logpost-
          log.post1(d=ff$d, odut=ff$odut,tem.dup=ff$tem.dup,y,lambda, gamma)$logpost
        u<-runif(1)
        if(log(u)<com.gamma){
          gamma[k]<-gamma.cond[k]
          accept.gamma[1,k]=accept.gamma[1,k]+1
        }
      }
    }
    mcmc.ma.lambda[i,]<-lambda
    mcmc.ma.gamma[i,]<-gamma
  }


  for(i in 2:2){
    for(k in 1:dp){
      lambda.cond<-NULL
      lambda.cond<-lambda
      temlambda<-0.95*rnorm(1, mean=lambda[k], sd=2.38*lambda.w.ini[k])+
        0.05*rnorm(1, mean=lambda[k], sd=0.1)
      lambda.cond[k]<-temlambda
      com.phi<-log.post1(d=ff$d, odut=ff$odut,tem.dup=ff$tem.dup, y,lambda.cond, gamma)$logpost-
        log.post1(d=ff$d, odut=ff$odut,tem.dup=ff$tem.dup, y,lambda, gamma)$logpost
      u<-runif(1)
      if(log(u)<com.phi){
        lambda[k]<-lambda.cond[k]
        accept.lambda[1,k]=accept.lambda[1,k]+1
      }
      reasonable.lambda[1, k]=reasonable.lambda[1,k]+1
      ###########gamma
      gamma.cond<-NULL
      gamma.cond<-gamma
      temgamma<-0.95*rnorm(1, mean=gamma[k], sd=2.38*gamma.w.ini[k])+
        0.05*rnorm(1, mean=gamma[k], sd=0.1)
      temg<-(temgamma)
      if((temg >=-20) & (temg <=20)){
        gamma.cond[k]<-temg
        com.gamma<-log.post1(d=ff$d, odut=ff$odut,tem.dup=ff$tem.dup, y,lambda, gamma.cond)$logpost-
          log.post1(d=ff$d, odut=ff$odut,tem.dup=ff$tem.dup,y,lambda, gamma)$logpost
        u<-runif(1)
        if(log(u)<com.gamma){
          gamma[k]<-gamma.cond[k]
          accept.gamma[1,k]=accept.gamma[1,k]+1
        }
      }
    }
    mcmc.ma.lambda[i,]<-lambda
    mcmc.ma.gamma[i,]<-gamma
  }

  for(i in 3:nmcmc){
    ###### lambda
    for(k in 1:dp){
      lambda.cond<-NULL
      lambda.cond<-lambda
      tg<-mcmc.ma.lambda[1:(i-1), k]
      temlambda<-0.95*rnorm(1, mean=lambda[k], sd=2.38*sd(tg))+
        0.05*rnorm(1, mean=lambda[k], sd=0.1)
      lambda.cond[k]<-temlambda
      com.phi<-log.post1(d=ff$d, odut=ff$odut,tem.dup=ff$tem.dup, y,lambda.cond, gamma)$logpost-
        log.post1(d=ff$d, odut=ff$odut,tem.dup=ff$tem.dup,y,lambda, gamma)$logpost
      u<-runif(1)
      if(log(u)<com.phi){
        lambda[k]<-lambda.cond[k]
        accept.lambda[1,k]=accept.lambda[1,k]+1
      }
      reasonable.lambda[1, k]=reasonable.lambda[1,k]+1
      #########gamma
      gamma.cond<-NULL
      gamma.cond<-gamma
      tgg<-mcmc.ma.gamma[1:(i-1), k]
      temgamma<-0.95*rnorm(1, mean=gamma[k], sd=2.38*sd(tgg))+
        0.05*rnorm(1, mean=gamma[k], sd=0.1)
      temg<-(temgamma)
      if((temg >=-20) & (temg <=20)){
        gamma.cond[k]<-temg
        com.gamma<-log.post1(d=ff$d, odut=ff$odut,tem.dup=ff$tem.dup, y,lambda, gamma.cond)$logpost-
          log.post1(d=ff$d, odut=ff$odut,tem.dup=ff$tem.dup,y,lambda, gamma)$logpost
        u<-runif(1)
        if(log(u)<com.gamma){
          gamma[k]<-gamma.cond[k]
          accept.gamma[1,k]=accept.gamma[1,k]+1
        }
      }
    }
    mcmc.ma.lambda[i,]<-lambda
    mcmc.ma.gamma[i,]<-gamma
    if(i>burn&&((i-burn)%%thin==0)){
      j=j+1
      fit10<-pred1(ff$d,ff$odut,ff$tem.dup,cov.dd$dis.tem, y, lambda, gamma)
      res[j,]=fit10$res
      v.term2[j,]=fit10$v.term2
    }
    if ((i%%(0.1*nmcmc))==0){
      print(c(i/nmcmc))
    }
  }
  m<-list(reasonable.lambda=reasonable.lambda, accept.lambda=accept.lambda,
          mcmc.ma.lambda=mcmc.ma.lambda, accept.gamma=accept.gamma,
          mcmc.ma.gamma=mcmc.ma.gamma, res=res, v.term2=v.term2)
  return(m)
}
