MPoWmle<-function (data,start)
{  library(goftest)

  #total gradient of log-likelihood function of NEWE
  grlnlMPoW <- function(par, obs, ...)
    {
    #gradient of log-likelihood function of NEWE "individual contribution"
    grdMPoW<- function(par,x)
      {
      a<-par[1]
      k<-par[2]
      lambda<-par[3]
      G<-pweibull(x,shape = k, scale = lambda) #CDF of Weibull
      g<-dweibull(x,shape = k, scale = lambda)  #pdf of Weibull
      dG_k<- (x/lambda)^k*log(x/lambda)*exp(-(x/lambda)^k) #partial derivative of G wrt k
      dG_lambda<- -(x/lambda)^k*k*exp(-(x/lambda)^k)/lambda #partial derivative of G wrt lambda
      dg_k<-(x/lambda)^(k-1)*exp(-(x/lambda)^k)/lambda+k*(x/lambda)^(k-1)*log(x/lambda)*exp(-(x/lambda)^k)/lambda-k*(x/lambda)^(k-1)*(x/lambda)^k*log(x/lambda)*exp(-(x/lambda)^k)/lambda#partial derivative of g wrt k
      dg_lambda<- -k*(x/lambda)^(k-1)*exp(-(x/lambda)^k)/lambda^2-k*(x/lambda)^(k-1)*(k-1)*exp(-(x/lambda)^k)/lambda^2+k^2*(x/lambda)^(k-1)*(x/lambda)^k*exp(-(x/lambda)^k)/lambda^2#partial derivative of g wrt lambda
      d_a<--1/a+1/a*G+(1/a)*(G/(1+G*log(a)))#partial derivative of log-lik wrt a
      d_k<-log(a)*dG_k+dg_k/g+log(a)/(1+G*log(a))*dG_k #partial derivative of log-lik wrt k
      d_lambda<-log(a)*dG_lambda+dg_lambda/g+log(a)/(1+G*log(a))*dG_lambda #partial derivative of log-lik wrt lambda
      res<-c(d_a,d_k,d_lambda)
      return(res)
      }
    -rowSums(sapply(obs, function(x) grdMPoW(par=par,x)))
    }
  # -log-likelihood function
  objfn <- function(parm,obs,...)
    {
    n<-length(obs)
    a<-parm[1]
    k<-parm[2]
    lambda<-parm[3]
    G<-pweibull(obs,shape = k, scale = lambda) #CDF of Weibull
    g<-dweibull(obs,shape = k, scale = lambda)  #pdf of Weibull
    LL <--n*log(a)+log(a)*sum(G)+sum(log(g))+sum(log(1+G*log(a)))
    return(-LL)
    }
  npar <- length(start)
  Mat <- diag(npar)
  colnames(Mat) <- c("a","k","lambda")
  rownames(Mat) <- paste0("constr", 1:3)
  initconstr <- Mat %*% start - c(exp(-1),0,0)
  if (any(initconstr < 0))
    stop("Starting values must be in the feasible region.")
  opt <- constrOptim(theta=start,f=objfn,grad = grlnlMPoW, obs = data , ui=Mat,ci = c(exp(-1),0,0),hessian = T)
  if (is.null(names(opt$par)))
    names(opt$par) <- c("a","k","lambda")
  #---
  hessiana = opt$hessian
  prop_sigma<-sqrt(sqrt(diag(solve(hessiana))))
  #prop_sigma<-sqrt(diag(A))
  upper<-opt$par+qnorm(0.975)*prop_sigma
  lower<-opt$par+qnorm(0.025)*prop_sigma
  #---
  loglik <- -opt$value
  kpar<-npar
  n<-length(data)
  AIC<--2*loglik+2*kpar
  #CAIC<--2*loglik+kpar*(log(n)+1)
  AICc = -2 * loglik + 2 * kpar + 2 * (kpar * (kpar + 1))/(n - kpar - 1) #goodness.fit formula
  BIC<- -2*loglik+ kpar*log(n)
  HQIC<- -2*loglik+2*kpar*log(log(n))
  res=cbind(opt$par,prop_sigma,lower,upper)
  colnames(res)=c("MLE","Std. Err", "Inf. 95% CI","Sup. 95% CI")

  res1=cbind(AIC,AICc,BIC, HQIC,opt$value)
  colnames(res1)=c("AIC","CAIC","BIC","HQIC", "-log(Likelihood)")
  rownames(res1)=c("")

  cvm<-cvm.test(data,"pMPo",par=opt$par,distr="weibull")
  res2=cbind(cvm$statistic,cvm$p.value)
  colnames(res2)=c("CVM Statistic","CVM p-value")
  rownames(res2)=c("")

  ad<-ad.test(data,"pMPo",par=opt$par,distr="weibull")
  res3=cbind(ad$statistic,ad$p.value)
  colnames(res3)=c("AD Statistic","AD p-value")
  rownames(res3)=c("")

  KS<-suppressWarnings(ks.test(data,"pMPo",par=opt$par,distr="weibull"))
  res4=cbind(KS$statistic,KS$p.value)
  colnames(res4)=c("KS Statistic","KS p-value")
  rownames(res4)=c("")
  res5=cbind(if(opt$convergence==0){"Algorithm Converged"} else {"Algorithm Not Converged"})
  colnames(res5)=c("")
  rownames(res5)=c("")
  list("Estimates"=res,"Measures"=res1,"Cramér–von Mises Test"=res2,"Anderson-Darling Test"=res3, "Kolmogorov-Smirnov Test"=res4,"Convergence Status"=res5)
  }






