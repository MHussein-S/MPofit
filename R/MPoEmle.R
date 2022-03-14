MPoEmle<-function (data,start)
{
  library(goftest)
  #total gradient of log-likelihood function
  grlnlNEW <- function(par, obs, ...)
    {
    #gradient of log-likelihood function
    grdNEWE<- function(par,x)
      {
      a<-par[1]
      rate<-par[2]
      G<-pexp(x,rate = rate) #CDF of exponential
      g<-dexp(x,rate = rate )  #pdf of exponential
      dG<-x*exp(-rate*x)#partial derivative of G wrt lambda
      dg<--rate*x*exp(-rate*x)+exp(-rate*x)#partial derivative of g wrt lambda
      d_a<--1/a+1/a*G+(1/a)*(G/(1+G*log(a)))
      d_rate<-log(a)*dG+dg/g+log(a)/(1+G*log(a))*dG
      res<-c(d_a,d_rate)
      return(res)
      }
    -rowSums(sapply(obs, function(x) grdNEWE(par=par,x)))
    }
  # -log-likelihood function
  objfn <- function(parm,obs,...)
    {
    n<-length(obs)
    a<-parm[1]
    rate<-parm[2]
    G<-pexp(obs,rate = rate) #CDF of exponential
    g<-dexp(obs,rate = rate )  #pdf of exponential
    LL <--n*log(a)+log(a)*sum(G)+sum(log(g))+sum(log(1+G*log(a)))
    return(-LL)
    }
  npar <- length(start)
  Mat <- diag(npar)
  colnames(Mat) <- c("a","rate")
  rownames(Mat) <- paste0("constr", 1:2)
  initconstr <- Mat %*% start - c(exp(-1),0)
  if (any(initconstr < 0))
    stop("Starting values must be in the feasible region.")
  opt <- constrOptim(theta=start,f=objfn,grad = grlnlNEW, obs = data , ui=Mat,ci = c(exp(-1),0),hessian = T)
  if (is.null(names(opt$par)))
    names(opt$par) <- c("a","rate")
  #---
  hessiana = opt$hessian
  prop_sigma<-sqrt(sqrt(diag(solve(hessiana))))
  upper<-opt$par+qnorm(0.975)*prop_sigma
  lower<-opt$par+qnorm(0.025)*prop_sigma
  #---
  loglik <- -opt$value
  k<-npar
  n<-length(data)
  AIC<--2*loglik+2*k
  AICc = -2 * loglik + 2 * k + 2 * (k * (k + 1))/(n - k - 1) #goodness.fit formula
  BIC<- -2*loglik+ k*log(n)
  HQIC<- -2*loglik+2*k*log(log(n))
  res=cbind(opt$par,prop_sigma,lower,upper)
  colnames(res)=c("MLE","Std. Err", "Inf. 95% CI","Sup. 95% CI")

  res1=cbind(AIC,AICc,BIC, HQIC,opt$value)
  colnames(res1)=c("AIC","CAIC","BIC","HQIC", "-log(Likelihood)")
  rownames(res1)=c("")
  cvm<-cvm.test(data,"pMPo",par=opt$par,distr="exp")


  res2=cbind(cvm$statistic,cvm$p.value)
  colnames(res2)=c("CVM Statistic","CVM p-value")
  rownames(res2)=c("")

  ad<-ad.test(data,"pMPo",par=opt$par,distr="exp")
  res3=cbind(ad$statistic,ad$p.value)
  colnames(res3)=c("AD Statistic","AD p-value")
  rownames(res3)=c("")

  KS<-suppressWarnings(ks.test(data,"pMPo",par=opt$par,distr="exp"))
  res4=cbind(KS$statistic,KS$p.value)
  colnames(res4)=c("KS Statistic","KS p-value")
  rownames(res4)=c("")
  res5=cbind(if(opt$convergence==0){"Algorithm Converged"} else {"Algorithm Not Converged"})
  colnames(res5)=c("")
  rownames(res5)=c("")
  list("Estimates"=res,"Measures"=res1,"Cramér–von Mises Test"=res2,"Anderson-Darling Test"=res3, "Kolmogorov-Smirnov Test"=res4,"Convergence Status"=res5)
}






