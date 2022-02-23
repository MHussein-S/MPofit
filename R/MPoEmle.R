#MLE of the MPo exponential distribution with constrOptim
#with CAIC AD  and cvm
MPoEmle<-function (data,start)
{
  library(goftest)
  #total gradient of log-likelihood function of MPoE
  grlnlMPo <- function(par, obs, ...)
    {
    #gradient of log-likelihood function of MPoE "individual contribution"
    grdMPoE<- function(par,x)
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
    -rowSums(sapply(obs, function(x) grdMPoE(par=par,x)))
    }
  # -log-likelihood function of MPoE distribution
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
  #variance covarince matrix
  vCmatrix <- function(par, obs)
    {
    n<-length(obs)
    a<-par[1]
    rate<-par[2]
    G<-pexp(obs,rate = rate) #CDF of exponential
    g<-dexp(obs,rate = rate )  #pdf of exponential
    dg<-(1-rate*obs)*exp(-rate*obs) #partial derivative of g wrt lambda
    d2g<--obs*(2-rate*obs)*exp(-rate*obs)
    dG<-obs*exp(-rate*obs) #partial derivative of G wrt lambda
    d2G<--obs^2*exp(-rate*obs)
    domn<-1+G*log(a)

    d2a<-n/a^2-1/a^2*sum(G)-1/a^2*sum(G^2/domn^2)-1/a^2*sum(G/domn)

    d2a_lambda<-1/a*sum(dG)+1/a*sum(dG/domn^2)

    d2lambda<-log(a)*sum(d2G)+sum(d2g/g-dg^2/g^2)+log(a)*sum(d2G/domn)-log(a)^2*sum(dG^2/domn^2)

    I<-matrix(c(-d2a,-d2a_lambda,-d2a_lambda,-d2lambda),nrow = 2,ncol = 2,byrow = T)
    solve(I)
    }
  npar <- length(start)
  Mat <- diag(npar)
  colnames(Mat) <- c("a","rate")
  rownames(Mat) <- paste0("constr", 1:2)
  initconstr <- Mat %*% start - c(exp(-1),0)
  if (any(initconstr < 0))
    stop("Starting values must be in the feasible region.")
  opt <- constrOptim(theta=start,f=objfn,grad = grlnlMPo, obs = data , ui=Mat,ci = c(exp(-1),0))
  if (is.null(names(opt$par)))
    names(opt$par) <- c("a","rate")
  #---
  A<-vCmatrix(opt$par,data)
  prop_sigma<-sqrt(diag(A))
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
  cvm<-cvm.test(data,"pMPo",par=opt$par,distr="exp")
  W<-cvm$statistic
  ad<-ad.test(data,"pMPo",par=opt$par,distr="exp")
  A<-ad$statistic
  res1=cbind(AIC,AICc,BIC, HQIC,W,A, opt$value)
  colnames(res1)=c("AIC","CAIC","BIC","HQIC","W","A", "-log(Likelihood)")
  rownames(res1)=c("")
  KS<-suppressWarnings(ks.test(data,"pMPo",par=opt$par,distr="exp"))
  res2=cbind(KS$statistic,KS$p.value)
  colnames(res2)=c("KS Statistic","KS p-value")
  rownames(res2)=c("")
  res3=cbind(if(opt$convergence==0){"Algorithm Converged"} else {"Algorithm Not Converged"})
  colnames(res3)=c("")
  rownames(res3)=c("")
  list("Estimates"=res,"Measures"=res1,"Kolmogorov-Smirnov Test"=res2,"Convergence Status"=res3)
}






