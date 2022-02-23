#MLE of the MPoF distribution with constrOptim
#with CAIC AD  and cvm
MPoFmle<-function (data,start)
{
  library(evd)
  #total gradient of log-likelihood function of MPoF
  grlnlMPoF <- function(par, obs, ...)
    {
    #gradient of log-likelihood function of MPoF "individual contribution"
    grdMPoF<- function(par,x)
      {
      a<-par[1]
      k<-par[2]
      lambda<-par[3]
      G<-pfrechet(x,shape = k, scale = lambda) #CDF of Frechet
      g<-dfrechet(x,shape = k, scale = lambda)  #pdf of Weibull
      dG_k<- (x/lambda)^(-k)*log(x/lambda)*exp(-(x/lambda)^(-k))#partial derivative of G wrt k
      dG_lambda<- -(x/lambda)^(-k)*k*exp(-(x/lambda)^(-k))/lambda #partial derivative of G wrt lambda
      dg_k<-(x/lambda)^(-k-1)*exp(-(x/lambda)^(-k))/lambda-k*(x/lambda)^(-k-1)*log(x/lambda)*exp(-(x/lambda)^(-k))/lambda+k*(x/lambda)^(-k-1)*(x/lambda)^(-k)*log(x/lambda)*exp(-(x/lambda)^(-k))/lambda#partial derivative of g wrt k
      dg_lambda<- -k*(x/lambda)^(-k-1)*exp(-(x/lambda)^(-k))/lambda^2-k*(x/lambda)^(-k-1)*(-k-1)*exp(-(x/lambda)^(-k))/lambda^2-k^2*(x/lambda)^(-k-1)*(x/lambda)^(-k)*exp(-(x/lambda)^(-k))/lambda^2#partial derivative of g wrt lambda
      d_a<--1/a+1/a*G+(1/a)*(G/(1+G*log(a)))#partial derivative of log-lik wrt a
      d_k<-log(a)*dG_k+dg_k/g+log(a)/(1+G*log(a))*dG_k #partial derivative of log-lik wrt k
      d_lambda<-log(a)*dG_lambda+dg_lambda/g+log(a)/(1+G*log(a))*dG_lambda #partial derivative of log-lik wrt lambda
      res<-c(d_a,d_k,d_lambda)
      return(res)
      }
    -rowSums(sapply(obs, function(x) grdMPoF(par=par,x)))
    }
  # -log-likelihood function of MPoF distribution
  objfn <- function(parm,obs,...)
    {
    n<-length(obs)
    a<-parm[1]
    k<-parm[2]
    lambda<-parm[3]
    G<-pfrechet(obs,shape = k, scale = lambda) #CDF of frechet
    g<-dfrechet(obs,shape = k, scale = lambda)  #pdf of frechet
    LL <--n*log(a)+log(a)*sum(G)+sum(log(g))+sum(log(1+G*log(a)))
    -LL    }
  #variance covarince matrix
  vCmatrix <- function(par, obs)
    {
    n<-length(obs)
    a<-par[1]
    k<-par[2]
    lambda<-par[3]
    G<-pfrechet(obs,shape = k, scale = lambda) #CDF of frechet
    g<-dfrechet(obs,shape = k, scale = lambda)  #pdf of frechet
    dg_k<-(obs/lambda)^(-k-1)*exp(-(obs/lambda)^(-k))/lambda-k*(obs/lambda)^(-k-1)*log(obs/lambda)*exp(-(obs/lambda)^(-k))/lambda+k*(obs/lambda)^(-k-1)*(obs/lambda)^(-k)*log(obs/lambda)*exp(-(obs/lambda)^(-k))/lambda#partial derivative of g wrt k
    dg_lambda<- -k*(obs/lambda)^(-k-1)*exp(-(obs/lambda)^(-k))/lambda^2-k*(obs/lambda)^(-k-1)*(-k-1)*exp(-(obs/lambda)^(-k))/lambda^2-k^2*(obs/lambda)^(-k-1)*(obs/lambda)^(-k)*exp(-(obs/lambda)^(-k))/lambda^2#partial derivative of g wrt lambda

    d2g_k<--2*(obs/lambda)^(-k-1)*log(obs/lambda)*exp(-(obs/lambda)^(-k))/lambda+2*(obs/lambda)^(-k-1)*(obs/lambda)^(-k)*log(obs/lambda)*exp(-(obs/lambda)^(-k))/lambda+k*(obs/lambda)^(-k-1)*log(obs/lambda)^2*exp(-(obs/lambda)^(-k))/lambda-3*k*(obs/lambda)^(-k-1)*log(obs/lambda)^2*(obs/lambda)^(-k)*exp(-(obs/lambda)^(-k))/lambda+k*(obs/lambda)^(-k-1)*((obs/lambda)^(-k))^2*log(obs/lambda)^2*exp(-(obs/lambda)^(-k))/lambda
    d2g_lambda<-2*k*(obs/lambda)^(-k-1)*exp(-(obs/lambda)^(-k))/lambda^3+3*k*(obs/lambda)^(-k-1)*(-k-1)*exp(-(obs/lambda)^(-k))/lambda^3+3*k^2*(obs/lambda)^(-k-1)*(obs/lambda)^(-k)*exp(-(obs/lambda)^(-k))/lambda^3+k*(obs/lambda)^(-k-1)*(-k-1)^2*exp(-(obs/lambda)^(-k))/lambda^3+2*k^2*(obs/lambda)^(-k-1)*(-k-1)*(obs/lambda)^(-k)*exp(-(obs/lambda)^(-k))/lambda^3-k^3*(obs/lambda)^(-k-1)*(obs/lambda)^(-k)*exp(-(obs/lambda)^(-k))/lambda^3+k^3*(obs/lambda)^(-k-1)*((obs/lambda)^(-k))^2*exp(-(obs/lambda)^(-k))/lambda^3
    d2g_k_lambda<- -(obs/lambda)^(-k-1)*exp(-(obs/lambda)^(-k))/lambda^2-(obs/lambda)^(-k-1)*(-k-1)*exp(-(obs/lambda)^(-k))/lambda^2-2*(obs/lambda)^(-k-1)*(obs/lambda)^(-k)*k*exp(-(obs/lambda)^(-k))/lambda^2+k*(obs/lambda)^(-k-1)*log(obs/lambda)*exp(-(obs/lambda)^(-k))/lambda^2+k*(obs/lambda)^(-k-1)*(-k-1)*log(obs/lambda)*exp(-(obs/lambda)^(-k))/lambda^2+k*(obs/lambda)^(-k-1)*exp(-(obs/lambda)^(-k))/lambda^2+2*k^2*(obs/lambda)^(-k-1)*log(obs/lambda)*(obs/lambda)^(-k)*exp(-(obs/lambda)^(-k))/lambda^2-k*(obs/lambda)^(-k-1)*(obs/lambda)^(-k)*log(obs/lambda)*exp(-(obs/lambda)^(-k))/lambda^2-k*(obs/lambda)^(-k-1)*(-k-1)*(obs/lambda)^(-k)*log(obs/lambda)*exp(-(obs/lambda)^(-k))/lambda^2-k^2*(obs/lambda)^(-k-1)*((obs/lambda)^(-k))^2*log(obs/lambda)*exp(-(obs/lambda)^(-k))/lambda^2


    dG_k<- (obs/lambda)^(-k)*log(obs/lambda)*exp(-(obs/lambda)^(-k)) #partial derivative of G wrt k
    dG_lambda<- -(obs/lambda)^(-k)*k*exp(-(obs/lambda)^(-k))/lambda #partial derivative of G wrt lambda

    d2G_k<--(obs/lambda)^(-k)*log(obs/lambda)^2*exp(-(obs/lambda)^(-k))+((obs/lambda)^(-k))^2*log(obs/lambda)^2*exp(-(obs/lambda)^(-k))
    d2G_lambda<--(obs/lambda)^(-k)*k^2*exp(-(obs/lambda)^(-k))/lambda^2+(obs/lambda)^(-k)*k*exp(-(obs/lambda)^(-k))/lambda^2+((obs/lambda)^(-k))^2*k^2*exp(-(obs/lambda)^(-k))/lambda^2
    d2G_k_lambda<-(obs/lambda)^(-k)*k*log(obs/lambda)*exp(-(obs/lambda)^(-k))/lambda-(obs/lambda)^(-k)*exp(-(obs/lambda)^(-k))/lambda-((obs/lambda)^(-k))^2*log(obs/lambda)*k*exp(-(obs/lambda)^(-k))/lambda
    domn<-1+G*log(a)

    d2a<-n/a^2-1/a^2*sum(G)-1/a^2*sum(G^2/domn^2)-1/a^2*sum(G/domn)

    d2a_k<-1/a*sum(dG_k)+1/a*sum(dG_k/domn^2)

    d2a_lambda<-1/a*sum(dG_lambda)+1/a*sum(dG_lambda/domn^2)

    d2k<- log(a)*sum(d2G_k)+sum(d2g_k/g-dg_k^2/g^2)+log(a)*sum(d2G_k/domn)-log(a)^2*sum(dG_k^2/domn^2)

    d2k_lambda<-log(a)*sum(d2G_k_lambda)+sum(d2g_k_lambda/g-dg_k*dg_lambda/g^2)+log(a)*sum(d2G_k_lambda/domn)-log(a)^2*sum(dG_k*dG_lambda/domn^2)

    d2lambda<-log(a)*sum(d2G_lambda)+sum(d2g_lambda/g-dg_lambda^2/g^2)+log(a)*sum(d2G_lambda/domn)-log(a)^2*sum(dG_lambda^2/domn^2)

    I<-matrix(c(-d2a,-d2a_k,-d2a_lambda,-d2a_k,-d2k,-d2k_lambda,-d2a_lambda,-d2k_lambda,-d2lambda),nrow = 3,ncol = 3,byrow = T)
    solve(I)
    }
  npar <- length(start)
  Mat <- diag(npar)
  colnames(Mat) <- c("a","k","lambda")
  rownames(Mat) <- paste0("constr", 1:3)
  initconstr <- Mat %*% start - c(exp(-1),0,0)
  if (any(initconstr < 0))
    stop("Starting values must be in the feasible region.")
  opt <- constrOptim(theta=start,f=objfn,grad = grlnlMPoF, obs = data , ui=Mat,ci = c(exp(-1),0,0))
  if (is.null(names(opt$par)))
    names(opt$par) <- c("a","k","lambda")
  #---
  A<-vCmatrix(opt$par,data)
  prop_sigma<-sqrt(diag(A))
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
  cvm<-cvm.test(data,"pMPo",par=opt$par,distr="frechet")
  W<-cvm$statistic
  ad<-ad.test(data,"pMPo",par=opt$par,distr="frechet")
  A<-ad$statistic
  res1=cbind(AIC,AICc,BIC, HQIC,W,A, opt$value)
  colnames(res1)=c("AIC","CAIC","BIC","HQIC","W","A", "-log(Likelihood)")
  rownames(res1)=c("")
  # yemp<-ecdf(data)
  #  ytheo<-pMPo(sort(data),par=opt$par,'weibull' )
  #KS<-ks.test(yemp(data),ytheo)
  KS<-suppressWarnings(ks.test(data,"pMPo",par=opt$par,distr="frechet"))
  res2=cbind(KS$statistic,KS$p.value)
  colnames(res2)=c("KS Statistic","KS p-value")
  rownames(res2)=c("")
  res3=cbind(if(opt$convergence==0){"Algorithm Converged"} else {"Algorithm Not Converged"})
  colnames(res3)=c("")
  rownames(res3)=c("")
  list("Estimates"=res,"Measures"=res1,"Kolmogorov-Smirnov Test"=res2,"Convergence Status"=res3)
}






