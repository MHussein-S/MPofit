# Generate random deviates from the new distribution for any base distribution

rMPo<-function(n,par,distr)
  {
  U<-runif(n)
  qMPo(U,par,distr,lower=0,upper=100000)
  }
