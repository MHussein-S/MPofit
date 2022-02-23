# Generate random deviates from the new distribution for any base distribution

rMPo<-function(n,par,distr)
  {
  U<-runif(n)
  x<-qMPo(U,par,distr,lower=0,upper=1000)
  }
