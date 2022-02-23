\name{qMPo}
\alias{qMPo}
\title{Modified Power Exponential Distribution}
\description{Quantile function of MPo distribution for any base distribution}
\value{An object of the same length as p, giving the quantile function values computed at p}
\usage{qMPo(p,par,distr,lower=0,upper,lower.tail=TRUE,log.p=FALSE)}
\arguments{
  \item{p}{scalar or vector of probabilities.}
  \item{par}{named vector or list of the distribution parameters.}
  \item{distr}{character string specifying the base distribution (for example, "exp" for exponential base distribution).}
  \item{lower}{the lower end point of the interval to be searched.}
  \item{upper}{the upper end point of the interval to be searched.}
  \item{lower.tail}{logical; if TRUE, probability are P(X <= x), otherwise, P(X > x).}
  \item{log.p}{logical; if TRUE, probabilities p are given as log(p).}
}
\references{}
\author{Mohamed Hussein}
\examples{qMPo(p=0.5,c(rate=0.2,a=2),"exp",upper=100)
  qMPo(0.25,c(rate=0.2,a=2),"exp",upper=100)
  qMPo(0.75,c(a=2),"norm",upper=100)
  qMPo(log(0.5),c(a=5,rate=0.5),"exp",lower=0,upper=100,log.p=TRUE)
  qMPo(0.10,c(a=1.5,rate=0.5),"exp",lower=0,upper=100,lower.tail=FALSE)
  qMPo(p=0.5,c(rate=0.2,a=2),"exp",upper=100)
}
\keyword{modified power distribution, family of distributions}
