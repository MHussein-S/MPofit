\name{dMPo}
\alias{dMPo}
\title{Modified Power Distribution}
\description{Probability density function of MPo distribution for any base distribution}
\value{An object of the same length as x, giving the pdf values computed at x }
\usage{dMPo(x,par,distr, log = FALSE)}
\arguments{

  \item{x}{vector of quantiles}
  \item{par}{named vector or list of the distribution parameters.}
  \item{distr}{a character string specifying the base distribution (for example, "exp" for exponential base distribution).}
  \item{log}{logical; if TRUE, the values of log(pdf), at x, are given.}
  }
\references{}

\author{Mohamed Hussein}
\examples{dMPo(x=2,par=list(a=1,rate=3),distr="exp")
dMPo(c(2,3,4),par=c(scale=2,shape=2,a=3),distr="weibull")}

\keyword{modified power distribution, family of distributions}


