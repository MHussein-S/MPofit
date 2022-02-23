\name{rMPo}
\alias{rMPo}
\title{Modified Power Distribution}
\description{Generates random deviates from MPo distribution for any base distribution }
\value{Random deviates from the modified power distribution, length of the result is determined by n}
\usage{rMPo(n,par,distr)}
\arguments{

  \item{n}{number of random variates to be generated}
  \item{par}{named vector or list of the distribution parameters}
  \item{distr}{a character string specifying the base distribution (for example, "exp" for exponential base distribution).}
  }

\references{}
\author{Mohamed Hussein}
\examples{x<-rMPo(n=500,par=c(rate=2,a=3),distr="exp")
rMPo(n=100,par=c(scale=2,shape=2,a=3),distr="weibull")}
\keyword{modified power distribution, family of distributions}
