\name{mMPo}
\alias{mMPo}
\title{Modified Power Distribution}
\description{The rth moment about the origin of the MPo distribution for any base distribution}
\value{An object of the same length as x, giving the rth moment}
\usage{mMPo(r,par,distr)}
\arguments{

  \item{r}{moment of order r}
  \item{par}{named vector or list of the distribution parameters.}
  \item{distr}{a character string specifying the base distribution (for example, "exp" for exponential base distribution).}
  }

\references{}
\author{Mohamed Hussein}
\examples{mMPo(r=2,par=c(rate=2,a=3),distr="exp")
mMPo(r=3,par=c(rate=0.2,a=1),distr="exp")
mMPo(r=1,par=c(scale=2,shape=2,a=3),distr="weibull")
}
\keyword{modified power distribution, family of distributions}
