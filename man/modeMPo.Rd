\name{modeMPo}
\alias{modeMPo}
\title{Modified Power Distribution}
\description{Mode of MPo distribution by direct maximization of log(g(x))}
\value{Approximation of the Mode of MPo distribution}
\usage{modeMPo(start=0,parm,lower=0,distr,upper)}
\arguments{
  \item{start}{Initial values for the parameters to be maximized over.}
  \item{parm}{named vector or list of the distribution parameters.}
  \item{lower}{the lower end point of the interval to be searched.}
  \item{distr}{a character string specifying the base distribution (for example, "exp" for exponential base distribution).}
  \item{upper}{the upper end point of the interval to be searched.}
  }
\references{}
\author{Mohamed Hussein}
\examples{modeMPo(parm=c(a=2,rate=0.2),distr="exp",upper=1000)
modeMPo(parm=c(a=3,rate=1),distr="exp",upper=1000)
}

\keyword{modified power distribution, family of distributions}


