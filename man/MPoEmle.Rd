\name{MPoEmle}
\alias{MPoEmle}
\title{Modified Power Distribution}
\description{MLE of parameters, KS test statistic and p-value, AIC, CAIC, BIC, HQIC, -log-likelihood function and convergence status for MPo distribution}
\value{Maximum likelihood estimates, standard error, 95\% CI, Cramer-von Misses statistic, Anderson Darling statistic, Kolmogorov Smirnov test statistic and p-value,  Akaike Information Criterion, Consistent Akaikes Information Criterion, Bayesian Information Criterion, Hannan-Quinn information criterion,  minimum value of the negative log-likelihood function and convergence status when the MPoE distribution is fitted to some data.}

\usage{MPoEmle(data,start)}
\arguments{
  \item{data}{a vector of data values for which the distribution is to be fitted}
  \item{start}{named vector or list of the starting values of the distribution parameters}
}

\references{}
\author{Mohamed Hussein}
\examples{
#{
data<-rMPo(n=500,par=c(rate=2,a=3),distr="exp")
MPoEmle(data,start=c(a=0.1,rate=0.1))
#}
}
\keyword{modified power distribution, family of distributions}
