\name{plotMPo}
\alias{plotMPo}
\title{Modified Power Distribution}
\description{Plots the empirical and theoretical density, CDF, P-P plot and Q-Q plot}
\value{Graph of empirical versus theoretical cdf, pdf, P-P plot and Q-Q plot}
\usage{plotMPo(data, distr, para, histo = TRUE, breaks = "default", demp = TRUE)}
\arguments{

  \item{data}{a vector of data values for which the distribution was fitted}
  \item{distr}{a character string specifying the base distribution (for example, "exp" for exponential base distribution).}
  \item{para}{named vector or list of the distribution parameters}
  \item{histo}{A logical to plot the histogram using the hist function}
  \item{breaks}{If "default" the histogram is plotted with the function hist with its default
breaks definition. Else breaks is passed to the function hist.}
  \item{demp}{A logical to plot the empirical density on the first plot (alone or superimposed on
the histogram depending of the value of the argument histo) using the density
function.}
  }
\references{}
\author{Mohamed Hussein}
\examples{
#{
data<-rMPo(n=500,par=c(rate=2,a=2),distr="exp")
plotMPo(data,distr = "exp",para=c(rate=2,a=3))
plotMPo(data,distr = "exp",para=c(rate=2,a=3),histo=FALSE)
#}
}
\keyword{modified power distribution, family of distributions}
