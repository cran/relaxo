


\name{plot.relaxo}
\alias{plot.relaxo}
\title{Plot of a Relaxed Lasso Object }
\description{
Plot of Relaxed Lasso solutions for various relaxation parameters
}
\usage{
\method{plot}{relaxo}(x, type = "l", lty = 1, main = NULL, xlab = "|beta|/max|beta| (phi=1)", ylab = expression("coefficients  " * beta[j]), plotphi = unique(x$phi), ...)
}
\arguments{
  \item{x}{ A object of class \code{relaxo}}
  \item{type}{ Plot type }
  \item{lty}{ Line style }
  \item{main}{ Title of plot }
  \item{xlab}{ Label of x-axis }
  \item{ylab}{ Label of y-axis }
  \item{plotphi}{ Which values of the relaxation parameter phi should be plotted? }
  \item{\dots}{ Arguments passed on to the plot function }
}
\details{
}
\value{
No return value
}
\references{ N. Meinshausen, "Relaxed Lasso", Computational Statistics and Data Analysis, to appear.  
\url{http://www.stat.berkeley.edu/~nicolai} }
\author{ Nicolai Meinshausen \email{nicolai@stat.berkeley.edu} }

\keyword{ models }


