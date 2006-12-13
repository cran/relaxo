\name{relaxo}
\alias{relaxo}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Relaxed Lasso (relaxo)}
\description{
   Computes all "Relaxed Lasso" solutions.
}
\usage{
relaxo(X, Y, phi = seq(0, 1, length = 4), fraclambda = NULL, max.steps = min(length(Y), 2 * ncol(X)), fast = TRUE, scale = TRUE, keep.data = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{ n x p-dimensional matrix or data frame containing the predictor variables }
   \item{Y}{ n-dimensional numerical response vector }
  \item{phi}{ Relaxation parameter in [0,1]. A value of phi=1 corresponds to the regular Lasso solutions; a value of phi=0 computes the OLS estimates on the set of variables selected by the Lasso.  }
  \item{fraclambda}{ The values of |beta|/max|beta|, for which the solutions should be computed. For a NULL value, the natural Lasso solutions are computed; if given a numeric vector with entries in [0,1], solutions are forced to be computed at these values of \code{fraclambda}.  }
  \item{max.steps}{ Maximal number of steps the LARS algorithm is run. }
  \item{fast}{ Should the estimates be computed in approx. the same time as the LARS algorithm? If fast=TRUE, minor deviations from the original Relaxed Lasso solution can occur. }
  \item{scale}{ Should predictor variables be scaled (recommended)? }
  \item{keep.data}{ Should the data be kept for later usage e.g. (when computing predicted values for the training data) ?  }
}
\details{
  It is safe to leave the option for \code{fraclambda} at its default value; the option is mainly of use internally for \code{cvrelaxo}.
}
\value{
  An object of class \code{relaxo}, for which plot and predict methods are available.
}
\references{ N. Meinshausen, "Relaxed Lasso", Computational Statistics and Data Analysis, to appear. \url{http://www.stat.berkeley.edu/~nicolai} }
\author{ Nicolai Meinshausen \email{nicolai@stat.berkeley.edu}}

\seealso{ See also  \code{\link{cvrelaxo}} for computation of the cross-validated solution with optimal predictive performance }



\examples{
	data(diabetes)
	attach(diabetes)
	
## Compute "Relaxed Lasso" solution and plot results
 	object <- relaxo(x,y)
	plot(object)   
	
## Compute cross-validated solution with optimal 
## predictive performance and print relaxation parameter phi and 
## penalty parameter lambda of the found solution
	cvobject <- cvrelaxo(x,y)
	print(cvobject$phi)
	print(cvobject$lambda)

## The found penalty lambda re-expressed as fraclambda: 
## the l1-norm of the found solution divided by
## the maximal l1-norm of all solutions
        print(cvobject$fraclambda)


## Compute fitted values and plot them versus actual values     
	fitted.values <- predict(cvobject)
	plot(fitted.values,y)
	abline(c(0,1))

	detach(diabetes)       
}
\keyword{ models }
