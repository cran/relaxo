\name{relaxo}
\alias{relaxo}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Relaxed Lasso (relaxo)}
\description{
   Computes all "Relaxed Lasso" solutions.
}
\usage{
relaxo(X, Y, phi = seq(0, 1, length = 4), max.steps = min(2*length(Y), 2 * ncol(X)), fast = TRUE,  keep.data = TRUE, warn=TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{ n x p-dimensional matrix or data frame containing the predictor variables; columns are supposed to be scaled and centered. }
   \item{Y}{ n-dimensional numerical response vector; supposed to be centered to mean 0.  }
  \item{phi}{ Relaxation parameter in [0,1]. A value of phi=1 corresponds to the regular Lasso solutions; a value of phi=0 computes the OLS estimates on the set of variables selected by the Lasso.  }
  \item{max.steps}{ Maximal number of steps the LARS algorithm is run. }
  \item{fast}{ Should the estimates be computed in approx. the same time as the LARS algorithm? If fast=TRUE, minor deviations from the original Relaxed Lasso solution can occur. }
  \item{keep.data}{ Should the data be kept for later usage e.g. (when computing predicted values for the training data) ?  }
  \item{warn}{ If TRUE, warnings are given if the predictor variables X are not centered and scaled or if the reponse variable is not centered) ?  }
}
\value{
  An object of class \code{relaxo}, for which plot and predict methods are available.
}
\references{ N. Meinshausen, "Relaxed Lasso", Computational Statistics and Data Analysis, to appear. \url{http://www.stat.berkeley.edu/~nicolai} }
\author{ Nicolai Meinshausen \email{nicolai@stat.berkeley.edu}}

\seealso{ See also  \code{\link{cvrelaxo}} for computation of the cross-validated solution with optimal predictive performance }



\examples{
	data(diabetes)

## Center and scale variables
   	x <- scale(diabetes$x)
	y <- as.numeric(scale(diabetes$y))
	
## Compute "Relaxed Lasso" solution and plot results
 	object <- relaxo(x,y)
	plot(object)   
	
## Compute cross-validated solution with optimal 
## predictive performance and print relaxation parameter phi and 
## penalty parameter lambda of the found solution
	cvobject <- cvrelaxo(x,y)
	print(cvobject$phi)
	print(cvobject$lambda)

## Compute fitted values and plot them versus actual values     
	fitted.values <- predict(cvobject)
	plot(fitted.values,y)
	abline(c(0,1))

}
\keyword{ models }
