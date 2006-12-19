\name{cvrelaxo}
\alias{cvrelaxo}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Cross validation for "Relaxed Lasso" }
\description{
	Compute the "Relaxed Lasso" solution with minimal cross-validated L2-loss.
}
\usage{
cvrelaxo(X, Y, K = 5, phi = seq(0, 1, length = 10), max.steps = min( 2* length(Y), 2 * ncol(X)), fast = TRUE, keep.data = TRUE, warn=TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{ as in function \code{relaxo} }
  \item{Y}{ as in function \code{relaxo} }
  \item{K}{ Number of folds. Defaults to 5. }
  \item{phi}{  as in function \code{relaxo} } 
  \item{max.steps}{ as in function \code{relaxo} }
  \item{fast}{ as in function \code{relaxo} }
  \item{keep.data}{ as in function \code{relaxo} }
  \item{warn}{ as in function \code{relaxo} }
}
\value{
	An object of class \code{relaxo}, for which print and predict methods exist
}

\details{
 The plot method is not useful for result of \code{cvrelaxo} (as no path of solutions exists).
}
\references{ N. Meinshausen, "Relaxed Lasso", Computational Statistics and Data Analysis, to appear. \url{http://www.stat.berkeley.edu/~nicolai} }
\author{ Nicolai Meinshausen \email{nicolai@stat.berkeley.edu}}

\seealso{ See also \code{\link{relaxo}} for computation of the entire solution path }



\examples{
	data(diabetes)

## Center and scale variables
   	x <- scale(diabetes$x)
	y <- scale(diabetes$y)
	
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

