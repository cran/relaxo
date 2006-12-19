\name{predict.relaxo}
\alias{predict.relaxo}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ predict method for class relaxo }
\description{
	Prediction with Relaxed Lasso
}
\usage{
\method{predict}{relaxo}(object, newX = NULL, lambda = NULL, phi = NULL,...)
}
\arguments{
  \item{object}{ An object of class \code{relaxo} }
  \item{newX}{ A data frame or matrix containing new data. If not given, the prediction for the original training data is returned. }
  \item{lambda}{ The penalty parameter for variable selection. }
  \item{phi}{ The relaxation parameter in the interval [0,1].  }
  \item{...}{ other arguments (currently ignored). }
}

\value{
    A numerical vector, containing the predictions for the new data points (or the fitted values if newX=NULL).	
}

\references{ N. Meinshausen, "Relaxed Lasso", Computational Statistics and Data Analysis, to appear.  \url{http://www.stat.berkeley.edu/~nicolai}}
\author{ Nicolai Meinshausen \email{nicolai@stat.berkeley.edu} }

\keyword{ models }

