`predict.relaxo` <-
function(object, newX=NULL, fraclambda=NULL, phi=NULL,...) {

  if(is.null(fraclambda) | is.null(phi)){
    fraclambda <- object$fraclambda[1]
    phi <- object$phi[1]
  }
    
  ## fraclambda has to be to between 0 and 1
  stopifnot(0 <= fraclambda, fraclambda <= 1)
  if(is.null(newX)) newX <- object$X
  if(ncol(newX) != ncol(object$beta) ) stop(" number of predictors must match ")

  ## scale and center
  newX <- sweep(newX,2, object$shiftX)
  newX <- sweep(newX,2, object$scaleX, FUN="/")

  unlam <- unique(object$lambda)
  unphi <- unique(object$phi)
  lambda <- fraclambda * max(unlam)

  matchlambda <- unlam[ which.min( abs( unlam - lambda ) ) ]
  matchphi    <- unphi[ which.min( abs( unphi - phi ) ) ]

  Y <- c( newX %*% object$beta[object$lambda == matchlambda &
                               object$phi == matchphi , ] )
  Y <- Y + object$shiftY

  return(Y)
}

