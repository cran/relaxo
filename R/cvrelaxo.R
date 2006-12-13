`cvrelaxo` <-
function(X, Y, K=5, phi= seq(0,1, length=10), fraclambda=seq(0,1,length=50), max.steps= min(length(Y),2*ncol(X)), fast = TRUE, keep.data=TRUE){

  
  n <- length(Y)
  index <- sample( rep(1:K,each=ceiling(n/K)), n, replace=FALSE ) 

  loss <- rep(0,length=length(phi)*(length(fraclambda)-1))
  for (k in 1:K){
    rel <- relaxo(X[index!=k,],Y[index!=k],phi=phi,fraclambda=fraclambda,fast=fast,scale=TRUE,keep.data=FALSE)
    
    pred <-  X[index==k,]%*%t(rel$beta)
    loss <- loss + apply( sweep(pred,1,Y[index==k])^2,2, mean )
  }

  relall <- relaxo(X,Y,phi=phi,fraclambda=fraclambda,fast=fast,scale=TRUE,keep.data=keep.data)
  select <- which.min(loss)

  relall$beta <- relall$beta[select, ,drop=FALSE]
  relall$lambda <- relall$lambda[select]
  relall$fraclambda <- relall$fraclambda[select]
  relall$phi <- relall$phi[select]
  
  
  return(relall)
}

