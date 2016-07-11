score.variance.linear <- function(k, Y, X = matrix(1, length(Y)), K=NULL, ...) {
  if ( !is.null(K) & !is.matrix(K) & !is.list(K) ) stop("K must be a matrix, a list of matrix or 'NULL'") 
  if (!is.matrix(k)) stop("'k' must to be a matrix") 
  if(length(Y) != nrow(k) | length(Y) != ncol(k)  ) stop("Dimensions of Y and k mismatch")
  if (is.matrix(K)) if ( length(Y) != nrow(K) | length(Y) != ncol(K) ) stop("Dimensions of Y and K mismatch")
  if (is.list(K)) if (length( unique(unlist(lapply(K, function(x) c(nrow(x),ncol(x)) )) ))!=1 ) stop("Matrix in K must be square and of same dimensions")
  if (is.list(K)) if  (length(Y) != unique(unlist(lapply(K, function(x) c(nrow(x),ncol(x)) )) )) stop("Dimensions of Y and K mismatch")
  if(!is.matrix(X) & !is.null(X)) stop("X should be a matrix or NULL")
  
  if( any(is.na(Y)) )
  {
    w <- !is.na(Y)
    if (!is.null(X)) X <- as.matrix(X[w,])
    Y <- Y[w]
    if (is.list(K)) K <- lapply(K, function(x) x[w,w]) else if (is.matrix(K)) K <- K[w,w]
    k <- k[w,w]
    warning(sum(!w), ' individuals with missing phenotype are ignored.\n')
  }
  
  if (!is.null(X))
  {
    r <- ncol(X)
    if(length(Y) != nrow(X)) stop("Dimensions of Y and X mismatch")
    
    if ( any( rowSums(is.na(X))>0 ) )
    {
      w <- (rowSums(is.na(X))==0)
      X <- as.matrix(X[w,])
      Y <- Y[w]
      if (is.list(K)) K <- lapply(K, function(x) x[w,w]) else if (is.matrix(K)) K <- K[w,w]
      k <- k[w,w]
      warning(sum(!w), ' individuals with missing covariates are ignored.\n')
    }
  } else r <- 0 
    
  if (is.list(K) | is.matrix(K))
  {
    model <- lmm.aireml(Y, X = X, K = K, get.P = TRUE, ...)
    eig <- eigen(model$P)
    eig$values[ eig$values<0 ] <- 0
    PP <- eig$vectors%*%( sqrt(eig$values)*t(eig$vectors) )
    T <- PP%*%k%*%PP
    s <- (t(Y)%*%model$P)%*%k%*%(model$P%*%Y)
  } else {
    if (!is.null(X)) P <- diag(1, length(Y)) - X%*%solve(t(X)%*%X)%*%t(X) else P <- diag(1, length(Y))
    P <- P/( sum( (P%*%Y)**2 )/(length(Y)-r) )
    T <- P%*%k%*%P*sum( (P%*%Y)**2 )/(length(Y)-r)
    s <- (t(Y)%*%P)%*%k%*%(P%*%Y)
  }  
  l <- eigen(T, symmetric=TRUE)$values
  l <- l[l>1e-8]
  p <- CompQuadForm::davies(s, lambda = l, acc = 1e-6)$Qq
  
  return(list(score=s,p=p))
}

score.variance.logistic <- function(k, Y, X = matrix(1, length(Y)), K, ...) {
  if ( !is.null(K) & !is.matrix(K) & !is.list(K) ) stop("K must be a matrix, a list of matrix or 'NULL'") 
  if (!is.matrix(k)) stop("'k' must to be a matrix") 
  if(length(Y) != nrow(k) | length(Y) != ncol(k)  ) stop("Dimensions of Y and k mismatch")
  if (is.matrix(K)) if ( length(Y) != nrow(K) | length(Y) != ncol(K) ) stop("Dimensions of Y and K mismatch")
  if (is.list(K)) if (length( unique(unlist(lapply(K, function(x) c(nrow(x),ncol(x)) )) ))!=1 ) stop("Matrix in K must be square and of same dimensions")
  if (is.list(K)) if  (length(Y) != unique(unlist(lapply(K, function(x) c(nrow(x),ncol(x)) )) )) stop("Dimensions of Y and K mismatch")
  if(!is.matrix(X) & !is.null(X)) stop("X should be a matrix or NULL")
  
  if( any(is.na(Y)) )
  {
    w <- !is.na(Y)
    if (!is.null(X)) X <- as.matrix(X[w,])
    Y <- Y[w]
    if (is.list(K)) K <- lapply(K, function(x) x[w,w]) else if (is.matrix(K)) K <- K[w,w]
    k <- k[w,w]
    warning(sum(!w), ' individuals with missing phenotype are ignored.\n')
  }
  
  if (!is.null(X))
  {
    if(length(Y) != nrow(X)) stop("Dimensions of Y and X mismatch")
    
    if ( any( rowSums(is.na(X))>0 ) )
    {
      w <- ( rowSums(is.na(X))==0 )
      X <- as.matrix(X[w,])
      Y <- Y[w]
      if (is.list(K)) K <- lapply(K, function(x) x[w,w]) else if (is.matrix(K)) K <- K[w,w]
      k <- k[w,w]
      warning(sum(!w), ' individuals with missing covariates are ignored.\n')
    } }

  if (is.list(K) | is.matrix(K))
  {
    model <- logistic.mm.aireml(Y, X=X, K=K, get.P = TRUE, ...)
    omega <- model$BLUP_omega
    if (!is.null(X)) omega <- omega + X%*%model$BLUP_beta
    pi <- 1/(1+exp(-omega))	
 
    eig <- eigen(model$P)
    eig$values[ eig$values<0 ] <- 0
    PP <- eig$vectors%*%( sqrt(eig$values)*t(eig$vectors) )
    T <- PP%*%k%*%PP
    s <- t(Y-pi)%*%k%*%(Y-pi)
  } else {
    if (!is.null(X))
    {
      model <- glm(Y~X, family=binomial())
      pi <- model$fitted.values
      V <- diag( sqrt(pi*(1-pi)) )
      P <- V - V%*%X%*%solve(t(X)%*%(V**2)%*%X)%*%t(X)%*%V**2
    } else {
      pi <- rep(1/2, length(Y)) 
      V <- diag( sqrt(pi*(1-pi)) )
      P <- V
    }
    T <- P%*%k%*%t(P)
    s <- t(Y-pi)%*%k%*%(Y-pi)
  }
  l <- eigen(T, symmetric=TRUE)$values
  l <- l[l>1e-8]
  p <- CompQuadForm::davies(s, lambda = l, acc = 1e-6)$Qq
  
  return(list(score=s,p=p))
}


score.fixed.linear <- function(x, Y, X = matrix(1, length(Y)), K, ...) {
  if ( !is.matrix(K) & !is.list(K) ) stop("K must be a matrix or a list of matrix") 
  if (!is.matrix(x)) stop("'x' must to be a matrix") 
  if(length(Y) != nrow(x)) stop("Dimensions of Y and x mismatch")
  if (is.matrix(K)) if ( length(Y) != nrow(K) | length(Y) != ncol(K) ) stop("Dimensions of Y and K mismatch")
  if (is.list(K)) if (length( unique(unlist(lapply(K, function(x) c(nrow(x),ncol(x)) )) ))!=1 ) stop("Matrix in K must be square and of same dimensions")
  if (is.list(K)) if  (length(Y) != unique(unlist(lapply(K, function(x) c(nrow(x),ncol(x)) )) )) stop("Dimensions of Y and K mismatch")
  if(!is.matrix(X) & !is.null(X)) stop("X should be a matrix or NULL")
  
  if( any(is.na(Y)) | any( rowSums(is.na(x))>0 ) )
  {
    w <- ( !is.na(Y) & rowSums(is.na(x))==0 )
    if (!is.null(X)) X <- as.matrix(X[w,])
    Y <- Y[w]
    if (is.list(K)) K <- lapply(K, function(x) x[w,w]) else K <- K[w,w]
    x <- as.matrix(x[w,])
    warning(sum(!w), ' individuals with missing phenotype are ignored.\n')
  }
  
  if (!is.null(X))
  {
    if(length(Y) != nrow(X)) stop("Dimensions of Y and X mismatch")
    
    if ( any( rowSums(is.na(X))>0 ) )
    {
      w <- ( rowSums(is.na(X))==0 )
      X <- as.matrix(X[w,])
      Y <- Y[w]
      if (is.list(K)) K <- lapply(K, function(x) x[w,w]) else if (is.matrix(K)) K <- K[w,w]
      x <- as.matrix(x[w,])
      warning(sum(!w), ' individuals with missing covariates are ignored.\n')
   } }

  model <- lmm.aireml(Y, X = X, K = K, get.P = TRUE, ...)
  XtP <- t(x)%*%model$P
  V <- XtP%*%x 
  T <- XtP%*%Y 
  
  s <- t(T)%*%solve(V)%*%T
  logp <- pchisq( s, df = ncol(x), lower.tail=FALSE, log.p=TRUE)
  
  return(list(score=s,p=exp(logp), log.p=logp))
}


score.fixed.logistic <- function(x, Y, X = matrix(1,  length(Y)), K, ...) {
  if ( !is.matrix(K) & !is.list(K) ) stop("K must be a matrix, a list of matrix") 
  if (!is.matrix(x)) stop("'x' must to be a matrix") 
  if(length(Y) != nrow(x)) stop("Dimensions of Y and x mismatch")
  if (is.matrix(K)) if ( length(Y) != nrow(K) | length(Y) != ncol(K) ) stop("Dimensions of Y and K mismatch")
  if (is.list(K)) if (length( unique(unlist(lapply(K, function(x) c(nrow(x),ncol(x)) )) ))!=1 ) stop("Matrix in K must be square and of same dimensions")
  if (is.list(K)) if  (length(Y) != unique(unlist(lapply(K, function(x) c(nrow(x),ncol(x)) )) )) stop("Dimensions of Y and K mismatch")
  if(!is.matrix(X) & !is.null(X)) stop("X should be a matrix or NULL")
  
  if( any(is.na(Y)) | any( rowSums(is.na(x))>0 ) )
  {
    w <- ( !is.na(Y) & rowSums(is.na(x))==0 )
    if (!is.null(X)) X <- as.matrix(X[w,])
    Y <- Y[w]
    if (is.list(K)) K <- lapply(K, function(x) x[w,w]) else K <- K[w,w]
    x <- as.matrix(x[w,])
    warning(sum(!w), ' individuals with missing phenotype are ignored.\n')
  }
  
  if (!is.null(X))
  {
    if(length(Y) != nrow(X)) stop("Dimensions of Y and X mismatch")
    
    if ( any( rowSums(is.na(X))>0 ) )
    {
      w <- ( rowSums(is.na(X))==0 )
      X <- as.matrix(X[w,])
      Y <- Y[w]
      if (is.list(K)) K <- lapply(K, function(x) x[w,w]) else if (is.matrix(K)) K <- K[w,w]
      x <- as.matrix(x[w,])
      warning(sum(!w), ' individuals with missing covariates are ignored.\n')
   } }
  
  model <- logistic.mm.aireml(Y, X=X, K=K, get.P = TRUE, ...)
    
  omega <- model$BLUP_omega
  if (!is.null(X)) omega <- omega + X%*%model$BLUP_beta
  pi <- 1/(1+exp(-omega))	
  
  T <- t(x)%*%(Y-pi)
  V <- t(x)%*%model$P%*%x
  
  s <- t(T)%*%solve(V)%*%T
  logp <- pchisq( s, df = ncol(x), lower.tail=FALSE, log.p=TRUE)
  
  return(list(score=s,p=exp(logp), log.p=logp))
}
