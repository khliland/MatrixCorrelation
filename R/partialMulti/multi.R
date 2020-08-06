library(pracma)
GCA <- function(X, tol=10^-12){
  n_block <- length(X)
  n <- nrow(X[[1]])
  X <- lapply(X, function(i) scale(i, scale = FALSE))
  minrank <- min(min(unlist(lapply(X,ncol))),n)
  T <- list()
  for(i in 1:n_block){
    udv <- svd(X[[i]])
    thisrank <- ifelse(sum(udv$d>tol) > 1, sum(udv$d>tol), 1)
    minrank <- min(thisrank, minrank)
    T[[i]] <- udv$u[,1:thisrank,drop=FALSE]/rep(apply(udv$u[,1:thisrank,drop=FALSE],2,sd),each=n)
  }
  
  udv <- svd(do.call(cbind,T))
  C <- udv$u[,1:minrank,drop=FALSE]
  A <- U <- list()
  for(i in 1:n_block){
    A[[i]] <- pracma::pinv(X[[i]]) %*% C
    # A[[i]] <- tcrossprod(solve(crossprod(X[[i]])), X[[i]]) %*% C
    U[[i]] <- X[[i]] %*% A[[i]]
  }
  
  ii <- 0
  R <- numeric(minrank)
  for(i in 1:(n_block-1)){
    for(j in (i+1):n_block){
      ii <- ii+1
      R <- R + diag(cor(U[[i]],U[[j]]))
    }
  }
  R <- R/ii
  list(A=A, U=U, T=T, C=C, R=R)
}

# Multiblock RV coefficient
multiRV <- function(X1,X2,X3){
  # Handle inputs
  X1 <- center(unclass(as.matrix(X1)))
  X2 <- center(unclass(as.matrix(X2)))
  X3 <- center(unclass(as.matrix(X3)))
  
  # Vectorized configuration matrices
  S1 <- c(tcrossprod(X1))
  S2 <- c(tcrossprod(X2))
  S3 <- c(tcrossprod(X3))
  
  # Generalized canonical analysis
  g <- GCA(list(as.matrix(S1), as.matrix(S2), as.matrix(S3)))
  return(g$R)
}

# Multiblock RV2 coefficient
multiRV2 <- function(X1,X2,X3){
  # Handle inputs
  X1 <- center(unclass(as.matrix(X1)))
  X2 <- center(unclass(as.matrix(X2)))
  X3 <- center(unclass(as.matrix(X3)))
  
  # Remove diagonal of a matrix  
  remDiag <- function(x){
    diag(x) <- 0
    x
  }
  
  # Vectorized configuration matrices
  S1 <- c(remDiag(tcrossprod(X1)))
  S2 <- c(remDiag(tcrossprod(X2)))
  S3 <- c(remDiag(tcrossprod(X3)))
  
  # Generalized canonical analysis
  g <- GCA(list(as.matrix(S1), as.matrix(S2), as.matrix(S3)))
  return(g$R)
}

# Multiblock RV coefficient
multiSMI <- function(X1,X2,X3, ncomp1,ncomp2,ncomp3){
  # Handle inputs
  X1 <- center(unclass(as.matrix(X1)))
  X2 <- center(unclass(as.matrix(X2)))
  X3 <- center(unclass(as.matrix(X3)))
  
  # PCAs for each input matrix
  usv1 <- svd(X1); X1 <- usv1$u[, 1:ncomp1, drop=FALSE]
  usv2 <- svd(X2); X2 <- usv2$u[, 1:ncomp2, drop=FALSE]
  usv3 <- svd(X3); X3 <- usv3$u[, 1:ncomp3, drop=FALSE]
  
  # Vectorized configuration matrices
  S1 <- c(tcrossprod(X1))
  S2 <- c(tcrossprod(X2))
  S3 <- c(tcrossprod(X3))
  
  # Generalized canonical analysis
  g <- GCA(list(as.matrix(S1), as.matrix(S2), as.matrix(S3)))
  return(g$R)
}

multiRV(M1,M2,M3)
multiRV(M1,M2,M4)
multiRV2(M1,M2,M3)
multiRV2(M1,M2,M4)
multiSMI(M1,M2,M3, 2,2,2)
multiSMI(M1,M2,M4, 2,2,2)
