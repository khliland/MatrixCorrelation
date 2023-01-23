# Center matrix
center <- function(X){
  X - rep(colMeans(X), each=nrow(X))
}

# Partial RV coefficient
partialRV <- function(X1,X2,X3){
  # Handle inputs
  X1 <- center(unclass(as.matrix(X1)))
  X2 <- center(unclass(as.matrix(X2)))
  X3 <- center(unclass(as.matrix(X3)))

  # Vectorized configuration matrices
  S1 <- c(tcrossprod(X1))
  S2 <- c(tcrossprod(X2))
  S3 <- c(tcrossprod(X3))

  # Orthogonalize S/X1 and S/X2 on S/X3
  S1 <- resid(lm(S1 ~ S3))
  S2 <- resid(lm(S2 ~ S3))

  cor(S1,S2)
}

partialRV2 <- function(X1,X2,X3){
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

  S1 <- resid(lm(S1 ~ S3))
  S2 <- resid(lm(S2 ~ S3))

  cor(S1,S2)
}

partialSMI <- function(X1,X2,X3, ncomp1,ncomp2,ncomp3){
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

  S1 <- resid(lm(S1 ~ S3))
  S2 <- resid(lm(S2 ~ S3))

  cor(S1,S2)
}

pairWiseSMI <- function(X, ncomp, B=0){
  nMat <- length(X)
  n    <- dim(X[[1]])[1]

  # Center and "score" all matrices
  C <- list()
  for(i in 1:nMat){
    X[[i]] <- center(unclass(as.matrix(X[[i]])))
    usv    <- svd(X[[i]])
    X[[i]] <- usv$u[, 1:ncomp[i], drop=FALSE]
    C[[i]] <- c(tcrossprod(X[[i]]))
  }

  # Pair-wise SMI
  pws  <- matrix(0.0, nMat, nMat)
  for(i in 1:(nMat-1)){
    for(j in (i+1):nMat){
      pws[i,j] <- cor(C[[i]],C[[j]])
    }
  }
  pws <- pws + t(pws) + diag(1, nMat, nMat)

  if(B>0){
    pwsPerm  <- array(0, c(nMat, nMat, B))
    for(b in 1:B){
      for(i in 1:nMat){
        X[[i]] <- X[[i]][sample(n),, drop=FALSE]
        C[[i]] <- c(tcrossprod(X[[i]]))
      }
      for(i in 1:(nMat-1)){
        for(j in (i+1):nMat){
          pwsPerm[i,j,b] <- cor(C[[i]],C[[j]])
        }
      }
    }
    for(b in 1:B){
      pwsPerm[,,b] <- pwsPerm[,,b] + t(pwsPerm[,,b]) + diag(1, nMat, nMat)
    }
    return(list(cors=pws, cors_perm=pwsPerm))
  } else {
    return(pws)
  }
}
# pairWiseSMI(list(M1,M2,M3), c(2,2,2))
# suffStat2 <- pairWiseSMI(list(M1,M2,M3), c(2,2,2), 4)
# library(pcalg)
# pc.fit2 <- pc(suffStat=suffStat2, indepTest=rv.link.significance, labels=c("M1","M2","M3"),
#             alpha=0.05, conservative=TRUE, solve.confl=TRUE)
# plot(pc.fit2, main="")

permutedPartialSMI <- function(X1,X2,X3, ncomp1,ncomp2,ncomp3, B){
  # Handle inputs
  X1 <- center(unclass(as.matrix(X1)))
  X2 <- center(unclass(as.matrix(X2)))
  X3 <- center(unclass(as.matrix(X3)))

  # PCAs for each input matrix
  usv1 <- svd(X1); X1 <- usv1$u[, 1:ncomp1, drop=FALSE]
  usv2 <- svd(X2); X2 <- usv2$u[, 1:ncomp2, drop=FALSE]
  usv3 <- svd(X3); X3 <- usv3$u[, 1:ncomp3, drop=FALSE]

  cors <- numeric(B)
  n <- dim(X1)[1]

  S1 <- c(tcrossprod(X1))

  for(i in 1:B){
    X2 <- X2[sample(n),]
    X3 <- X3[sample(n),]
    # Vectorized configuration matrices
    S2 <- c(tcrossprod(X2))
    S3 <- c(tcrossprod(X3))

    S1r <- resid(lm(S1 ~ S3))
    S2r <- resid(lm(S2 ~ S3))

    cors[i] <- cor(S1r,S2r)
  }
  cors
}

M1 <- matrix(rnorm(20*10),20,10)
M2 <- matrix(rnorm(20*10),20,10)
M3 <- matrix(rnorm(20*10),20,10)
M4 <- cbind(M1[,6:10], M2[,6:10])
partialRV(M1,M2,M3)
partialRV(M1,M2,M4)
partialRV2(M1,M2,M3)
partialRV2(M1,M2,M4)
partialSMI(M1,M2,M3,2,2,2)
partialSMI(M1,M2,M4,2,2,2)

library(pcalg)
library(iTOP)
pairWiseSMI(list(M1,M2,M3), c(2,2,2))
suffStat2 <- pairWiseSMI(list(M1,M2,M3), c(2,2,2), 4)
pc.fit2 <- pc(suffStat=suffStat2, indepTest=rv.link.significance, labels=c("M1","M2","M3"),
              alpha=0.05, conservative=TRUE, solve.confl=TRUE)
plot(pc.fit2, main="")


#ncomp <- lapply(1:7,function(i)which.min(PCAcv(data[[i]], 5)))
data(potato, package="multiblock")
ncomp <- rep(5,7)
suffStat2 <- pairWiseSMI(data, ncomp, 100)
pc.fit2 <- pc(suffStat=suffStat2, indepTest=rv.link.significance, labels=names(data),
              alpha=0.05, conservative=TRUE, solve.confl=TRUE)
plot(pc.fit2, main="")
