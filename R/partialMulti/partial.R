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

# Permuterting i parvis SMI
pairWiseSMI <- function(X, ncomp, B=0){
  nMat <- length(X)
  n    <- dim(X[[1]])[1]

  # Center and "score" all matrices
  C <- list()
  for(i in 1:nMat){
    X[[i]] <- center(unclass(as.matrix(X[[i]])))
    usv    <- svd(X[[i]])
    X[[i]] <- usv$u[, 1:ncomp[i], drop=FALSE]
    C[[i]] <- tcrossprod(X[[i]])
  }

  # Pair-wise SMI
  pws  <- matrix(0.0, nMat, nMat)
  for(i in 1:(nMat-1)){
    for(j in (i+1):nMat){
      pws[i,j] <- cor(c(C[[i]]),c(C[[j]]))
    }
  }
  pws <- pws + t(pws) + diag(1, nMat, nMat)

  if(B>0){
    pwsPerm  <- array(0, c(nMat, nMat, B))
    for(b in 1:B){
      for(i in 1:nMat){
#        X[[i]] <- X[[i]][sample(n),, drop=FALSE]
#        C[[i]] <- c(tcrossprod(X[[i]]))
        ind <- sample(n)
        C[[i]] <- C[[i]][ind,ind]
      }
      for(i in 1:(nMat-1)){
        for(j in (i+1):nMat){
          pwsPerm[i,j,b] <- cor(c(C[[i]]),c(C[[j]]))
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


################
# Wines of Loire
data(wine, package="multiblock")
library(pcalg)
library(iTOP)
nmax <- unlist(lapply(wine, ncol))
ncomp <- unlist(lapply(1:5,function(i){
  if(nmax[i]>=3)
    return(which.min(PCAcv(wine[[i]], nmax[i])))
  else
    return(2)}))
pairWiseSMI(wine, ncomp)
suffStat2 <- pairWiseSMI(wine, ncomp, 1000)
pc.fit2 <- pc(suffStat=suffStat2, indepTest=rv.link.significance, labels=names(wine),
              alpha=0.3, conservative=TRUE, solve.confl=TRUE)
plot(pc.fit2, main="")

# RV
config_matrices <- compute.config.matrices(wine)
cors       <- rv.cor.matrix(config_matrices)
cors_perm  <- run.permutations(config_matrices, nperm=1000)
suffStatRV <- list(cors=cors, cors_perm=cors_perm)
pc.fitRV   <- pc(suffStat=suffStatRV, indepTest=rv.link.significance, labels=names(wine),
              alpha=0.05, conservative=TRUE, solve.confl=TRUE)
plot(pc.fitRV, main="")

##############
# iTOP data
load("~/Documents/GitHub/MatrixCorrelation/R/partialMulti/iTOP.RData")

#ncomp <- c(2,2,2,2,2,2,2)
ncomp <- c(1,10,10,2,10,10,10)
#ncomp <- lapply(1:7,function(i)which.min(PCAcv(data[[i]], 10)))
#ncomp <- unlist(ncomp)
pairWiseSMI(data, ncomp)
suffStat2 <- pairWiseSMI(data, ncomp, 100)
pc.fit2 <- pc(suffStat=suffStat2, indepTest=rv.link.significance, labels=names(data),
              alpha=0.05, conservative=TRUE, solve.confl=TRUE)
plot(pc.fit2, main="")


ncomp <- pmin(nrow(data[[1]])-1, unlist(lapply(data, ncol)))
#ncomp <- lapply(1:7,function(i)which.min(PCAcv(data[[i]], 10)))
#ncomp <- unlist(ncomp)
pairWiseSMI(data, ncomp)
suffStat2RV <- pairWiseSMI(data, ncomp, 100)
pc.fit2RV <- pc(suffStat=suffStat2RV, indepTest=rv.link.significance, labels=names(data),
              alpha=0.05, conservative=TRUE, solve.confl=TRUE)
plot(pc.fit2RV, main="")

#################
# iTOP simulated
set.seed(2)
n = 100
p = 100
x1 = matrix(rnorm(n*p), n, p)
x2 = x1 + matrix(rnorm(n*p), n, p)
x3 = x2 + matrix(rnorm(n*p), n, p)
data = list(x1=x1, x2=x2, x3=x3)
config_matrices = compute.config.matrices(data)
cors = rv.cor.matrix(config_matrices)
cors_perm = run.permutations(config_matrices, nperm=1000)

## Not run:
library(pcalg)
suffStat = list(cors=cors, cors_perm=cors_perm)
pc.fit = pc(suffStat=suffStat, indepTest=rv.link.significance, labels=names(data),
            alpha=0.05, conservative=TRUE, solve.confl=TRUE)
plot(pc.fit, main="")


ncomp <- c(1,2,1)
pairWiseSMI(data, ncomp)
suffStat2 <- pairWiseSMI(data, ncomp, 100)
pc.fit2 <- pc(suffStat=suffStat2, indepTest=rv.link.significance, labels=names(data),
              alpha=0.05, conservative=TRUE, solve.confl=TRUE)
plot(pc.fit2, main="")

##########################################
# Hvor mange komponenter bør man velge?
# - PCA-CV per matrise?
ncomp <- pmin(nrow(data[[1]])-2, unlist(lapply(data, ncol)))
RMSECV <- as.data.frame(lapply(1:3, function(x)PCAcv(data[[x]],ncomp[x])))
matplot(RMSECV, type='l', xlab="#comp", ylim=c(0,10^5))

# - Plotte SMI og signifikans?
smi <- (SMI(x1,x2))
sig <- significant(SMI(x1,x2))
image(sig)
contour(smi, add=TRUE, col="gray")

# - Kryssvalidert SMI?

# - Kryssvalidert "PCR" x1~x2 / x2~x1
mv <- R2(pcr(x1~x2, scale=TRUE, validation="LOO"))$val[1,,]
matplot(t(mv), type="l", ylab = "R2", xlab="#comp.", ylim = c(-0.1,1), panel.first=grid())
lines(colMeans(mv), col=2, lwd=3)
##########################################


#################
# Potato

# ncomp <- lapply(1:7,function(i)which.min(PCAcv(data[[i]], 10)))
data(potato, package="multiblock")
ncomp <- rep(5,7)
suffStat2 <- pairWiseSMI(data, ncomp, 100)
pc.fit2 <- pc(suffStat=suffStat2, indepTest=rv.link.significance, labels=names(data),
              alpha=0.05, conservative=TRUE, solve.confl=TRUE)
plot(pc.fit2, main="")
