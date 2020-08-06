#' Coxhead's coefficient
#'
#' @param X1 first \code{matrix} to be compared (\code{data.frames} are also accepted).
#' @param X2 second \code{matrix} to be compared (\code{data.frames} are also accepted).
#' @param weighting \code{string} indicating if weighting should be \code{sqrt(p*q)} or \code{min(p,q)} (default = 'sqrt').
#'
#' @return A single value measuring the similarity of two matrices. For diagnostic purposes it is accompanied by an attribute "canonical.correlation".
#' @export
#'
#' @references Coxhead, P; 1974. "Measuring the releationship between two sets of variables". British Journal of Mathematical and Statistical Psychology 27: 205-212.
#'
#' @importFrom stats cancor
#'
#' @seealso \code{\link{SMI}}, \code{\link{RV}} (RV2/RVadj), \code{\link{Rozeboom}}, \code{\link{r1}} (r2/r3/r4/GCD).
#'
#' @examples
#' X <- matrix(rnorm(100*13),nrow=100)
#' X1 <- X[, 1:5]  # Random normal
#' X2 <- X[, 6:12] # Random normal
#' X2[,1] <- X2[,1] + X[,5] # Overlap in one variable
#' Coxhead(X1, X2)
Coxhead <- function(X1, X2, weighting = c('sqrt','min')){
  X1 <- as.matrix(X1)
  X2 <- as.matrix(X2)
  lambda <- cancor(X1, X2)$cor # Canonical correlations
  if(weighting[1] == 'sqrt'){
    s <- sqrt(ncol(X1)*ncol(X2))
  } else {
    s <- min(ncol(X1), ncol(X2))
  }
  C <- 1 - s/sum(1/(1-lambda))

  attr(C, 'canonical.correlations') <- lambda
  C

  # X1 <- as.matrix(X1)
  # X2 <- as.matrix(X2)
  # if(center){
  #   X1 <- X1 - rep(colMeans(X1), each = nrow(X1))
  #   X2 <- X2 - rep(colMeans(X2), each = nrow(X1))
  # }
  # if(ncol(X2)>ncol(X1)){ # Reverse order if dimensionality demands it
  #   Xt <- X2
  #   X2 <- X1
  #   X1 <- Xt
  # }
  # Rxx <- crossprod(X1)
  # Ryy <- crossprod(X2)
  # Rxy <- crossprod(X1, X2)
  # Ryx <- crossprod(X2, X1)
  # RRRRinv <- solve(Ryy - Ryx%*%solve(Rxx)%*%Rxy)
  # Trace(RRRRinv %*% Ryx%*%solve(Rxx)%*%Rxy) /
  #   Trace(RRRRinv %*% Ryy)
}
