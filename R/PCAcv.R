#' @title Principal Component Analysis cross-validation error
#'
#' @description PRESS values for PCA as implemented by Eigenvector and described by Bro et al. (2008).
#'
#' @param X \code{matrix} object to perform PCA on.
#' @param ncomp \code{integer} number of components.
#'
#' @details For each number of components predicted residual sum of squares are calculated
#' based on leave-one-out cross-validation. The implementation ensures no over-fitting or
#' information bleeding.
#'
#' @return A vector of PRESS-values.
#'
#' @author Kristian Hovde Liland
#'
#' @references R. Bro, K. Kjeldahl, A.K. Smilde, H.A.L. Kiers, Cross-validation of component models: A critical look at current methods. Anal Bioanal Chem (2008) 390: 1241-1251.
#'
#' @seealso \code{\link{plot.SMI}} (print.SMI/summary.SMI), \code{\link{RV}} (RV2/RVadj), \code{\link{r1}} (r2/r3/r4/GCD), \code{\link{allCorrelations}} (matrix correlation comparison).
#'
#' @examples
#' X1  <- scale( matrix( rnorm(100*300), 100,300), scale = FALSE)
#' PCAcv(X1,10)
#'
#' @export
PCAcv <- function(X, ncomp){

  X <- as.matrix(X)
  class(X) <- "matrix"
  N <- dim(X)

  # Center X
  X <- X - rep(colMeans(X), each = N[1])

  # Set ncomp
  if(missing(ncomp)){
    ncomp <- min(N[1]-1,N[2])
  } else {
    ncomp <- min(ncomp,min(N[1]-1,N[2]))
  }

  # Prepare storage
  Xhat <- array(0, dim = c(N[1],N[2],ncomp))

  # Cross-validation (leave-one-out)
  pb <- progress_bar$new(total = N[1], format = "  [:bar] :percent (:eta)")
  for(i in 1:N[1]){
    Xi  <- X[-i, , drop = FALSE]
    Pi  <- svds(Xi, k = ncomp,  nv = ncomp, nu = 0)$v
    Xii <- matrix(rep(X[i,], N[2]), N[2], N[2], byrow = TRUE)
    diag(Xii) <- 0

    # Magic to avoid information bleed
    PiP  <- apply(Pi^2, 1, cumsum)
    PiP1 <- t(PiP/(1-PiP)+1)
    PihP <- t(Pi*(Xii%*%Pi))
    for(j in 1:N[2]){
      PP <- PihP[,j, drop = FALSE] %*% PiP1[j,, drop = FALSE]
      PP[lower.tri(PP)] <- 0
      Xhat[i,j, ] <- colSums(PP)
    }
    pb$tick()
  }

  error <- numeric(ncomp)
  for(i in 1:ncomp){
    error[i] <- sum((X-Xhat[,,i])^2)
  }
  error
}

sumNA <- function(X, dims=1){
  if(is.null(dim(X))){
    return(sum(X, na.rm=TRUE)*length(X)/sum(is.finite(X)))
  } else {
    if(dims == 1)
      return(colSums(X, na.rm=TRUE)*nrow(X)/colSums(is.finite(X)))
    else
      return(rowSums(X, na.rm=TRUE)*ncol(X)/rowSums(is.finite(X)))
  }
}
prodNA <- function(X,Y){
  actives <- (1-is.na(X)) %*% (1-is.na(Y))
  X[is.na(X)] <- 0
  Y[is.na(Y)] <- 0
  X%*%Y * ncol(X) / actives
}

meanNA <- function(X, dims=1){
  sumX <- sumNA(X, dims)
  return(sumX / dim(X)[dims])
}

#' @title Principal Component Analysis based imputation
#'
#' @description Imputation of missing data, NA, using Principal Component Analysis with
#' iterative refitting and mean value updates. The chosen number of components
#' and convergence parameters (iterations and tolerance) influence the
#' precision of the imputation.
#'
#' @param X \code{matrix} object to perform PCA on.
#' @param ncomp \code{integer} number of components.
#' @param center \code{logical} indicating if centering (default) should be performed.
#' @param max_iter \code{integer} number of iterations of PCA if sum of squared
#' change in imputed values is above \code{tol}.
#' @param tol \code{numeric} tolerance for sum of squared cange in imputed values.
#'
#' @return Final singular value decomposition, imputed \code{X} matrix and
#' convergence metrics (sequence of sum of squared change and number of iterations).
#' @export
#'
#' @examples
#' X <- matrix(rnorm(12),3,4)
#' X[c(2,6,10)] <- NA
#' PCAimpute(X, 3)
PCAimpute <- function(X, ncomp, center=TRUE, max_iter=20, tol=10^-5){
  n_row <- nrow(X)
  NAs   <- is.na(X)
  X_imp <- X
  imp_diffs <- numeric(max_iter)

  # First impute with 0 in centred X
  mean_X  <- meanNA(X)
  if(center){
    Xc  <- X-rep(mean_X, each=n_row)
  } else {
    Xc <- X
  }
  Xc[NAs] <- 0

  # SVD
  udv <- svd(Xc, ncomp, ncomp)

  # Impute
  X_hat        <- tcrossprod(udv$u%*% diag(udv$d[1:ncomp]), udv$v) + rep(mean_X, each=n_row)
  imp_diffs[1] <- sum(X_hat[NAs]^2)
  X_imp[NAs]   <- X_hat[NAs]

  # Loop
  i <- 1
  while(i < max_iter && imp_diffs[i]>tol){
    # SVD
    if(center){
      mean_X   <- colMeans(X_imp)
      udv      <- svd(X_imp-rep(mean_X, each=n_row), ncomp, ncomp)
    } else {
      udv      <- svd(X_imp, ncomp, ncomp)
    }

    # Impute
    i <- i+1
    if(center){
      X_hat <- tcrossprod(udv$u%*% diag(udv$d[1:ncomp]), udv$v) + rep(mean_X, each=n_row)
    } else {
      X_hat <- tcrossprod(udv$u%*% diag(udv$d[1:ncomp]), udv$v)
    }
    imp_diffs[i] <- sum((X_imp[NAs] - X_hat[NAs])^2)
    X_imp[NAs]   <- X_hat[NAs]
  }
  return(list(u=udv$u, d=udv$d, v=udv$v, X_hat=X_hat, imp_diffs=imp_diffs, iter=i))
}

