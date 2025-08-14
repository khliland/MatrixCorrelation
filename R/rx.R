#' @aliases r2 r3 r4
#' @title Correlational Measures for Matrices
#'
#' @description Matrix similarity as described by Ramsey et al. (1984).
#'
#' @param X1 first \code{matrix} to be compared (\code{data.frames} are also accepted).
#' @param X2 second \code{matrix} to be compared (\code{data.frames} are also accepted).
#' @param ncomp1 (GCD) number of subspace components from the first \code{matrix} (default: full subspace).
#' @param ncomp2 (GCD) number of subspace components from the second \code{matrix} (default: full subspace).
#' @param center \code{logical} indicating if input matrices should be centered (default = TRUE).
#' @param impute \code{logical} indicating if missing values are expected in \code{X1} or \code{X2}.
#' @param impute_par named \code{list} of imputation parameters in case of NAs in X1/X2.
#'
#' @details Details can be found in Ramsey's paper:
#' \describe{
#'  \item{r1:}{ inner product correlation}
#'  \item{r2:}{ orientation-independent inner product correlation}
#'  \item{r3:}{ spectra-independent inner product correlations (including orientation)}
#'  \item{r4:}{ Spectra-Independent inner product Correlations}
#'  \item{GCD:}{ Yanai's Generalized Coefficient of Determination (GCD) Measure. To reproduce the original GCD, use all components. When \code{X1} and \code{X2} are dummy variables, GCD is proportional with Pillai's criterion: tr(W^-1(B+W)).}
#' }
#'
#' @return A single value measuring the similarity of two matrices.
#'
#' @author Kristian Hovde Liland
#'
#' @references Ramsay, JO; Berg, JT; Styan, GPH; 1984. "Matrix Correlation". Psychometrica 49(3): 403-423.
#'
#' @seealso \code{\link{SMI}}, \code{\link{RV}} (RV2/RVadj), \code{\link{Rozeboom}}, \code{\link{Coxhead}},
#' \code{\link{allCorrelations}} (matrix correlation comparison), \code{\link{PCAcv} (cross-validated PCA)}, \code{\link{PCAimpute} (PCA based imputation)}.
#'
#' @examples
#' X1  <- matrix(rnorm(100*300),100,300)
#' usv <- svd(X1)
#' X2  <- usv$u[,-3] %*% diag(usv$d[-3]) %*% t(usv$v[,-3])
#'
#' r1(X1,X2)
#' r2(X1,X2)
#' r3(X1,X2)
#' r4(X1,X2)
#' GCD(X1,X2)
#' GCD(X1,X2, 5,5)
#'
#' # Missing data
#' X1[c(1, 50, 400, 900)] <- NA
#' X2[c(10, 200, 450, 1200)] <- NA
#' r1(X1,X2, impute = TRUE)
#' r2(X1,X2, impute = TRUE)
#' r3(X1,X2, impute = TRUE)
#' r4(X1,X2, impute = TRUE)
#' GCD(X1,X2, impute = TRUE)
#' GCD(X1,X2, 5,5, impute = TRUE)
#'
#'
#' @export
r1 <- function(X1, X2, center = TRUE, impute = FALSE){
  X1 <- as.matrix(X1)
  X2 <- as.matrix(X2)
  if(impute){ # NA robust centring and inner product
    if(center){
      X1 <- X1 - rep(sumNA(X1), each = nrow(X1))
      X2 <- X2 - rep(sumNA(X2), each = nrow(X1))
    }
    return(Trace(prodNA(X1,t(X2)))/sqrt(Trace(prodNA(X1,t(X1)))*Trace(prodNA(X2,t(X2)))))
  } else {
    if(center){
    X1 <- X1 - rep(colMeans(X1), each = nrow(X1))
    X2 <- X2 - rep(colMeans(X2), each = nrow(X1))
    }
    return(Trace(crossprod(X1,X2))/sqrt(Trace(crossprod(X1))*Trace(crossprod(X2))))
  }
}

#' @rdname r1
#' @export
r2 <- function(X1, X2, center = TRUE, impute = FALSE,
               impute_par = list(max_iter=20, tol=10^-5)){
  X1 <- as.matrix(X1)
  X2 <- as.matrix(X2)
  if(impute){ # NA robust centring and inner product
    if(center){
      X1 <- X1 - rep(sumNA(X1), each = nrow(X1))
      X2 <- X2 - rep(sumNA(X2), each = nrow(X1))
    }
    udv_x1 <- PCAimpute(X1, ncomp = min(dim(X1)), center = center,
                        max_iter = impute_par$max_iter, tol = impute_par$tol)
    udv_x2 <- PCAimpute(X2, ncomp = min(dim(X2)), center = center,
                        max_iter = impute_par$max_iter, tol = impute_par$tol)
  } else {
    if(center){
      X1 <- X1 - rep(colMeans(X1), each = nrow(X1))
      X2 <- X2 - rep(colMeans(X2), each = nrow(X1))
    }
    udv_x1 <- svd(X1)
    udv_x2 <- svd(X2)
  }
  return(r1(udv_x1$u %*% diag(udv_x1$d),udv_x2$u %*% diag(udv_x2$d)))
}

#' @rdname r1
#' @export
r3 <- function(X1, X2, center = TRUE, impute = FALSE,
               impute_par = list(max_iter=20, tol=10^-5)){
  X1 <- as.matrix(X1)
  X2 <- as.matrix(X2)
  if(impute){ # NA robust centring and inner product
    if(center){
      X1 <- X1 - rep(sumNA(X1), each = nrow(X1))
      X2 <- X2 - rep(sumNA(X2), each = nrow(X1))
    }
    udv_x1 <- PCAimpute(X1, ncomp = min(dim(X1)), center = center,
                        max_iter = impute_par$max_iter, tol = impute_par$tol)
    udv_x2 <- PCAimpute(X2, ncomp = min(dim(X2)), center = center,
                        max_iter = impute_par$max_iter, tol = impute_par$tol)
  } else {
    if(center){
      X1 <- X1 - rep(colMeans(X1), each = nrow(X1))
      X2 <- X2 - rep(colMeans(X2), each = nrow(X1))
    }
    udv_x1 <- svd(X1)
    udv_x2 <- svd(X2)
  }
  r1(tcrossprod(udv_x1$u, udv_x1$v),tcrossprod(udv_x2$u, udv_x2$v))
}

#' @rdname r1
#' @export
r4 <- function(X1, X2, center = TRUE, impute = FALSE,
               impute_par = list(max_iter=20, tol=10^-5)){
  X1 <- as.matrix(X1)
  X2 <- as.matrix(X2)
  if(impute){ # NA robust centring and inner product
    if(center){
      X1 <- X1 - rep(sumNA(X1), each = nrow(X1))
      X2 <- X2 - rep(sumNA(X2), each = nrow(X1))
    }
    udv_x1 <- PCAimpute(X1, ncomp = min(dim(X1)), center = center,
                        max_iter = impute_par$max_iter, tol = impute_par$tol)
    udv_x2 <- PCAimpute(X2, ncomp = min(dim(X2)), center = center,
                        max_iter = impute_par$max_iter, tol = impute_par$tol)
  } else {
    if(center){
      X1 <- X1 - rep(colMeans(X1), each = nrow(X1))
      X2 <- X2 - rep(colMeans(X2), each = nrow(X1))
    }
    udv_x1 <- svd(X1)
    udv_x2 <- svd(X2)
  }
  r1(udv_x1$u, udv_x2$u)
}

#' @rdname r1
#' @export
GCD <- function(X1, X2, ncomp1 = min(dim(X1)), ncomp2 = min(dim(X2)),
                center = TRUE, impute = FALSE, impute_par = list(max_iter=20, tol=10^-5)){
  X1 <- as.matrix(X1)
  X2 <- as.matrix(X2)
  if(impute){ # NA robust centring and inner product
    if(center){
      X1 <- X1 - rep(sumNA(X1), each = nrow(X1))
      X2 <- X2 - rep(sumNA(X2), each = nrow(X1))
    }
    udv_x1 <- PCAimpute(X1, ncomp = min(dim(X1)), center = center,
                        max_iter = impute_par$max_iter, tol = impute_par$tol)
    udv_x2 <- PCAimpute(X2, ncomp = min(dim(X2)), center = center,
                        max_iter = impute_par$max_iter, tol = impute_par$tol)
  } else {
    if(center){
      X1 <- X1 - rep(colMeans(X1), each = nrow(X1))
      X2 <- X2 - rep(colMeans(X2), each = nrow(X1))
    }
    udv_x1 <- svd(X1)
    udv_x2 <- svd(X2)
  }
  r1(tcrossprod(udv_x1$u[,1:ncomp1,drop=FALSE]),
     tcrossprod(udv_x2$u[,1:ncomp2,drop=FALSE]))
}
