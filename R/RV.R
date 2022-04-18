#' @aliases RV2 RVadj RVadjMaye RVadjGhaziri
#' @title RV coefficients
#'
#' @description Three different RV coefficients: RV, RV2 and adusted RV.
#'
#' @param X1 first \code{matrix} to be compared (\code{data.frames} are also accepted).
#' @param X2 second \code{matrix} to be compared (\code{data.frames} are also accepted).
#' @param center \code{logical} indicating if input matrices should be centered (default = TRUE).
#' @param impute \code{logical} indicating if missing values are expected in \code{X1} or \code{X2} (only for RV and RV2).
#' @param version Which version of RV adjusted to apply: "Maye" (default) or "Ghaziri"
#' RV adjusted is run using the \code{RVadj} function.
#'
#' @details For each of the four coefficients a single scalar is computed to describe
#' the similarity between the two input matrices.
#'
#' @return A single value measuring the similarity of two matrices.
#'
#' @author Kristian Hovde Liland, Benjamin Leutner (RV2)
#'
#' @references
#' \itemize{
#'  \item{RV:}{ Robert, P.; Escoufier, Y. (1976). "A Unifying Tool for Linear Multivariate
#'   Statistical Methods: The RV-Coefficient". Applied Statistics 25 (3): 257-265.}
#'  \item{RV2:}{ Smilde, AK; Kiers, HA; Bijlsma, S; Rubingh, CM; van Erk, MJ (2009). "Matrix correlations
#'  for high-dimensional data: the modified RV-coefficient". Bioinformatics 25(3): 401-5.}
#'  \item{Adjusted RV:}{ Maye, CD; Lorent, J; Horgan, GW. (2011). "Exploratory analysis of multiple omics
#'  datasets using the adjusted RV coefficient". Stat Appl Genet Mol Biol. 10(14).}
#'  \item{Adjusted RV:}{ El Ghaziri, A; Qannari, E.M. (2015) "Measures of association between
#'  two datasets; Application to sensory data", Food Quality and Preference 40 (A): 116-124.}
#' }
#'
#' @seealso \code{\link{SMI}}, \code{\link{r1}} (r2/r3/r4/GCD), \code{\link{Rozeboom}}, \code{\link{Coxhead}},
#' \code{\link{allCorrelations}} (matrix correlation comparison), \code{\link{PCAcv} (cross-validated PCA)}, \code{\link{PCAimpute} (PCA based imputation)}.
#'
#' @examples
#' X1  <- matrix(rnorm(100*300),100,300)
#' usv <- svd(X1)
#' X2  <- usv$u[,-3] %*% diag(usv$d[-3]) %*% t(usv$v[,-3])
#'
#' RV(X1,X2)
#' RV2(X1,X2)
#' RVadj(X1,X2)
#'
#' # Missing data
#' X1[c(1, 50, 400, 900)] <- NA
#' X2[c(10, 200, 450, 1200)] <- NA
#' RV(X1,X2, impute = TRUE)
#' RV2(X1,X2, impute = TRUE)
#'
#' @export
RV <- function(X1, X2, center = TRUE, impute = FALSE){
  X1 <- as.matrix(X1)
  X2 <- as.matrix(X2)
  if(impute){ # NA robust centring and inner product
    if(center){
      X1 <- X1 - rep(sumNA(X1), each = nrow(X1))
      X2 <- X2 - rep(sumNA(X2), each = nrow(X1))
    }
    AA <- prodNA(X1,t(X1))
    BB <- prodNA(X2,t(X2))
  } else {
    if(center){
      X1 <- X1 - rep(colMeans(X1), each = nrow(X1))
      X2 <- X2 - rep(colMeans(X2), each = nrow(X1))
    }
    AA <- tcrossprod(X1)
    BB <- tcrossprod(X2)
  }

  RV <- Trace(AA%*%BB) / (Trace(AA%*%AA)*Trace(BB%*%BB))^0.5
  RV
}

#' @rdname RV
#' @export
RV2 <- function(X1, X2, center = TRUE, impute = FALSE){
  X1 <- as.matrix(X1)
  X2 <- as.matrix(X2)
  if(impute){ # NA robust centring and inner product
    if(center){
      X1 <- X1 - rep(sumNA(X1), each = nrow(X1))
      X2 <- X2 - rep(sumNA(X2), each = nrow(X1))
    }
    AA <- prodNA(X1,t(X1))
    BB <- prodNA(X2,t(X2))
    AA0 <- AA; diag(AA0) <- 0
    BB0 <- BB; diag(BB0) <- 0
    return(Trace(AA0%*%BB0) / (sum(AA0^2)^0.5*sum(BB0^2)^0.5))
  } else {
    if(center){
      X1 <- X1 - rep(colMeans(X1), each = nrow(X1))
      X2 <- X2 - rep(colMeans(X2), each = nrow(X1))
    }
    return(RV2cpp(X1, X2))
  }
  # AA  <- tcrossprod(X1)
  # BB  <- tcrossprod(X2)
  # AA0 <- AA; diag(AA0) <- 0
  # BB0 <- BB; diag(BB0) <- 0
  #
  # RV2 <- Trace(AA0%*%BB0) / (sum(AA0^2)^0.5*sum(BB0^2)^0.5)
  # RV2
}

#' @rdname RV
#' @export
RVadjMaye <- function(X1, X2, center = TRUE){
  X1 <- as.matrix(X1)
  X2 <- as.matrix(X2)
  if(center){
    X1 <- X1 - rep(colMeans(X1), each = nrow(X1))
    X2 <- X2 - rep(colMeans(X2), each = nrow(X1))
  }
  n <- dim(X1)[1]
  p <- dim(X1)[2]
  q <- dim(X2)[2]
  pq   <- p*q; pp <- p*p; qq <- q*q
  AA   <- tcrossprod(X1)
  BB   <- tcrossprod(X2)
  sx1  <- std(X1); sx2 <- std(X2)
  msxy <- c(min(sx1), max(sx1), min(sx2), max(sx2))

  if( any(msxy > 1+10^-12) || any(msxy < 1-10^-12)){       # Not standardized X/Y
    X1s <- X1/rep(sx1, each=p); X2s <- X2/rep(sx2, each=q) # Standardize
    AAs <- tcrossprod(X1s)
    BBs <- tcrossprod(X2s)

    # Find scaling between R2 and R2adj
    xy <- Trace(AAs %*% BBs) / (pq-(n-1) / (n-2)*(pq-Trace(AAs %*% BBs) / (n-1)^2))
    xx <- Trace(AAs %*% AAs) / (pp-(n-1) / (n-2)*(pp-Trace(AAs %*% AAs) / (n-1)^2))
    yy <- Trace(BBs %*% BBs) / (qq-(n-1) / (n-2)*(qq-Trace(BBs %*% BBs) / (n-1)^2))

    # Apply scaling to non-standarized data
    RVadj <- (Trace(AA %*% BB) / xy) / (Trace(AA %*% AA) / xx*Trace(BB %*% BB) / yy)^0.5
  } else {
    RVadj <- (pq-(n-1)/(n-2)*(pq-Trace(AA %*% BB)/(n-1)^2)) /
      sqrt((pp-(n-1)/(n-2)*(pp-Trace(AA %*% AA) / (n-1)^2)) *
             (qq-(n-1)/(n-2)*(qq-trace(BB %*% BB) / (n-1)^2)))
  }
  RVadj
}


#' @rdname RV
#' @export
RVadjGhaziri <- function(X1, X2, center = TRUE){
  X1 <- as.matrix(X1)
  X2 <- as.matrix(X2)
  if(center){
    X1 <- X1 - rep(colMeans(X1), each = nrow(X1))
    X2 <- X2 - rep(colMeans(X2), each = nrow(X1))
  }
  n <- dim(X1)[1]

  AA <- tcrossprod(X1)
  BB <- tcrossprod(X2)

  rv  <- Trace(AA %*% BB) / sqrt(Trace(AA %*% AA) * Trace(BB %*% BB))
  mrvB <- sqrt(Trace(AA)^2 / Trace(AA %*% AA)) * sqrt(Trace(BB)^2 / Trace(BB %*% BB)) / (n-1)
  aRV  <- (rv - mrvB) / (1 - mrvB)
  aRV
}


#' @rdname RV
#' @export
RVadj <- function(X1, X2, version = c("Maye","Ghaziri"), center = TRUE){
  if(length(version) == 2 || version == "Maye"){
    return( RVadjMaye(X1, X2, center = center) )
  } else {
    if(version == "Ghaziri"){
      return( RVadjGhaziri(X1, X2, center = center) )
    } else {
      stop("Unsupported or misspelled version of RV adjusted. Use \"Maye\" or \"Ghaziri\"")
    }
  }
}
