#' @encoding UTF-8
#' @title Similarity of Matrices Index (SMI)
#'
#' @description A similarity index for comparing coupled data matrices.
#'
#' @param X1 first \code{matrix} to be compared (\code{data.frames} are also accepted).
#' @param X2 second \code{matrix} to be compared (\code{data.frames} are also accepted).
#' @param ncomp1 maximum number of subspace components from the first \code{matrix}.
#' @param ncomp2 maximum number of subspace components from the second \code{matrix}.
#' @param projection type of projection to apply, defaults to "Orthogonal", alternatively "Procrustes".
#' @param Scores1 user supplied score-\code{matrix} to replace singular value decomposition of first \code{matrix}.
#' @param Scores2 user supplied score-\code{matrix} to replace singular value decomposition of second \code{matrix}.
#' @param impute \code{logical} for activation of PCA based imputation for X1/X2.
#' @param impute_par named \code{list} of imputation parameters in case of NAs in X1/X2.
#'
#' @details A two-step process starts with extraction of stable subspaces using
#' Principal Component Analysis or some other method yielding two orthonormal bases. These bases
#' are compared using Orthogonal Projection (OP / ordinary least squares) or Procrustes
#' Rotation (PR). The result is a similarity measure that can be adjusted to various
#' data sets and contexts and which includes explorative plotting and permutation based testing
#' of matrix subspace equality.
#'
#' @return A matrix containing all combinations of components. Its class is "SMI" associated
#' with print, plot, summary methods.
#'
#' @author Kristian Hovde Liland
#'
#' @references Ulf Geir Indahl, Tormod Næs, Kristian Hovde Liland; 2018. A similarity index for comparing coupled matrices. Journal of Chemometrics; e3049.
#'
#' @seealso \code{\link{plot.SMI}} (print.SMI/summary.SMI), \code{\link{RV}} (RV2/RVadj), \code{\link{r1}} (r2/r3/r4/GCD), \code{\link{Rozeboom}}, \code{\link{Coxhead}},
#' \code{\link{allCorrelations}} (matrix correlation comparison), \code{\link{PCAcv} (cross-validated PCA)}, \code{\link{PCAimpute} (PCA based imputation)}.
#'
#' @examples
#' # Simulation
#' X1  <- scale( matrix( rnorm(100*300), 100,300), scale = FALSE)
#' usv <- svd(X1)
#' X2  <- usv$u[,-3] %*% diag(usv$d[-3]) %*% t(usv$v[,-3])
#'
#' (smi <- SMI(X1,X2,5,5))
#' plot(smi, B = 1000 ) # default B = 10000
#'
#' # Sensory analysis
#' data(candy)
#' plot( SMI(candy$Panel1, candy$Panel2, 3,3, projection = "Procrustes"),
#'     frame = c(2,2), B = 1000, x1lab = "Panel1", x2lab = "Panel2" ) # default B = 10000
#'
#' # Missing data (100 missing completely at random points each)
#' X1[sort(round(runif(100)*29999+1))] <- NA
#' X2[sort(round(runif(100)*29999+1))] <- NA
#' (smi <- SMI(X1,X2,5,5, impute = TRUE))
#'
#' @importFrom pracma Rank subspace Trace std
#' @importFrom RSpectra svds
#' @importFrom plotrix color.legend
#' @importFrom Rcpp sourceCpp
#' @importFrom stats cor
#' @useDynLib MatrixCorrelation, .registration = TRUE
#' @export
SMI <- function(X1, X2, ncomp1 = Rank(X1)-1, ncomp2 = Rank(X2)-1,
                projection = "Orthogonal", Scores1 = NULL, Scores2 = NULL,
                impute = FALSE, impute_par = list(max_iter=20, tol=10^-5)){

  # Initialize
  smi <- matrix(0, ncomp1, ncomp2)

  # Check inputs
  projection <- match.arg(projection, c("Orthogonal", "Procrustes"))

  # Compute scores by PCA if not supplied
  has.scores <- TRUE
  if(is.null(Scores1) && is.null(Scores2)){
    mat.names <- c(substitute(X1), substitute(X2))
    X1  <- as.matrix(X1)
    X2  <- as.matrix(X2)
    nobj <- dim(X1)[1]
    if(dim(X2)[1] != nobj) stop('The number of objects (rows) must be equal.')
    if(dim(X1)[2] < ncomp1 || dim(X2)[2] < ncomp2) stop('The number of components cannot exceed the number of data columns.')
    if(min(dim(X1)) < 3 || min(dim(X2)) < 3 || dim(X1)[2] == ncomp1 ||  dim(X2)[2] == ncomp2){
      if(impute){
        X1imp <- PCAimpute(X1, ncomp1, max_iter=impute_par$max_iter, tol=impute_par$tol)
        X2imp <- PCAimpute(X2, ncomp2, max_iter=impute_par$max_iter, tol=impute_par$tol)
        Scores1 <- X1imp$u
        Scores2 <- X2imp$u
      } else {
        Scores1 <- svd(X1 - rep(colMeans(X1), each = nobj))$u
        Scores2 <- svd(X2 - rep(colMeans(X2), each = nobj))$u
        Scores1 <- Scores1[ ,1:ncomp1, drop = FALSE]
        Scores2 <- Scores2[ ,1:ncomp2, drop = FALSE]
      }
    } else {
      if(impute){
        X1imp <- PCAimpute(X1, ncomp1, max_iter=impute_par$max_iter, tol=impute_par$tol)
        X2imp <- PCAimpute(X2, ncomp2, max_iter=impute_par$max_iter, tol=impute_par$tol)
        Scores1 <- X1imp$u
        Scores2 <- X2imp$u
      } else {
        Scores1 <- svds(X1 - rep(colMeans(X1), each = nobj), k = ncomp1, nu = ncomp1, nv = 0)$u
        Scores2 <- svds(X2 - rep(colMeans(X2), each = nobj), k = ncomp2, nu = ncomp2, nv = 0)$u
      }
    }
    has.scores <- FALSE
  } else {
    mat.names <- c(substitute(Scores1), substitute(Scores2))
    nobj <- dim(Scores1)[1]
  }

  # Compute SMI and tests
  if(projection == "Orthogonal"){
    if(ncomp1 > 1 && ncomp2 > 1){
      smi <- apply(
        apply(
          crossprod(Scores1[, 1:ncomp1, drop=FALSE], Scores2[, 1:ncomp2, drop=FALSE])^2,
          1,cumsum),
        1,cumsum) / matrix(pmin(rep(1:ncomp1,ncomp2), rep(1:ncomp2,each=ncomp1)), ncomp1,ncomp2)
    }
    if(ncomp1 == 1 && ncomp2 == 1){
      smi <- crossprod(Scores1[, 1:ncomp1, drop=FALSE], Scores2[, 1:ncomp2, drop=FALSE])^2
    } else {
      if(ncomp1 == 1){
        smi <- t(
          apply(
            crossprod(Scores1[, 1:ncomp1, drop=FALSE], Scores2[, 1:ncomp2, drop=FALSE])^2,
            1,cumsum)) / matrix(pmin(rep(1:ncomp1,ncomp2), rep(1:ncomp2,each=ncomp1)), ncomp1,ncomp2)
      } else {
        if(ncomp2 == 1){
          smi <- apply(
            crossprod(Scores1[, 1:ncomp1, drop=FALSE], Scores2[, 1:ncomp2, drop=FALSE])^2,
          1,cumsum) / matrix(pmin(rep(1:ncomp1,ncomp2), rep(1:ncomp2,each=ncomp1)), ncomp1,ncomp2)
        }
      }
    }
    attr(smi, "orthogonal") <- TRUE
    smi[smi > 1] <- 1; smi[smi < 0] <- 0
  } else { # Procrustes Rotation
    TU <- crossprod(Scores1,Scores2)
    for(p in 1:ncomp1){
      for(q in 1:ncomp2){
        s <- svd(TU[1:p,1:q, drop=FALSE])$d
        smi[p,q] <- mean(s)^2
      }
    }
    attr(smi, "orthogonal") <- FALSE
  }
  if(!has.scores){
    attr(smi,'PCA') <- TRUE
  } else {
    attr(smi,'PCA') <- FALSE
  }
  attr(smi, "n")         <- nobj
  attr(smi, "mat.names") <- mat.names
  attr(smi, "scores")    <- list(Scores1=Scores1, Scores2=Scores2)
  dn            <- list(1:ncomp1, 1:ncomp2)
  names(dn)     <- mat.names
  dimnames(smi) <- dn
  class(smi)    <- c("SMI", "matrix")
  smi
}
