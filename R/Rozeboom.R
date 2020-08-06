#' Rozeboom's squared vector correlation
#'
#' @param X1 first \code{matrix} to be compared (\code{data.frames} are also accepted).
#' @param X2 second \code{matrix} to be compared (\code{data.frames} are also accepted).
#'
#' @return A single value measuring the similarity of two matrices. For diagnostic purposes it is accompanied by an attribute "canonical.correlation".
#'
#' @author Korbinian Strimmer and Kristian Hovde Liland
#'
#' @references Rozeboom, WW; 1965. "Linear correlations between sets of variables". Psychometrika 30(1): 57-71.
#'
#' @seealso \code{\link{SMI}}, \code{\link{RV}} (RV2/RVadj), \code{\link{Coxhead}}, \code{\link{r1}} (r2/r3/r4/GCD).
#'
#' @examples
#' X <- matrix(rnorm(100*13),nrow=100)
#' X1 <- X[, 1:5]  # Random normal
#' X2 <- X[, 6:12] # Random normal
#' X2[,1] <- X2[,1] + X[,5] # Overlap in one variable
#' Rozeboom(X1, X2)
#'
#' @export
Rozeboom <- function(X1, X2){
  X1 <- as.matrix(X1)
  X2 <- as.matrix(X2)
  lambda <- cancor(X1, X2)$cor # Canonical correlations
  alien  <- prod(1-lambda^2) # Hotelling (1936) vector alienation
  R <- 1-alien
  attr(R, 'canonical.correlations') <- lambda

  return( R )
}

#' @rdname Rozeboom
#' @export
sqveccor <- function(X1, X2){ # Synonym for Rozeboom
  lambda <- cancor(X1, X2)$cor # Canonical correlations
  alien  <- prod(1-lambda^2) # Hotelling (1936) vector alienation

  return( 1-alien )
}
