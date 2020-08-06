#' Procrustes Similarity Index
#' 
#' An index based on the RV coefficient with Procrustes rotation.
#'
#' @param X1 first \code{matrix} to be compared (\code{data.frames} are also accepted).
#' @param X2 second \code{matrix} to be compared (\code{data.frames} are also accepted).
#' @param center \code{logical} indicating if input matrices should be centered (default = TRUE).
#'
#' @return The Procrustes Similarity Index
#' 
#' @references Sibson, R; 1978. "Studies in the Robustness of Multidimensional Scaling: Procrustes Statistics". Journal of the Royal Statistical Society. Series B (Methodological), Vol. 40, No. 2, pp. 234-238.
#'
#' @examples
#' X1  <- scale( matrix( rnorm(100*300), 100,300), scale = FALSE)
#' usv <- svd(X1)
#' X2  <- usv$u[,-3] %*% diag(usv$d[-3]) %*% t(usv$v[,-3])
#' PSI(X1,X2)
#' 
#' @importFrom pracma procrustes
#' @export
PSI <- function(X1, X2, center = TRUE){
  X1 <- as.matrix(X1)
  X2 <- as.matrix(X2)
  if(center){
    X1 <- X1 - rep(colMeans(X1), each = nrow(X1))
    X2 <- X2 - rep(colMeans(X2), each = nrow(X1))
  }
  tr <- procrustes(as.matrix(X1),as.matrix(X2))$Q
  Trace(crossprod(X1,X2)%*%tr)/sqrt(Trace(crossprod(X1))*Trace(crossprod(X2)))
}
