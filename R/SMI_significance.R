#' @title Significance estimation for Similarity of Matrices Index (SMI)
#'
#' @description Permutation based hypothesis testing for SMI. The nullhypothesis is that a linear function
#' of one matrix subspace is included in the subspace of another matrix.
#'
#' @param smi \code{smi} object returned by call to \code{SMI}.
#' @param B integer number of permutations, default = 10000.
#'
#' @details For each combination of components significance is estimated by sampling from a null distribution
#' of no similarity, i.e. when the rows of one matrix is permuted B times and corresponding SMI values are
#' computed.
#'
#' @return A matrix containing P-values for all combinations of components.
#'
#' @author Kristian Hovde Liland
#'
#' @references Similarity of Matrices Index - Ulf G. Indahl, Tormod Næs, Kristian Hovde Liland
#'
#' @seealso \code{\link{plot.SMI}} (print.SMI/summary.SMI), \code{\link{RV}} (RV2/RVadj), \code{\link{r1}} (r2/r3/r4/GCD), \code{\link{allCorrelations}} (matrix correlation comparison).
#'
#' @examples
#' X1  <- matrix(rnorm(100*300),100,300)
#' usv <- svd(X1)
#' X2  <- usv$u[,-3] %*% diag(usv$d[-3]) %*% t(usv$v[,-3])
#'
#' (smi <- SMI(X1,X2,5,5))
#' significant(smi)
#'
#' @importFrom progress progress_bar
#' @export
significant <- function(smi, B = 10000){
  T <- attr(smi,'scores')$Scores1
  U <- attr(smi,'scores')$Scores2
  
  # Orthogonal Projection
  if(attr(smi, "orthogonal")){ 
    return(1-significantRcpp(smi, T, U, B)/B)
  
  # Procrustes Rotation
  } else { 
    N     <- dim(T)[1]
    ncomp <- dim(smi)
    smi_B <- array(0, dim = c(ncomp,B))
    pb <- progress_bar$new(total = B, format = "  [:bar] :percent (:eta)")
    for(b in 1:B){
      ind = sample(N)
      TU <- crossprod(T[ind,],U)
      for(p in 1:ncomp[1]){
        for(q in 1:ncomp[2]){
          s <- svd(TU[1:p,1:q, drop=FALSE], nu = 0, nv = 0)$d
          smi_B[p,q,b] <- mean(s)^2
        }
      }
      pb$tick()
    }
    Pval <- matrix(0, ncomp[1], ncomp[2])
    for(p in 1:ncomp[1]){
      for(q in 1:ncomp[2]){
        Pval[p,q] <- 1-sum(smi_B[p,q,]<(1-smi[p,q]))/B
      }
    }
    return(Pval)
  }
}

