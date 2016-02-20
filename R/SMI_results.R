#' @aliases print.SMI summary.SMI is.signif
#' @title Result functions for the Similarity of Matrices Index (SMI)
#'
#' @description Plotting, printing and summary functions for SMI, plus significance testing.
#'
#' @param x object of class \code{SMI}.
#' @param object object of class \code{SMI}.
#' @param y not used.
#' @param x1lab optional label for first matrix.
#' @param x2lab optional label for second matrix.
#' @param main optional heading (default = SMI).
#' @param signif significance level for testing (default=0.05).
#' @param xlim optional plotting limits.
#' @param ylim optional plotting limits.
#' @param B number of permutations (for significant, default=10000).
#' @param cex optional text scaling (default = 1)
#' @param cex.sym optional scaling for significance symbols (default = 1)
#' @param ... additional arguments for \code{plot}.
#'
#' @details For plotting a diamonad plot is used. High SMI values are light and low SMI values
#' are dark. If orthogonal projections have been used for calculating SMIs, significance symbols
#' are included in the plot unless signif=NULL.
#'
#' @return \code{plot} silently returns NULL. \code{print} and \code{summary} return the printed matrix.
#'
#' @author Kristian Hovde Liland
#'
#' @references Similarity of Matrices Index - Ulf G. Indahl, Tormod Næs, Kristian Hovde Liland
#'
#' @seealso \code{\link{SMI}}.
#'
#' @examples
#' X1  <- matrix(rnorm(100*300),100,300)
#' usv <- svd(X1)
#' X2  <- usv$u[,-3] %*% diag(usv$d[-3]) %*% t(usv$v[,-3])
#'
#' smi <- SMI(X1,X2,5,5)
#' plot(smi)
#' print(smi)
#' summary(smi)
#' is.signif(smi)
#'
#' @export
plot.SMI <- function(x, y = NULL, x1lab = attr(x, "mat.names")[[1]], x2lab = attr(x, "mat.names")[[2]],
                     main = "SMI", signif = 0.05,
                     xlim = c(-(pq[1])/2,(pq[2]+1)/2), ylim = c(0.75,(sum(pq)+2.5)/2), B = 10000, cex = 1, cex.sym = 1, ...){
  pq <- dim(x)

  # Check inputs/defaults
  Procrustes <- ifelse(attr(x,"orthogonal"), FALSE, TRUE)

  dev.hold()
  on.exit(dev.flush())

  mar.old <- par("mar") # Outer plot margins
  par(mar = c(2,3,2,1))
  pty.old <- par("pty") # Squared plotting region
  par(pty="s")

  # Prepare plotting
  plot(0,1.5,col='white', xlim = xlim, ylim = ylim, main = main, axes = FALSE, xlab='', ylab='', cex=cex, ...)

  # Polygon loops
  for(p in 1:pq[1]){
    for(q in 1:pq[2]){
      x1 <- -(p-q)/2; x2 <- (p+q)/2;
      polygon(x1+c(0, 0.5, 0, -0.5), x2+c(0, 0.5, 1, 0.5), col = rgb(x[p,q],x[p,q],x[p,q]), border=NA, )
    }
  }

  # Plot significance symbols
  if(!Procrustes && !is.null(signif)){
    Pval <- significant(smi, B)
    n <- attr(x, "n")
    for(p in 1:pq[1]){
      for(q in 1:pq[2]){
        x1 <- -(p-q)/2; x2 <- (p+q)/2;
        if(Pval[p,q] > signif){
          if(p < q){
            text(x1, x2+0.5, "U", adj=0.5, srt=-90, cex=cex*cex.sym, ...)
          } else {
            if(p > q){
              text(x1, x2+0.5, "U", adj=0.5, srt=90, cex=cex*cex.sym, ...)
            } else { # p == q
              text(x1, x2+0.5, "=", adj=0.5, cex=cex*cex.sym, ...)
            }
          }
        } else {
          if(Pval[p,q] < 0.001){
            text(x1, x2+0.5, "***", adj=0.5, cex=cex*cex.sym, ...)
          } else {
            if(Pval[p,q] < 0.01){
              text(x1, x2+0.5, "**", adj=0.5, cex=cex*cex.sym, ...)
            } else {
              if(Pval[p,q] < 0.05){
                text(x1, x2+0.5, "*", adj=0.5, cex=cex*cex.sym, ...)
              } else {
                if(Pval[p,q] < 0.1){
                  text(x1, x2+0.5, ".", adj=0.5, cex=cex*cex.sym, ...)
                }
              }
            }
          }
        }
      }
    }
  }

  # Lines and numbers
  for(p in 0:pq[2]){
    lines( 0.5*p+c(0,-0.5*pq[1]), 1+0.5*p+c(0,0.5*pq[1]))
    if(p > 0){
      text(-0.15+p*0.5, 0.5+0.5*p, p, adj=0, cex=cex, ...)
    }
  }
  for(q in 0:pq[1]){
    lines(-0.5*q+c(0, 0.5*pq[2]), 1+0.5*q+c(0,0.5*pq[2]))
    if(q > 0){
      text(0.15-q*0.5,  0.5+0.5*q, q, adj=1, cex=cex, ...)
    }
  }

  # Legend
  cols <- seq(0,1,length.out=200)
  color.legend(xlim[2]+diff(xlim)*0.06, ylim[1], xlim[2]+diff(xlim)*0.11, ylim[2],
               seq(0,1,by=0.1), rgb(cols,cols,cols), align="lt", gradient="y", cex=cex*0.8)

  # Axis labels
  text(3*pq[2]/8, 1+pq[2]/8, x2lab, adj=0.5, cex=cex, ...)
  text(-3*pq[1]/8,1+pq[1]/8, x1lab, adj=0.5, cex=cex, ...)
  par(pty=pty.old)
  par(mar=mar.old)
  invisible(NULL)
}

#' @rdname plot.SMI
#' @export
print.SMI <- function(x, ...){
  class(x) <- "matrix"
  attr(x, "mat.names")  <- NULL
  attr(x, "orthogonal") <- NULL
  attr(x, "n")          <- NULL
  attr(x, "PCA")        <- NULL
  attr(x, "scores")     <- NULL
  print(x)
}

#' @rdname plot.SMI
#' @export
summary.SMI <- function(object, ...){
  class(object) <- "matrix"
  attr(object, "mat.names")  <- NULL
  PCA <- ifelse(attr(object, "PCA"), 'PCA', 'custom')
  if(!attr(object, "orthogonal")){
    cat(paste("SMI with ", PCA, " subspace and Procrustes Rotation\n", sep=""))
  } else {
    cat(paste("SMI with ", PCA, " subspace and Orthogonal Projection\n", sep=""))
  }
  attr(object, "orthogonal") <- NULL
  attr(object, "n")          <- NULL
  attr(object, "PCA")        <- NULL
  attr(object, "scores")     <- NULL
  print(object)
}

#' @rdname plot.SMI
#' @export
is.signif <- function(x, signif = 0.05, B = 10000){
  # Check inputs/defaults
  if(!attr(x,"orthogonal"))
    stop('Significance testing is only supported for Orthogonal Projections')
  pq <- dim(x)
  sig <- matrix(FALSE, pq[1], pq[2])
  n <- attr(x, "n")

  Pval <- significant(x, B)
  for(p in 1:pq[1]){
    for(q in 1:pq[2]){
      if(Pval[p,q] > signif){
        sig[p,q] <- TRUE
      }
    }
  }
  dimnames(sig) <- dimnames(x)
  sig
}

