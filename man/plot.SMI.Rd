% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/SMI_results.R
\name{plot.SMI}
\alias{plot.SMI}
\alias{print.SMI}
\alias{summary.SMI}
\alias{is.signif}
\title{Result functions for the Similarity of Matrices Index (SMI)}
\usage{
\method{plot}{SMI}(
  x,
  y = NULL,
  x1lab = attr(x, "mat.names")[[1]],
  x2lab = attr(x, "mat.names")[[2]],
  main = "SMI",
  signif = 0.05,
  xlim = c(-(pq[1] + 1)/2, (pq[2] + 1)/2),
  ylim = c(0.5, (sum(pq) + 3)/2),
  B = 10000,
  cex = 1,
  cex.sym = 1,
  frame = NULL,
  frame.col = "red",
  frame.lwd = 2,
  replicates = NULL,
  ...
)

\method{print}{SMI}(x, ...)

\method{summary}{SMI}(object, ...)

is.signif(x, signif = 0.05, B = 10000, ...)
}
\arguments{
\item{x}{object of class \code{SMI}.}

\item{y}{not used.}

\item{x1lab}{optional label for first matrix.}

\item{x2lab}{optional label for second matrix.}

\item{main}{optional heading (default = SMI).}

\item{signif}{significance level for testing (default=0.05).}

\item{xlim}{optional plotting limits.}

\item{ylim}{optional plotting limits.}

\item{B}{number of permutations (for significant, default=10000).}

\item{cex}{optional text scaling (default = 1)}

\item{cex.sym}{optional scaling for significance symbols (default = 1)}

\item{frame}{two element integer vector indicating framed components.}

\item{frame.col}{color for framed components.}

\item{frame.lwd}{line width for framed components.}

\item{replicates}{vector of replicates for significance testing.}

\item{...}{additional arguments for \code{plot}.}

\item{object}{object of class \code{SMI}.}
}
\value{
\code{plot} silently returns NULL. \code{print} and \code{summary} return the printed matrix.
}
\description{
Plotting, printing and summary functions for SMI, plus significance testing.
}
\details{
For plotting a diamonad plot is used. High SMI values are light and low SMI values
are dark. If orthogonal projections have been used for calculating SMIs, significance symbols
are included in the plot unless signif=NULL.
}
\examples{
X1  <- scale( matrix( rnorm(100*300), 100,300), scale = FALSE)
usv <- svd(X1)
X2  <- usv$u[,-3] \%*\% diag(usv$d[-3]) \%*\% t(usv$v[,-3])

smi <- SMI(X1,X2,5,5)
plot(smi, B = 1000) # default B = 10000
print(smi)
summary(smi)
is.signif(smi, B = 1000) # default B = 10000

}
\references{
Similarity of Matrices Index - Ulf G. Indahl, Tormod Næs, Kristian Hovde Liland
}
\seealso{
\code{\link{SMI}}, \code{\link{PCAcv} (cross-validated PCA)}.
}
\author{
Kristian Hovde Liland
}
