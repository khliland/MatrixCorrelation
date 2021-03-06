#' @title Similiarity of Matrices Coefficients
#' @description Computation and visualization of matrix correlation coefficients.
#' The main method is the Similarity of Matrices Index, while various related
#' measures like r1, r2, r3, r4, Yanai's GCD, RV, RV2, adjusted RV, Rozeboom's
#' linear correlation and Coxhead's coefficient are included
#' for comparison and flexibility.
#' @references
#' \itemize{
#'  \item{SMI:}{ Indahl, U.G.; Næs, T.; Liland, K.H.; 2018. A similarity index for comparing coupled matrices. Journal of Chemometrics; e3049.}
#'  \item{RV:}{ Robert, P.; Escoufier, Y. (1976). "A Unifying Tool for Linear Multivariate
#'   Statistical Methods: The RV-Coefficient". Applied Statistics 25 (3): 257-265.}
#'  \item{RV2:}{ Smilde, AK; Kiers, HA; Bijlsma, S; Rubingh, CM; van Erk, MJ (2009). "Matrix correlations
#'  for high-dimensional data: the modified RV-coefficient". Bioinformatics 25(3): 401-5.}
#'  \item{Adjusted RV:}{ Mayer, CD; Lorent, J; Horgan, GW. (2011). "Exploratory analysis of multiple omics
#'  datasets using the adjusted RV coefficient". Stat Appl Genet Mol Biol. 10(14).}
#'  \item{PSI:}{ Sibson, R; 1978. "Studies in the Robustness of Multidimensional Scaling: 
#'  Procrustes Statistics". Journal of the Royal Statistical Society. Series B (Methodological), Vol. 40, No. 2, pp. 234-238.}
#'  \item{Rozeboom:}{ Rozeboom, WW; 1965. "Linear correlations between sets of variables". 
#'  Psychometrika 30(1): 57-71.}
#'  \item{Coxhead:}{ Coxhead, P; 1974. "Measuring the releationship between two sets of variables". British Journal of 
#'  Mathematical and Statistical Psychology 27: 205-212.}
#' }
#' @seealso \code{\link{SMI}}, \code{\link{plot.SMI}} (print.SMI/summary.SMI), \code{\link{RV}} (RV2/RVadj), \code{\link{r1}} (r2/r3/r4/GCD), \code{\link{Rozeboom}}, \code{\link{Coxhead}}, \code{\link{allCorrelations}} (matrix correlation comparison).
#' @name MatrixCorrelation
"_PACKAGE"
