// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;

// [[Rcpp::export]]
double RV2cpp(arma::mat X, arma::mat Y){
  arma::mat AA = X * X.t();
  arma::mat BB = Y * Y.t();
  AA.diag() *= 0.0;
  BB.diag() *= 0.0;
  return trace(AA*BB) / (sqrt(accu(pow(AA,2))) * sqrt(accu(pow(BB,2))));
}
