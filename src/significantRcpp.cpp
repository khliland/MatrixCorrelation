// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;

// inline int randWrapper(const int n) { return R::ftrunc(R::runif(0,n)); }

// [[Rcpp::export]]
NumericVector significantRcpp(SEXP smi1, SEXP T1, SEXP U1, SEXP B1) {
  arma::mat T   = Rcpp::as<arma::mat>(T1);
  arma::mat U   = Rcpp::as<arma::mat>(U1);
  arma::mat smi = Rcpp::as<arma::mat>(smi1);
  arma::uword B = as<int>(B1);
  arma::uword ncomp1 = T.n_cols;
  arma::uword ncomp2 = U.n_cols;
  double csum = 0;
  arma::mat P(ncomp1,ncomp2, arma::fill::zeros);
  arma::mat m(ncomp1,ncomp2, arma::fill::zeros);
  arma::vec TU(ncomp1);

  for(arma::uword i=0; i<ncomp1; ++i){
    for(arma::uword j=0; j<ncomp2; ++j){
      m(i,j) = (double)i+1;
      if(j<i){
        m(i,j) = (double)j+1;
      }
    }
  }

  for(arma::uword b=0; b<B; ++b){
    auto rng = std::default_random_engine {};
    std::shuffle ( T.begin(), T.end(), rng );
    TU.zeros();
    for(arma::uword j=0; j<ncomp2; ++j){
      csum = 0;
      for(arma::uword i=0; i<ncomp1; ++i){
        csum += pow(arma::as_scalar(T.col(i).t()*U.col(j)),2);
        TU(i) += csum;
        if(std::max(TU(i)/m(i,j),1-TU(i)/m(i,j)) < smi(i,j)){
          P(i,j) += 1;
        }
      }
    }
  }
  return(wrap(P));
}

