// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>  
using namespace Rcpp;
using arma::mat;
using arma::vec;

// [[Rcpp::export]]
arma::mat cmr(arma::mat psi, arma::mat X, int n, arma::mat K, int reps) {
  RNGScope scope; // ensure RNG gets set/reset
  int nr = X.n_rows, nc = X.n_cols;
  
  // SVD of covariance matrix
  mat U;
  vec s;
  mat V;
  svd(U, s, V, X);
  mat tU = U.t();
  // covariance structure
  mat R = trans(V * (tU.each_col() % sqrt(s)));
  
  int nk = K.n_rows;
  mat stat(reps, nk);
  for (int r = 0; r < reps; r++) { // iterate number of random rotations
    // rotation matrix
    mat rnm(n, nc);
    for (int i = 0; i < n; i++) {
      for (int j = 0; j < nc; j++) {
        rnm(i,j) = rnorm(1,0,1)[0];
      }
    }
    // random rotation
    mat G = rnm * R;
    // calculate test statistics
    mat Gpsi = trans(G % psi.t());
    vec est = K * Gpsi * arma::ones(n);
    mat covi = K * ((Gpsi * Gpsi.t()) / n) * K.t();
    vec stati = (est / sqrt(n)) / sqrt(diagvec(covi));
    // store test statistics in output matrix
    for (int k = 0; k < nk; k++){
      stat(r,k) = stati(k);
    }
  }
  // return matrix with rotated test statistics
  return stat;  
}