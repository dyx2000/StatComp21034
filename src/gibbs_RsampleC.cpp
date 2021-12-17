#include <Rcpp.h>
using namespace Rcpp;
//' @title A Gibbs sampler using Rcpp
//' @description A Gibbs sampler using Rcpp
//' @param N the number of samples
//' @param thin the number of between-sample random numbers
//' @param a parameter in the distribution
//' @param b parameter in the distribution
//' @return a random sample of size \code{n}
//' @export


// [[Rcpp::export]]
NumericMatrix gibbs_RsampleC(int N, double a, double b, int n, int thin) {
  NumericMatrix mat(N, 2);
  double x = floor(n/2), y = 0.5;
  for(int i = 0; i < N; i++) {
    for(int j = 0; j < thin; j++) {
      x = rbinom(1, n, y )[0];
      y = rbeta(1, x+a, n-x+b)[0];
    }
    mat(i, 0) = x;
    mat(i, 1) = y;
  }
  return(mat);
}