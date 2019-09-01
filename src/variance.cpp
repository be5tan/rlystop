#include <Rcpp.h>
using namespace Rcpp;

//' Variance
//' 
//' Computes the variance of the filter estimator a given diagonal design
//' matrix \code{lambda} and noise level delta.
//'
//' @param lambda Vector of decreasing, strictly positive entries of the
//'   diagonal design matrix.
//' @param delta Numeric noise level.
//' @param alpha Numeric smoothing parameter.
//' @param filt Character string giving the filter to be used.
//'
//' @return Returns all values of the variancee in a double vector.
//'
//' @export
// [[Rcpp::export]]
NumericVector variance(NumericVector lambda, double delta, double alpha = -1.0,
        std::string filt = "cutoff")
{
    int D = lambda.length();
    NumericVector V(D);
    V[0] = pow(lambda[0], 2 * alpha) * pow(delta, 2);
    for (int m = 1; m < D; ++m) {
        V[m] = V[m - 1] + pow(lambda[m], 2 * alpha) * pow(delta, 2);
    }
    return(V);
}
