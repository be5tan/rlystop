#include <Rcpp.h>
using namespace Rcpp;

//' Squared bias
//' 
//' Computes the squared bias of the filter estimator a given diagonal design
//' matrix \code{lambda} and signal ' \code{mu}.
//'
//' @param lambda Vector of decreasing, strictly positive entries of the
//'   diagonal design matrix.
//' @param mu Vector valued input signal. 
//' @param alpha Numeric smoothing parameter.
//' @param filt Character string giving the filter to be used.
//'
//' @return Returns all values of the squared bias in a double vector.
//'
//' @export
// [[Rcpp::export]]
NumericVector bias2(NumericVector lambda, NumericVector mu, double alpha
        = -1.0, std::string filt = "cutoff")
{
    int D = mu.length();
    NumericVector B2(D);
    B2[0] = sum(pow(lambda, 2 + 2 * alpha) * pow(mu, 2)) -
        pow(lambda[0], 2 + 2 * alpha) * pow(mu[0], 2);
    for (int m = 1; m < D; ++m) {
        B2[m] = B2[m - 1] -
            pow(lambda[m], 2 + 2 * alpha) * pow(mu[m], 2);
    }
    return(B2);
}
