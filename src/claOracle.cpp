#include <Rcpp.h>
using namespace Rcpp;

//' Classical oracle
//'
//' Computes the classical oracle time.
//'
//' @param lambda Vector of decreasing, strictly positive entries of the
//'   diagonal design matrix.
//' @param mu Vector valued input signal.
//' @param delta Strictly positive level of random noise.
//' @param filt Character string giving the filter to be used. This must match
//'   one of "cutoff".
//' @param alpha Numeric smoothing parameter.
//'
//' @return Returns the classical oracle time as an integer.
//'
//' @export
// [[Rcpp::export]]
int claOracle(NumericVector lambda, NumericVector mu, double delta, double
        alpha = -1, std::string filt = "cutoff") {
    int m_amin = 0;
    int D = mu.length();
    double delta2 = pow(delta, 2);
    NumericVector B2 (D);
    NumericVector V (D);
    NumericVector MSE (D);
    B2[0] = sum(pow(lambda[Range(1, D - 1)], 2 + 2 * alpha) * pow(mu[Range(1, D - 1)], 2));
    V[0] = delta2 * pow(lambda[0], 2 * alpha);
    MSE[0] = B2[0] + V[0];
    for (int m = 1; m < D; m++) {
        V[m] = V[m - 1] + pow(lambda[m], 2 * alpha) * delta2;
        B2[m] = B2[m - 1] - pow(lambda[m], 2 + 2 * alpha) * pow(mu[m], 2);
        MSE[m] = B2[m] + V[m];
        if (MSE[m] < MSE[m_amin]) {
            m_amin = m;
        }
    }
    return(m_amin + 1);
}
