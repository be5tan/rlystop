#include <Rcpp.h>
using namespace Rcpp;

//' Balanced oracle 
//' 
//' Computes the balanced oracle time. 
//'
//' @param lambda Vector of decreasing, strictly positive entries of the
//'   diagonal design matrix.
//' @param mu Vector valued input signal. 
//' @param delta Strictly positive level of random noise.
//' @param alpha Numeric smoothing parameter.
//' @param filt Character string giving the filter to be used.
//'
//' @return Returns the strong balanced oracle time as an integer.
//'
//' @export
// [[Rcpp::export]]
int balOracle(NumericVector lambda, NumericVector mu, double delta, double alpha
        = -1.0, std::string filt = "cutoff")
{
    int D = lambda.length();
    int m = 0;
    double B2_m_alpha = sum(pow(lambda, 2 + 2 * alpha) * pow(mu, 2));
    double V_m_alpha = 0;
    while (V_m_alpha < B2_m_alpha) {
       V_m_alpha += pow(lambda[m], 2 * alpha) * pow(delta, 2);
       B2_m_alpha = sum(pow(lambda[Range(m + 1, D - 1)], 2 + 2 * alpha) *
               pow(mu[Range(m + 1, D - 1)], 2));
       m += 1;
    }
    return(m);
}
