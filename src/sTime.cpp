#include <Rcpp.h>
using namespace Rcpp;

//' Residual stopping time
//'
//' Computes the residual based stopping time for a given data
//' vector \code{Y} and a given diagonal design matrix \code{lambda}.
//'
//' @param Y Numeric vector of observed data.
//' @param lambda Numeric vector of decreasing, strictly positive entries of the
//'   diagonal design matrix.
//' @param alpha Numeric smoothing index for the residuals.
//' @param kappa Strictly positive stopping value.
//' @param filt Character string giving the filter to be used.
//'
//' @return Returns the integer value of the smoothed residual based stopping
//'   time.
//'
//' @export
// [[Rcpp::export]]
int sTime(NumericVector Y, NumericVector lambda, double alpha, double kappa,
        std::string filt = "cutoff")
{
    if (filt != "cutoff" && filt != "landw") {
      Rcout << "Error: filt should be one of \"cutoff\", \"landw\"" << std::endl;
    }

    int tau = 0;
    int D = Y.length();
    NumericVector muHat(D, 0.0);
    NumericVector smoothingTerm = pow(lambda, 2 * alpha);
    double residuals2 = sum(smoothingTerm * pow(Y, 2));
    // double residuals2 = sum(pow(lambda[Range(tau, D - 1)], 2 * alpha) *
    //         pow(Y[Range(tau, D - 1)], 2));

    if (filt == "cutoff") {
        while (residuals2 > kappa) {
            residuals2 -= pow(lambda[tau], 2 * alpha) * pow(Y[tau], 2);
            tau += 1;
        }
    }

    if (filt == "landw") {
        while (residuals2 > kappa) {
               muHat += lambda * (Y - lambda * muHat);
               residuals2 = sum(smoothingTerm * pow(Y - lambda * muHat, 2));
               tau += 1;
        }
    }

    return(tau);
}
