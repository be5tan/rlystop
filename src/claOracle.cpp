#include <Rcpp.h>
using namespace Rcpp;

//' Classical oracle
//'
//' Computes the classical oracle index.
//'
//' @param lambda Vector of decreasing, strictly positive entries of the
//'   diagonal design matrix.
//' @param mu Vector valued input signal.
//' @param delta Numeric noise level.
//' @param alpha Numeric smoothing parameter. \code{alpha} = -1 gives
//'   the strong classical oracle. \code{alpha} = 0 gives the weak classical
//'   oracle. 
//' @param filt Character string designating the filter to be used. filt should
//'   be one of "cutoff" or "landw".
//'
//' @return Returns the classical oracle index as an integer.
//'
//' @export
// [[Rcpp::export]]
int claOracle(NumericVector lambda, NumericVector mu, double delta, double
        alpha = -1, std::string filt = "cutoff") {
    if (filt != "cutoff" && filt != "landw") {
      Rcout << "Error: filt should be one of \"cutoff\", \"landw\"" << std::endl;
    }

    int m_amin = 0;
    int D = mu.length();
    double delta2 = pow(delta, 2);
    NumericVector B2 (D);
    NumericVector V (D);
    NumericVector MSE (D);

    if (filt == "cutoff") {
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
    }

    if (filt == "landw") {
        NumericVector auxFilterTerm = 1 - pow(lambda, 2);
        NumericVector biasFilterTerm = pow(auxFilterTerm, 2);
        NumericVector varFilterTerm  = pow(1 - auxFilterTerm, 2);
        NumericVector biasSmoothingTerm = pow(lambda, 2 + 2 * alpha);
        NumericVector varSmoothingTerm = pow(lambda, 2 * alpha);
        B2[0] = sum(biasFilterTerm * biasSmoothingTerm * pow(mu, 2));
        V[0]  = sum(varFilterTerm  * varSmoothingTerm) * pow(delta, 2);
        MSE[0] = B2[0] + V[0];
        for (int m = 1; m < D; m++) {
            auxFilterTerm  = auxFilterTerm * (1 - pow(lambda, 2));
            biasFilterTerm = pow(auxFilterTerm, 2);
            varFilterTerm  = pow(1 - auxFilterTerm, 2);
            B2[m]  = sum(biasFilterTerm * biasSmoothingTerm * pow(mu, 2));
            V[m]   = sum(varFilterTerm * varSmoothingTerm) * pow(delta, 2);
            MSE[m] = B2[m] + V[m];
            if (MSE[m] < MSE[m_amin]) {
                m_amin = m;
            }
        }
    }

    return(m_amin + 1);
}
