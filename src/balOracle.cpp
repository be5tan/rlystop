#include <Rcpp.h>
using namespace Rcpp;

//' Balanced oracle
//'
//' Computes the balanced oracle index.
//'
//' @param lambda Vector of decreasing, strictly positive entries of the
//'   diagonal design matrix.
//' @param mu Vector valued input signal.
//' @param delta Numeric noise level.
//' @param alpha Numeric smoothing parameter. \code{alpha} = -1 gives
//'   the strong balanced oracle. \code{alpha} = 0 gives the weak balanced
//'   oracle. 
//' @param filt Character string designating the filter to be used. filt should
//'   be one of "cutoff" or "landw".
//'
//' @return Returns the strong balanced oracle index as an integer.
//'
//' @export
// [[Rcpp::export]]
int balOracle(NumericVector lambda, NumericVector mu, double delta, double alpha
        = -1.0, std::string filt = "cutoff")
{
    if (filt != "cutoff" && filt != "landw") {
      Rcout << "Error: filt should be one of \"cutoff\", \"landw\"" << std::endl;
    }

    int D = lambda.length();
    int m = 0;

    if (filt == "cutoff") {
        double B2_m_alpha = sum(pow(lambda, 2 + 2 * alpha) * pow(mu, 2));
        double V_m_alpha = 0;
        while (V_m_alpha < B2_m_alpha) {
           V_m_alpha += pow(lambda[m], 2 * alpha) * pow(delta, 2);
           B2_m_alpha = sum(pow(lambda[Range(m + 1, D - 1)], 2 + 2 * alpha) *
                   pow(mu[Range(m + 1, D - 1)], 2));
           m += 1;
        }
    }

    if (filt == "landw") {
        NumericVector auxFilterTerm(D, 1.0);
        NumericVector biasFilterTerm(D);
        NumericVector varFilterTerm(D);
        NumericVector biasSmoothingTerm = pow(lambda, 2 + 2 * alpha);
        NumericVector varSmoothingTerm = pow(lambda, 2 * alpha);
        double B2_m_alpha = sum(biasSmoothingTerm * pow(mu, 2));
        double V_m_alpha = 0;
        while (V_m_alpha < B2_m_alpha) {
            auxFilterTerm  = auxFilterTerm * (1 - pow(lambda, 2));
            // auxFilterTerm  = pow(1 - pow(lambda, 2), m + 1);
            biasFilterTerm = pow(auxFilterTerm, 2);
            varFilterTerm  = pow(1 - auxFilterTerm, 2);
            V_m_alpha  = sum(varFilterTerm * varSmoothingTerm) * pow(delta, 2);
            B2_m_alpha = sum(biasFilterTerm * biasSmoothingTerm * pow(mu, 2));
            m += 1;
        }
    }

    return(m);
}
