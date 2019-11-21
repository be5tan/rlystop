#include <Rcpp.h>
using namespace Rcpp;

//' Landweber iteration
//'
//' Computes the Landweber iteration for a diagon
//' 
//' @param m Integer stopping index. 
//' @param Y Numeric vector of observed data.
//' @param lambda Numeric vector of decreasing, strictly positive entries of
//'   the diagonal design matrix.
//'
//' @return Returns the m-th iterate of the Landweber iteration.
//'
//' @export
// [[Rcpp::export]]
NumericVector landw(int m, NumericVector Y, NumericVector lambda)
{
    int D = Y.length();
    int iter = 0;
    NumericVector muHat(D);
    while (iter < m) {
        muHat += lambda * (Y - lambda * muHat); 
        iter += 1;
    }
    return(muHat);
}
