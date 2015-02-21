#include <Rcpp.h>

using namespace Rcpp;

// [[Rcpp::depends(Rcpp)]]

/*
*   dinvgamma1_Rcpp: log density of the inverse gamma distribution
*   ddinvgamma1_Rcpp: first derivative of the log density of the inverse gamma distribution
*   d2dinvgamma1_Rcpp: second derivative of the log density of the inverse gamma distribution
*/

// [[Rcpp::export]]
double dinvgamma1_Rcpp(const double x, const double a, const double b){
  return a * log( b ) - lgamma( a ) - ( a + 1 ) * log( x ) - b / x ;
}

// [[Rcpp::export]]
double ddinvgamma1_Rcpp(const double x, const double a, const double b){
   return -( a + 1 ) / x + b / ( x * x ) ;
}

// [[Rcpp::export]]
double d2dinvgamma1_Rcpp(const double x, const double a, const double b){
    return ( a + 1 ) / ( x * x ) - 2 * b / (x * x * x) ;
}


