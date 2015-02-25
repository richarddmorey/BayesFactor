#ifndef BFCOMMON_HPP_
#define BFCOMMON_HPP_

#include <RcppEigen.h>

using namespace Rcpp;

int RcppCallback(time_t *last, Rcpp::Function cb, double progress, double callbackInterval);

double dinvgamma1_Rcpp(const double x, const double a, const double b);
double ddinvgamma1_Rcpp(const double x, const double a, const double b);
double d2dinvgamma1_Rcpp(const double x, const double a, const double b);

double log_determinant_pos_def(Eigen::MatrixXd A);

List jzs_log_marginal_posterior_logg(const NumericVector q, const double sumSq, const int N, const NumericMatrix XtCnX, const NumericMatrix CnytCnX0, const NumericVector rscale, const IntegerVector gMap, const NumericVector gMapCounts, const NumericMatrix priorX, const int incCont, const bool limit, const NumericVector limits, const int which);

double jzs_mc_marg_like(const NumericVector g, const double sumSq, const int N, const NumericMatrix XtCnX, const NumericMatrix CnytCnX, const NumericVector rscale, const IntegerVector gMap, const NumericVector gMapCounts, const NumericMatrix priorX, const double logDetPriorX, const int incCont);
double jzs_mc_marg_like2(const NumericVector g, const double sumSq, const int N, const NumericMatrix XtCnX, const NumericMatrix XtCny, const IntegerVector gMap, const NumericMatrix priorX, const double logDetPriorX, const int incCont);

double jzs_importance_marg_like(const NumericVector q, const NumericVector mu, const NumericVector sig, const double sumSq, const int N, const NumericMatrix XtCnX, const NumericMatrix CnytCnX, const NumericVector rscale, const IntegerVector gMap, const NumericVector gMapCounts, const NumericMatrix priorX, const double logDetPriorX, const int incCont);
double jzs_importance_marg_like2(const NumericVector q, const NumericVector mu, const NumericVector sig, const double sumSq, const int N, const NumericMatrix XtCnX, const NumericMatrix XtCny, const NumericVector rscale, const IntegerVector gMap, const NumericMatrix priorX, const double logDetPriorX, const int incCont);


// sign function 
// from http://stackoverflow.com/questions/1903954/is-there-a-standard-sign-function-signum-sgn-in-c-c
template <typename T> int sgn(T val) {
    return (T(0) < val) - (val < T(0));
}

double log1pExp(double x);
double logExpXplusExpY( const double x, const double y );
double logExpXminusExpY( const double x, const double y );
Eigen::MatrixXd random_multivariate_normal(const Eigen::MatrixXd mu, const Eigen::MatrixXd Sigma);


#endif //BFCOMMON_HPP_
