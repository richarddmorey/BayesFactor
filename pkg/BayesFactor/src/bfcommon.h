#ifndef BFCOMMON_HPP_
#define BFCOMMON_HPP_

#include <RcppEigen.h>

int RcppCallback(time_t *last, Rcpp::Function cb, double progress, double callbackInterval);

double dinvgamma1_Rcpp(const double x, const double a, const double b);
double ddinvgamma1_Rcpp(const double x, const double a, const double b);
double d2dinvgamma1_Rcpp(const double x, const double a, const double b);

double log_determinant_pos_def(Eigen::MatrixXd A);


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
