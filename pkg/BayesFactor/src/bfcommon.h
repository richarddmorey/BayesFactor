#ifndef BFCOMMON_HPP_
#define BFCOMMON_HPP_

#include <Rcpp.h>

int RcppCallback(time_t *last, Rcpp::Function cb, double progress, double callbackInterval);

template <typename Type>
Type rowSums( Type X );

template <typename Type>
Type colSums( Type X );


// sign function 
// from http://stackoverflow.com/questions/1903954/is-there-a-standard-sign-function-signum-sgn-in-c-c
template <typename T> int sgn(T val) {
    return (T(0) < val) - (val < T(0));
}

double log1pExp(double x);
double logExpXplusExpY( const double x, const double y );
double logExpXminusExpY( const double x, const double y );


#endif //BFCOMMON_HPP_
