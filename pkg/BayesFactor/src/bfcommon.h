#ifndef BFCOMMON_HPP_
#define BFCOMMON_HPP_

#include <Rcpp.h>

int RcppCallback(time_t *last, Rcpp::Function cb, double progress, double callbackInterval);

#endif //BFCOMMON_HPP_
