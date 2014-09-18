#ifndef BFCOMMON_HPP_
#define BFCOMMON_HPP_

#include <Rcpp.h>

int RcppCallback(time_t *last, Rcpp::Function cb, double progress, double callbackInterval);

template <typename Type>
Type rowSums( Type X );

template <typename Type>
Type colSums( Type X );


#endif //BFCOMMON_HPP_
