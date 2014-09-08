#include <Rcpp.h>
#include <time.h>
#include "bfcommon.hpp"

using namespace Rcpp;

int RcppCallback(time_t *last, Rcpp::Function cb, double progress, double callbackInterval)
{
   
  IntegerVector callbackResult(1);
   
  time_t now = time(NULL);

  if( difftime( now , *last ) > callbackInterval ){
    callbackResult = cb( progress );
    *last = now;
    return callbackResult[0];
  }else{
    return (int) 0;
  }
}
