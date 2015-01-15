#include "bfcommon.h"
#include "logRepresentedReal.h"

using namespace Rcpp;

// http://www.johndcook.com/blog/standard_deviation/
// [[Rcpp::export]]
NumericVector logSummaryStatsRcpp(NumericVector x) 
{
    int N = x.size(), i=0;
    logRepresentedReal oldM;
    logRepresentedReal M = logRepresentedReal(x[0], 1);
    logRepresentedReal S = logRepresentedReal(0, 0);
    logRepresentedReal thisX;
    
    NumericVector ret(2);
    ret(0) = x[0];
    ret(1) = R_NegInf;

    // Ensure that for N=0 and N=1 we return something meaningful 
    if( N == 0 ){
      ret(0) = NA_REAL;
      ret(1) = NA_REAL;
      return ret;
    }
        
    if(N == 1)
      return ret;
    
    for( i = 1; i < N; i++ ){
      oldM = M;
      thisX = logRepresentedReal( x[i], 1);
      M = oldM + ( thisX - oldM) / double(i + 1);   
      S = S + ( thisX - oldM ) * ( thisX - M );
    }
  
    ret(0) = M.modulo();
    ret(1) = S.modulo() - log(N - 1);
    return ret;
}
