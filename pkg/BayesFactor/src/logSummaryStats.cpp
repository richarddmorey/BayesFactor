#include "bfcommon.h"
#include "logRepresentedReal.h"

using namespace Rcpp;

// Welford's method
// http://www.johndcook.com/blog/standard_deviation/
// [[Rcpp::export]]
List logSummaryStats(NumericVector x) 
{
    int N = x.size(), i=0;
    logRepresentedReal oldM;
    logRepresentedReal M = logRepresentedReal(x[0], 1);
    logRepresentedReal S = logRepresentedReal(0, 0);
    logRepresentedReal thisX;
    
    NumericVector retM(1);
    NumericVector retS(1);
    NumericVector retCumu( N == 0 ? 1 : N, NA_REAL);
    List ret;
    
    ret = List::create(Rcpp::Named("logMean") = retM, 
                 Rcpp::Named("logVar") =  retS, 
                 Rcpp::Named("cumLogMean") = retCumu);
                 
    // Ensure that for N=0 and N=1 we return something meaningful 
    if( N == 0 ){
      retM(0) = NA_REAL;
      retS(0) = NA_REAL;
      retCumu(0) = NA_REAL;
      return ret;
    }
    
    retM(0) = x[0];
    retS(0) = R_NegInf;
    retCumu(0) = x[0];
    
    if(N == 1)
      return ret;
    
    for( i = 1; i < N; i++ ){
      oldM = M;
      thisX = logRepresentedReal( x[i], 1);
      M = oldM + ( thisX - oldM) / double(i + 1);   
      S = S + ( thisX - oldM ) * ( thisX - M );
      
      retCumu(i) = M.modulo();
    }
  
    retM(0) = M.modulo();
    retS(0) = S.modulo() - log(N - 1);
    return ret;
}
