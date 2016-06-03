// [[Rcpp::interfaces(r, cpp)]]

#define _USE_MATH_DEFINES
#include <math.h>
#include "bfcommon.h"

using namespace Rcpp;

bool isgood( NumericVector x, double tol)
{ 
  int i=0;
  for( i=0 ; i < x.size() ; i++){
    if(x[i] != NA_REAL)
      if ( fabs(x[i]) > tol )
        return 0;
  }  
  return 1;
}

// [[Rcpp::export]]
NumericVector genhypergeo_series_pos( NumericVector U,
                                 NumericVector L,
                                 NumericVector z,
                                 const double tol,
                                 const int maxiter,
                                 const bool check_mod,
                                 const bool check_conds,
                                 const bool polynomial)
{
  
  NumericVector fac( z.size() );
  NumericVector temp( z.size() );
  NumericVector series( z.size() );
  LogicalVector greater( z.size() );
  int i=0,j=0;
  
  if(check_conds){
    if(is_true(any(U<0)) || is_true(any(L<0)) || is_true(any(z<0))){
      stop("All arguments must be positive.");
    }
  }

  if(check_mod){
    
    if( U.size() > (L.size()+1) ){
      greater = abs(z)>0;
    } else if(U.size() > L.size()) {
      greater = abs(z)>1;
    } else {
      greater  = abs(z)<0;
    }
    if( is_true( Rcpp::all(greater) ) ){
      return(z * NA_REAL);
    }else{
      for( i = 0 ; i < z.size() ; i++){
        if( greater[i] ) z[i] = NA_REAL;
      }
    } 
  }
  
  if(maxiter==0){
    return z*0+fac;
  }
  
  for ( i = 0; i < maxiter; i++ ) {

    fac = fac + sum(log( U + i) )  - sum( log( L + i ) ) +  log(z) - log( i + 1 );

    for( j = 0 ; j < z.size() ; j++ ){
      series[j] = logExpXplusExpY( temp[j], fac[j] );
    }
    
    if ( isgood( series - temp, tol ) ){
      return series;
    }
    temp = clone(series);
  }
  
  if(polynomial){
    return series;
  }else{
    Rcpp::warning("Series not converged.");
    return z * NA_REAL;  
  }
  
}


