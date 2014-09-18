#include <Rcpp.h>
#include "bfcommon.h"

using namespace Rcpp;

template <typename Type>
Type colSums( Type X ) {
   
   int I = X.nrow(), i = 0;
   int J = X.ncol(), j = 0;
   Type sums( 1, J );
   
   
   for( j = 0; j < J; j++ )
    for( i = 0; i < I; i++)
      sums[j] += X[i,j];
    
  return sums;  
}

template <typename Type>
Type rowSums( Type X ) {
   
   int I = X.nrow(), i = 0;
   int J = X.ncol(), j = 0;
   Type sums( I, 1 );
   
   
   for( j = 0; j < J; j++ )
    for( i = 0; i < I; i++)
      sums[i] += X[i,j];
    
  return sums;  
}
