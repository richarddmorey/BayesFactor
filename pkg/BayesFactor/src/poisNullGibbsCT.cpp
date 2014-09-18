// [[Rcpp::depends(RcppProgress)]]
//#include <progress.hpp>
#include <time.h>
#include "bfcommon.h"

/* NOTE: RcppProgress code is disabled until
   I can fix the header issue. */

using namespace Rcpp;

// [[Rcpp::export]]
NumericMatrix gibbsContTabPoissonNull(IntegerMatrix y, double a, double b, int iterations, int progress, Function callback, double callbackInterval) 
{
    // setting last_cb to the beginning of the epoch 
    // ensures that the callback is called once, first
    time_t last_cb = static_cast<time_t>(int(0));    
  
    double lambda_a = 0, lambda_b = 0;
    int k = 0, i = 0, j = 0;
    int I = y.nrow(), J = y.ncol(), sumY = sum(y);

    // starting values
    double lambda = sumY / ( I * J );
    NumericVector p(I), q(J);
    
    // create progress bar
    //Progress p(iterations, (bool) progress);

    // Create matrix for chains
    NumericMatrix chains(iterations, 1 + I + J + I*J);
    
    // Start Gibbs sampler
    for( k = 0 ; k < iterations ; k++ )
    {

      // Check interrupt
      //if (Progress::check_abort() )
      //  Rcpp::stop("Operation cancelled by interrupt.");
      
      //p.increment(); // update progress
      
      // Check callback
      if( RcppCallback( &last_cb, callback, ( 1000.0 * ( k + 1 ) ) / iterations, callbackInterval) )
        Rcpp::stop("Operation cancelled by callback function.");


      // sample lambda
      lambda_a = sumY + I*J*(a - 1) - I - J + 1;
      lambda_b = b;
      for( i = 0 ; i < I ; i++)
        for( j = 0 ; j < J ; j++ )
          lambda_b += p[i] * q[j];
      lambda = Rf_rgamma( lambda_a, 1 / lambda_b );      

      // sample p_i
      
  	  
      // sample q_j
      
      
      
      // copy to chains
      chains(k, 0) = lambda;
    } // end Gibbs sampler

    return chains;
}
