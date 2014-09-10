// [[Rcpp::depends(RcppProgress)]]
//#include <progress.hpp>
#include <time.h>
#include "bfcommon.h"

/* NOTE: RcppProgress code is disabled until
   I can fix the header issue. */

using namespace Rcpp;

// [[Rcpp::export]]
NumericMatrix gibbsOneSampleRcpp(double ybar, double s2, int N, double rscale, int iterations, bool doInterval, 
                      NumericVector interval, bool intervalCompl, bool nullModel, int progress, Function callback, double callbackInterval) 
{
    // setting last_cb to the beginning of the epoch 
    // ensures that the callback is called once, first
    time_t last_cb = static_cast<time_t>(int(0));    

    int i = 0, whichInterval = 0, signAgree = 1;
    double meanMu, varMu, scaleSig2, scaleg;
    double shapeSig2 = 0.5 * N  + 0.5;
    double sumy2 = (N - 1) * s2 + N * pow(ybar, 2);
    double rscaleSq = pow(rscale, 2);
    double intLower = 0, intUpper = 1, areaLower, areaUpper;
    

    // For intervals
    if( doInterval){
      signAgree = (interval[0] * interval[1]) >= 0;
      if( interval.size() == 0){
        doInterval = false;
      }else if( interval.size() != 2 ){
        Rcpp::stop("Incorrect number of interval points specified.");
      }
    }

    // starting values
    double mu = ybar, sig2 = s2, g = pow(ybar, 2) / s2 + 1;
    
    if(nullModel) mu = 0;
    
    // create progress bar
    //Progress p(iterations, (bool) progress);

    // Create matrix for chains
    NumericMatrix chains(iterations, 4);
    
    // Start Gibbs sampler
    for( i = 0 ; i < iterations ; i++ )
    {

      // Check interrupt
      //if (Progress::check_abort() )
      //  Rcpp::stop("Operation cancelled by interrupt.");
      
      //p.increment(); // update progress
      
      // Check callback
      if( RcppCallback( &last_cb, callback, ( 1000.0 * ( i + 1 ) ) / iterations, callbackInterval) )
        Rcpp::stop("Operation cancelled by callback function.");


      // sample mu
  	  varMu  = sig2 / ( 1.0 * N + 1/g );
      meanMu = ybar * N * varMu / sig2;
      
      if(doInterval && !nullModel){
        if( !intervalCompl ){
          // Interval as given
          intLower = Rf_pnorm5( sqrt(sig2) * interval[0], meanMu, sqrt(varMu), 1, 0 );
          intUpper = Rf_pnorm5( sqrt(sig2) * interval[1], meanMu, sqrt(varMu), 1, 0 );
        }else{
          // Complement of interval
          // Compute area of both sides and choose one
          areaLower = Rf_pnorm5( sqrt(sig2) * interval[0], meanMu, sqrt(varMu), 1, 1 );
          areaUpper = Rf_pnorm5( sqrt(sig2) * interval[1], meanMu, sqrt(varMu), 0, 1 ); 
          whichInterval = Rf_rlogis( areaUpper - areaLower, 1 ) > 0;           
          // Sample from chosen side
          if(whichInterval){
            intLower = Rf_pnorm5( sqrt(sig2) * interval[1], meanMu, sqrt(varMu), 1, 0 );
            intUpper = 1;
          }else{
            intLower = 0;
            intUpper = Rf_pnorm5( sqrt(sig2) * interval[0], meanMu, sqrt(varMu), 1, 0 );
          }
         } 
        mu = Rf_runif(intLower, intUpper);
        mu = Rf_qnorm5( mu, meanMu, sqrt(varMu), 1, 0 );
      }else{
        // no interval
        if(nullModel){
          mu = 0; 
        }else{
          mu = Rf_rnorm( meanMu, sqrt(varMu) );
        }
      }
      
      // sample sig2
		  scaleSig2 = 0.5 * ( sumy2 - 2.0 * N * ybar * mu + (N + 1/g)*pow(mu,2) );
      if(doInterval){
        if( !intervalCompl){
          // Interval as given
          if( signAgree ){
            // signs of endpoints of interval agree - lower and upper bound
            intLower = Rf_pgamma( pow( interval[0] / mu, 2), shapeSig2, 1/scaleSig2, 1, 0 );
            intUpper = Rf_pgamma( pow( interval[1] / mu, 2), shapeSig2, 1/scaleSig2, 1, 0 );
          }else{
           // signs of endpoints of interval do not agree - no lower bound
           intLower = 0;
            if( mu >= 0 ){
              intUpper = Rf_pgamma( pow( interval[1] / mu, 2), shapeSig2, 1/scaleSig2, 1, 0 );
            }else{
              intUpper = Rf_pgamma( pow( interval[0] / mu, 2), shapeSig2, 1/scaleSig2, 1, 0 );
            }
          }          
        }else{
          // Complement of interval
          if( signAgree ){
            // Signs of interval end points agree 
            if( (mu * interval[0]) < 0){
                // Unrestricted sampling 
                intLower = 0;
                intUpper = 1;              
            }else{
              // Compute area of both sides and choose one
              areaLower = Rf_pgamma( pow( interval[0] / mu, 2), shapeSig2, 1/scaleSig2, 1, 1 );
              areaUpper = Rf_pgamma( pow( interval[1] / mu, 2), shapeSig2, 1/scaleSig2, 0, 1 );
              whichInterval = Rf_rlogis( areaUpper - areaLower, 1 ) > 0;           
              // Sample from chosen side
              if(whichInterval){
                intLower = Rf_pgamma( pow( interval[1] / mu, 2), shapeSig2, 1/scaleSig2, 1, 0 );
                intUpper = 1;
              }else{
                intLower = 0;
                intUpper = Rf_pgamma( pow( interval[0] / mu, 2), shapeSig2, 1/scaleSig2, 1, 0 );
              }
            }
          }else{
            // signs of endpoints of interval do not agree - no upper bound
            intUpper = 1;
            if(mu >= 0){
              intLower = Rf_pgamma( pow( interval[1] / mu, 2), shapeSig2, 1/scaleSig2, 1, 0 );
            }else{
              intLower = Rf_pgamma( pow( interval[0] / mu, 2), shapeSig2, 1/scaleSig2, 1, 0 );   
            }
          } 
        } // end doInterval
        sig2 = Rf_runif(intLower, intUpper);
        sig2 = 1 / Rf_qgamma( sig2, shapeSig2, 1/scaleSig2, 1, 0 );
      }else{
        // No interval
        sig2 = 1 / Rf_rgamma( shapeSig2, 1/scaleSig2 );
      }
      
  	  // sample g
		  scaleg = 0.5 * ( pow(mu,2) / sig2 + rscaleSq );
		  g = 1 / Rf_rgamma( 1, 1/scaleg );
      
      // copy to chains
      chains(i, 0) = mu;
      chains(i, 1) = sig2;
      chains(i, 2) = mu / sqrt( sig2 );
      chains(i, 3) = g;
    } // end Gibbs sampler

    return chains;
}
