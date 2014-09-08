// [[Rcpp::depends(RcppProgress)]]
#include <progress.hpp>
#include <time.h>
#include "bfcommon.h"

using namespace Rcpp;

// [[Rcpp::export]]
NumericMatrix gibbsTwoSampleRcpp(NumericVector ybar, NumericVector s2, NumericVector N, double rscale, int iterations, bool doInterval, 
                      NumericVector interval, bool intervalCompl, bool nullModel, int progress, Function callback, double callbackInterval) 
{
    
    // setting last_cb to the beginning of the epoch 
    // ensures that the callback is called once, first
    time_t last_cb = static_cast<time_t>(int(0));
    
    int i = 0, whichInterval = 0, signAgree = 1;
    double meanMu, varMu, meanBeta, varBeta, scaleSig2, scaleg;
    double shapeSig2 = 0.5 * sum(N) + 0.5;
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
    double sumy1 = N[0] * ybar[0];
    double sumy2 = N[1] * ybar[1];
    double sumy1Sq = (N[0] - 1) * s2[0] + N[0] * pow(ybar[0], 2);
    double sumy2Sq = (N[1] - 1) * s2[1] + N[1] * pow(ybar[1], 2);
    double sumySq = sumy1Sq + sumy2Sq; 
    double sumy = sumy1 + sumy2; 
    double diffy = N[1] * ybar[1] - N[0] * ybar[0];
    double sumN = sum(N) * 1.0;
    double diffN = ( N[0] - N[1] ) * 1.0;
    
    double mu = sumy / sum(N);
    double beta = ybar[1] - ybar[0];
    double sig2 = ( s2[0] * (N[0] - 1) + s2[1] * (N[1] - 1) ) / ( sum(N) - 2 );
    double g = pow(beta, 2) / sig2 + 1;
    
    if(nullModel) beta = 0;

    
    // create progress bar
    Progress p(iterations, (bool) progress);

    // Create matrix for chains
    NumericMatrix chains(iterations, 5);
    
    // Start Gibbs sampler
    for( i = 0 ; i < iterations ; i++ )
    {

      // Check interrupt
      if (Progress::check_abort() )
        Rcpp::stop("Operation cancelled by interrupt.");
      
      p.increment(); // update progress
      
      // Check callback
      if( RcppCallback( &last_cb, callback, ( 1000.0 * ( i + 1 ) ) / iterations, callbackInterval) )
        Rcpp::stop("Operation cancelled by callback function.");


      // sample mu
      meanMu = (sumy + ( N[0] - N[1] ) * beta / 2) / sumN ;
      varMu = sig2 / sumN;
      mu = Rf_rnorm( meanMu, sqrt(varMu) );
      
      // sample beta
  	  varBeta  = sig2 / ( sumN/4 + 1/g );
      meanBeta = varBeta / sig2 * ( (sumy2 - sumy1) + mu * ( diffN ) ) / 2;
      
      if(doInterval && !nullModel){
        if( !intervalCompl ){
          // Interval as given
          intLower = Rf_pnorm5( sqrt(sig2) * interval[0], meanBeta, sqrt(varBeta), 1, 0 );
          intUpper = Rf_pnorm5( sqrt(sig2) * interval[1], meanBeta, sqrt(varBeta), 1, 0 );
        }else{
          // Complement of interval
          // Compute area of both sides and choose one
          areaLower = Rf_pnorm5( sqrt(sig2) * interval[0], meanBeta, sqrt(varBeta), 1, 1 );
          areaUpper = Rf_pnorm5( sqrt(sig2) * interval[1], meanBeta, sqrt(varBeta), 0, 1 ); 
          whichInterval = Rf_rlogis( areaUpper - areaLower, 1 ) > 0;           
          // Sample from chosen side
          if(whichInterval){
            intLower = Rf_pnorm5( sqrt(sig2) * interval[1], meanBeta, sqrt(varBeta), 1, 0 );
            intUpper = 1;
          }else{
            intLower = 0;
            intUpper = Rf_pnorm5( sqrt(sig2) * interval[0], meanBeta, sqrt(varBeta), 1, 0 );
          }
         } 
        beta = Rf_runif(intLower, intUpper);
        beta = Rf_qnorm5( beta, meanBeta, sqrt(varBeta), 1, 0 );
      }else{
        // no interval
        if(nullModel){
          beta = 0;
        }else{
          beta = Rf_rnorm( meanBeta, sqrt(varBeta) );
        }
      } // end sample beta
      
      
      // sample sig2
		  scaleSig2 = 0.5 * ( sumySq - 2.0 * mu * sumy - beta * diffy +
                          N[0] * pow(mu - beta/2, 2) + N[1] * pow(mu + beta/2, 2) +
                          pow(beta, 2) / g );
      if(doInterval){
        if( !intervalCompl){
          // Interval as given
          if( signAgree ){
            // signs of endpoints of interval agree - lower and upper bound
            intLower = Rf_pgamma( pow( interval[0] / beta, 2), shapeSig2, 1/scaleSig2, 1, 0 );
            intUpper = Rf_pgamma( pow( interval[1] / beta, 2), shapeSig2, 1/scaleSig2, 1, 0 );
          }else{
           // signs of endpoints of interval do not agree - no lower bound
           intLower = 0;
            if( beta >= 0 ){
              intUpper = Rf_pgamma( pow( interval[1] / beta, 2), shapeSig2, 1/scaleSig2, 1, 0 );
            }else{
              intUpper = Rf_pgamma( pow( interval[0] / beta, 2), shapeSig2, 1/scaleSig2, 1, 0 );
            }
          }          
        }else{
          // Complement of interval
          if( signAgree ){
            // Signs of interval end points agree 
            if( (beta * interval[0]) < 0){
                // Unrestricted sampling 
                intLower = 0;
                intUpper = 1;              
            }else{
              // Compute area of both sides and choose one
              areaLower = Rf_pgamma( pow( interval[0] / beta, 2), shapeSig2, 1/scaleSig2, 1, 1 );
              areaUpper = Rf_pgamma( pow( interval[1] / beta, 2), shapeSig2, 1/scaleSig2, 0, 1 );
              whichInterval = Rf_rlogis( areaUpper - areaLower, 1 ) > 0;           
              // Sample from chosen side
              if(whichInterval){
                intLower = Rf_pgamma( pow( interval[1] / beta, 2), shapeSig2, 1/scaleSig2, 1, 0 );
                intUpper = 1;
              }else{
                intLower = 0;
                intUpper = Rf_pgamma( pow( interval[0] / beta, 2), shapeSig2, 1/scaleSig2, 1, 0 );
              }
            }
          }else{
            // signs of endpoints of interval do not agree - no upper bound
            intUpper = 1;
            if(beta >= 0){
              intLower = Rf_pgamma( pow( interval[1] / beta, 2), shapeSig2, 1/scaleSig2, 1, 0 );
            }else{
              intLower = Rf_pgamma( pow( interval[0] / beta, 2), shapeSig2, 1/scaleSig2, 1, 0 );   
            }
          } 
        } // end doInterval
        sig2 = Rf_runif(intLower, intUpper);
        sig2 = 1 / Rf_qgamma( sig2, shapeSig2, 1/scaleSig2, 1, 0 );
      }else{
        // No interval
        sig2 = 1 / Rf_rgamma( shapeSig2, 1/scaleSig2 );
      } // end sample sig2
      
  	  // sample g
		  scaleg = 0.5 * ( pow(beta,2) / sig2 + rscaleSq );
		  g = 1 / Rf_rgamma( 1, 1/scaleg );
      
      // copy to chains
      chains(i, 0) = mu;
      chains(i, 1) = beta;
      chains(i, 2) = sig2;
      chains(i, 3) = beta / sqrt( sig2 );
      chains(i, 4) = g;

    } // end Gibbs sampler

    return chains;
}
