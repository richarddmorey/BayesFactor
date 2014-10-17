#include "progress.h"
#include <time.h>
#include "bfcommon.h"

double meta_t_like_Rcpp(double delta, Rcpp::NumericVector t, Rcpp::NumericVector n, Rcpp::NumericVector df, double rscale);

using namespace Rcpp;

// [[Rcpp::export]]
NumericMatrix metropMetaTRcpp(NumericVector t, NumericVector n1, NumericVector n2, bool twoSample, double rscale, int iterations, bool doInterval, 
                      NumericVector interval, bool intervalCompl, bool nullModel, int progress, Function callback, double callbackInterval) 
{
    // setting last_cb to the beginning of the epoch 
    // ensures that the callback is called once, first
    time_t last_cb = static_cast<time_t>(int(0));    

    int i = 0;
    double Ubounds[2];
    double candidate, z, transDelta;
    bool inInterval, validDelta = true;
    
    NumericVector d(clone(t));
    NumericVector eff_n(clone(n1));
    NumericVector nu(clone(n1));
  
    // For intervals
    if( doInterval){
      if( interval.size() == 0){
        doInterval = false;
      }else if( interval.size() != 2 ){
        Rcpp::stop("Incorrect number of interval points specified.");
      }
    }

    // effective sample size and degrees of freedom
    if(twoSample){
      eff_n = n1 * n2  / (n1 + n2);
      nu = n1 + n2 - 2.0;
    }else{
      eff_n = n1;
      nu = n1 - 1.0;
    }

    // starting values
    for( i = 0 ; i < t.size() ; i++ ){
      d[i] = d[i] / sqrt( eff_n[i] );
    }
    double delta0 = sum( d * eff_n ) / sum( eff_n ) ;
    double delta_sd = 1 / sqrt( sum( eff_n ) );
    double delta = delta0;
    
    // create progress bar
    Progress::Progress p(iterations, (bool) progress);

    // Create matrix for chains
    NumericMatrix chains(iterations, 1);
    
    if(nullModel)
      return chains;
     
    if(doInterval){
      Ubounds[0] = Rf_pnorm5( interval[0], delta0, delta_sd, 1, 0 );
      Ubounds[1] = Rf_pnorm5( interval[1], delta0, delta_sd, 1, 0 );
    }
 
    // Start sampler
    for( i = 0 ; i < iterations ; i++ )
    {

      // Check interrupt
      if (Progress::check_abort() )
        Rcpp::stop("Operation cancelled by interrupt.");
      
      p.increment(); // update progress
      
      // Check callback
      if( RcppCallback( &last_cb, callback, ( 1000.0 * ( i + 1 ) ) / iterations, callbackInterval) )
        Rcpp::stop("Operation cancelled by callback function.");

      // sample delta      
      if(doInterval){
        if(intervalCompl){
          candidate = Rf_runif(0, Ubounds[0] + 1 - Ubounds[1]);
          if( candidate > Ubounds[0]) 
            candidate = candidate - Ubounds[0] + Ubounds[1];
        }else{
          candidate = Rf_runif(Ubounds[0], Ubounds[1]);
        }
        candidate = Rf_qnorm5(candidate, delta0, delta_sd, 1, 0 );
      }else{
        candidate = Rf_rnorm( delta0, delta_sd );
      }
      
      // Metropolis-Hastings step
      z = meta_t_like_Rcpp(candidate, t, eff_n, nu, rscale) - 
        meta_t_like_Rcpp(delta, t, eff_n, nu, rscale) +
        Rf_dnorm4(delta, delta0, delta_sd, 1) - 
        Rf_dnorm4(candidate, delta0, delta_sd, 1);
      
      
      if(doInterval){
          transDelta = Rf_pnorm5(delta, delta0, delta_sd, 1, 0 );
          inInterval = ( Ubounds[0] > transDelta ) && ( Ubounds[1] < transDelta );
          if( (inInterval && intervalCompl) || (!inInterval && !intervalCompl) )
            validDelta = false;
      }
      
      if( ( Rf_rexp(1) > -z ) || !validDelta ){
        delta = candidate; 
      }

      // copy to chains
      chains(i, 0) = delta;
    } // end sampler

    return chains;
}

double meta_t_like_Rcpp(double delta, NumericVector t, NumericVector n, NumericVector df, double rscale)
{
  int i;
  double logdens = Rf_dcauchy( delta, 0, rscale, 1);
 
  for( i = 0; i < t.size() ; i++ ){
    logdens +=  Rf_dnt( t[i], df[i], delta*sqrt(n[i]), 1);
  }
  return logdens;
}


