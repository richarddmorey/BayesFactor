#include "progress.h"
#include <time.h>
#include "bfcommon.h"

double proptest_like_Rcpp(double lo, Rcpp::NumericVector y, Rcpp::NumericVector n, double p, double rscale);

using namespace Rcpp;

// [[Rcpp::export]]
NumericMatrix metropProportionRcpp(NumericVector y, NumericVector n, double p0, double rscale, int iterations, bool doInterval,
                      NumericVector interval, bool intervalCompl, bool nullModel, int progress, Function callback, double callbackInterval)
{
    RNGScope scope;

    // setting last_cb to the beginning of the epoch
    // ensures that the callback is called once, first
    time_t last_cb = static_cast<time_t>(int(0));

    int i = 0;
    double Ubounds[2], mu = log( p0 / ( 1 - p0 ) );
    double candidate, z, trans_lo;
    bool inInterval, valid_lo = true;


    // For intervals
    if( doInterval){
      if( interval.size() == 0){
        doInterval = false;
      }else if( interval.size() != 2 ){
        Rcpp::stop("Incorrect number of interval points specified.");
      }
    }

    // starting values
    double ln_p0 = log( sum(y) + 1.0 ) - log( sum(n) + 2.0 );
    double lo0 = Rf_qlogis( ln_p0, 0, 1, 1, 1 );
    double lo_sd = exp ( 0.5 * ( ln_p0 + Rf_pexp( -ln_p0, 1, 1, 1 )  -
                                 log( sum(n) + 2.0 ) - 2 * Rf_dlogis( lo0, 0, 1, 1 ) ) );
    double lo = lo0;

    // create progress bar
    class Progress p(iterations, (bool) progress);

    // Create matrix for chains
    NumericMatrix chains(iterations, 1);

    if(nullModel){
      std::fill(chains.begin(), chains.end(), mu);
      return chains;
    }

    if(doInterval){
      Ubounds[0] = Rf_pnorm5( interval[0], lo0, lo_sd, 1, 0 );
      Ubounds[1] = Rf_pnorm5( interval[1], lo0, lo_sd, 1, 0 );
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
        candidate = Rf_qnorm5(candidate, lo0, lo_sd, 1, 0 );
      }else{
        candidate = Rf_rnorm( lo0, lo_sd );
      }

      // Metropolis-Hastings step
      z = proptest_like_Rcpp(candidate, y, n, mu, rscale) -
        proptest_like_Rcpp(lo, y, n, mu, rscale) +
        Rf_dnorm4(lo, lo0, lo_sd, 1) -
        Rf_dnorm4(candidate, lo0, lo_sd, 1);


      if(doInterval){
          trans_lo = Rf_pnorm5(lo, lo0, lo_sd, 1, 0 );
          inInterval = ( Ubounds[0] > trans_lo ) && ( Ubounds[1] < trans_lo );
          if( (inInterval && intervalCompl) || (!inInterval && !intervalCompl) )
            valid_lo = false;
      }

      if( ( Rf_rexp(1) > -z ) || !valid_lo ){
        lo = candidate;
      }

      // copy to chains
      chains(i, 0) = lo;
    } // end sampler

    return chains;
}


double proptest_like_Rcpp(double lo, NumericVector y, NumericVector n, double mu, double rscale)
{
  int i;
  double theta = 1 / ( 1 + exp(-lo) );
  double logdens = Rf_dlogis( lo, mu, rscale, 1);

  for( i = 0; i < y.size() ; i++ ){
    logdens +=  Rf_dbinom( y[i], n[i], theta, 1);
  }
  return logdens;
}


