#include <time.h>
#include "bfcommon.h"
#include "progress.h"


using namespace Rcpp;


double aFunc(const double rho, const int n, const double r, const bool hg_checkmod, const int hg_iter)
{
  NumericVector U(2, (n-1)*0.5);
  NumericVector L(1, .5);
  NumericVector z(1, r * r * rho * rho);

  double hyper_term = genhypergeo_series_pos(U, L, z, hg_checkmod, hg_iter, 0, 0, 0)[0];
  return ( ( n - 1 ) * 0.5 ) * log1p( -(rho*rho) ) + hyper_term;
}


double bFunc(const double rho, int const n, const double r, const bool hg_checkmod, const int hg_iter)
{
  NumericVector U(2, n*0.5);
  NumericVector L(1, 1.5);
  NumericVector z(1, r * r * rho * rho);

  double hyper_term = genhypergeo_series_pos(U, L, z, hg_checkmod, hg_iter, 0, 0, 0)[0];
  double log_term = 2 * ( lgamma(n*0.5) - lgamma((n-1)*0.5) ) +
    (n-1) * 0.5 * log1p( -(rho*rho) ) + log(2.0);
  return r * rho * exp( log_term + hyper_term );

}

// [[Rcpp::export]]
double hFunc(const double rho, const int n, const double r, const bool hg_checkmod, const int hg_iter)
{
  return log( exp( aFunc(rho, n, r, hg_checkmod, hg_iter) ) + bFunc(rho, n, r, hg_checkmod, hg_iter) );
}

// [[Rcpp::export]]
double jeffreys_approx_corr(const double rho, const int n, const double r)
{
  return 0.5*(n - 1) * log1p( -(rho*rho) ) - (n - 1 - 0.5)*log1p( -(rho*r) );
}

double corrtest_like_Rcpp(double zeta, NumericVector r, NumericVector n, double a_prior, double b_prior, bool approx, const bool hg_checkmod, const int hg_iter)
{
  int i;
  double rho = tanh( zeta );
  double logdens = Rf_dbeta( (rho+1.0)/2.0, a_prior, b_prior, 1) + log1p( -(rho*rho) );

  for( i = 0; i < r.size() ; i++ ){
    if( approx ){
      logdens +=  jeffreys_approx_corr( rho, n[i], r[i] );
    }else{
      logdens +=  hFunc( rho, n[i], r[i], hg_checkmod, hg_iter );
    }
  }
  return logdens;
}


// [[Rcpp::export]]
NumericMatrix metropCorrRcpp_jeffreys(NumericVector r, NumericVector n, double a_prior, double b_prior, bool approx, int iterations, bool doInterval,
                      NumericVector intervalz, bool intervalCompl, bool nullModel, int progress, Function callback, double callbackInterval)
{
    RNGScope scope;

    // setting last_cb to the beginning of the epoch
    // ensures that the callback is called once, first
    time_t last_cb = static_cast<time_t>(int(0));

    int i = 0;
    double Ubounds[2];
    NumericVector fish_z = 0.5*log( ( 1 + r ) / ( 1 - r ) );
    double candidate, z, trans_zeta;
    bool inInterval, valid_zeta = true;


    // For intervals
    if( doInterval){
      if( intervalz.size() == 0){
        doInterval = false;
      }else if( intervalz.size() != 2 ){
        Rcpp::stop("Incorrect number of interval points specified.");
      }
    }

    // starting values
    double fish_z0 = sum ( fish_z * ( n - 2 ) ) / sum( n - 2 );
    double fish_sd = sqrt( 1 / sum( n - 2 ) );
    double zeta = fish_z0;

    // create progress bar
    class Progress p(iterations, (bool) progress);

    // Create matrix for chains
    NumericMatrix chains(iterations, 2);

    if(nullModel){
      std::fill(chains.begin(), chains.end(), 0);
      return chains;
    }

    if(doInterval){
      Ubounds[0] = Rf_pnorm5( intervalz[0], fish_z0, fish_sd, 1, 0 );
      Ubounds[1] = Rf_pnorm5( intervalz[1], fish_z0, fish_sd, 1, 0 );
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
        candidate = Rf_qnorm5(candidate, fish_z0, fish_sd, 1, 0 );
      }else{
        candidate = Rf_rnorm( fish_z0, fish_sd );
      }

      // Metropolis-Hastings step
      z = corrtest_like_Rcpp(candidate, r, n, a_prior, b_prior, approx, 0, 2000) -
        corrtest_like_Rcpp(zeta, r, n, a_prior, b_prior, approx, 0, 2000) +
        Rf_dnorm4(zeta, fish_z0, fish_sd, 1) -
        Rf_dnorm4(candidate, fish_z0, fish_sd, 1);


      if(doInterval){
          trans_zeta = Rf_pnorm5(zeta, fish_z0, fish_sd, 1, 0 );
          inInterval = ( Ubounds[0] > trans_zeta ) && ( Ubounds[1] < trans_zeta );
          if( (inInterval && intervalCompl) || (!inInterval && !intervalCompl) )
            valid_zeta = false;
      }

      if( ( Rf_rexp(1) > -z ) || !valid_zeta ){
        zeta = candidate;
      }

      // copy to chains
      chains(i, 1) = zeta;
      chains(i, 0) = tanh(zeta);
    } // end sampler

    colnames(chains) = CharacterVector::create("rho", "zeta");
    return chains;
}



