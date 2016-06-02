#include "progress.h"
#include "bfcommon.h"
#include <time.h>

using namespace Rcpp;
using Eigen::MatrixXd;
using Eigen::Map;
using Eigen::Lower;

// [[Rcpp::depends(RcppEigen)]]

void jzs_mc_sampler(NumericVector *logsamples, const int iterations, const double sumSq, const int N, const NumericMatrix XtCnX, const NumericMatrix CnytCnX, const NumericVector rscale, const IntegerVector gMap, const NumericVector gMapCounts, const NumericMatrix priorX, const double logDetPriorX, const int incCont, const int progress, const Function callback, const double callbackInterval)
{
  RNGScope scope;

  // setting last_cb to the beginning of the epoch 
  // ensures that the callback is called once, first
  time_t last_cb = static_cast<time_t>(int(0));    
  
  const int nGs = gMapCounts.size();
  NumericVector g(nGs);

  int i = 0, j = 0, P = XtCnX.nrow();
  
  NumericMatrix XtCny( P, 1 );
  for( i = 0 ; i < P ; i++ )
    XtCny( i, 0 ) = CnytCnX( 0, i );
  
  // create progress bar
  class Progress p(iterations, (bool) progress);

  
  // Sampler
  for( i = 0 ; i < iterations ; i++)
  {
    
    // Check interrupt
    if (Progress::check_abort() )
      Rcpp::stop("Operation cancelled by interrupt.");
      
    p.increment(); // update progress
      
    // Check callback
    if( RcppCallback( &last_cb, callback, ( 1000.0 * ( i + 1 ) ) / iterations, callbackInterval) )
        Rcpp::stop("Operation cancelled by callback function.");

    // do sampling
    for( j = 0 ; j < nGs ; j++ ){
      g(j) = 1 / Rf_rgamma( 0.5, 2.0 / ( rscale(j) * rscale(j) ) );
    }

    (*logsamples)(i) = jzs_mc_marg_like2(g, sumSq, N, XtCnX, XtCny, gMap, priorX, logDetPriorX, incCont);
    //(*logsamples)(i) = jzs_mc_marg_like(g, sumSq, N, XtCnX, CnytCnX, rscale, gMap, gMapCounts, priorX, logDetPriorX, incCont);
  }
}

void jzs_importance_sampler(NumericVector *logsamples, const int iterations, const NumericVector mu, const NumericVector sig, const double sumSq, const int N, const NumericMatrix XtCnX, const NumericMatrix CnytCnX, const NumericVector rscale, const IntegerVector gMap, const NumericVector gMapCounts, const NumericMatrix priorX, const double logDetPriorX, const int incCont, const int progress, const Function callback, const double callbackInterval)
{
  RNGScope scope;  
  
  // setting last_cb to the beginning of the epoch 
  // ensures that the callback is called once, first
  time_t last_cb = static_cast<time_t>(int(0));    
  
  const int nGs = gMapCounts.size();
  NumericVector q(nGs);

  int i = 0, j = 0, P = XtCnX.nrow();
  
  NumericMatrix XtCny( P, 1 );
  for( i = 0 ; i < P ; i++ )
    XtCny( i, 0 ) = CnytCnX( 0, i );
  
  // create progress bar
  class Progress p(iterations, (bool) progress);

  
  // Sampler
  for( i = 0 ; i < iterations ; i++)
  {
    
    // Check interrupt
    if (Progress::check_abort() )
      Rcpp::stop("Operation cancelled by interrupt.");
      
    p.increment(); // update progress
      
    // Check callback
    if( RcppCallback( &last_cb, callback, ( 1000.0 * ( i + 1 ) ) / iterations, callbackInterval) )
        Rcpp::stop("Operation cancelled by callback function.");

    // do sampling
    for( j = 0 ; j < nGs ; j++ ){
      q(j) = Rf_rnorm( mu(j), sig(j) );
    }

    (*logsamples)(i) = jzs_importance_marg_like2(q, mu, sig, sumSq, N, XtCnX, XtCny, rscale, gMap, priorX, logDetPriorX, incCont);
    //(*logsamples)(i) = jzs_importance_marg_like(q, mu, sig, sumSq, N, XtCnX, CnytCnX, rscale, gMap, gMapCounts, priorX, logDetPriorX, incCont);
  }
}

// [[Rcpp::export]]
NumericVector jzs_sampler(const int iterations, const NumericVector y, const NumericMatrix X, const NumericVector rscale, const IntegerVector gMap, const int incCont, const NumericVector importanceMu, const NumericVector importanceSig, const int progress, const Function callback, const double callbackInterval, const int which)
{
  // which = 0 for mc sampler 
  // which = 1 for importance sampler
  int i = 0, j = 0;
  const int N = X.nrow();
  const int P = X.ncol();
  const double ybar = mean(y);
  double sumSq = 0, logDetPriorX = 0;
  NumericVector logsamples( iterations );
  // gMapCounts is not needed for the sampler, but we need to pass something.
  NumericVector gMapCounts( max(gMap) + 1 );
  NumericMatrix priorX(incCont, incCont);
  NumericMatrix CnX(N, P);
  NumericMatrix CnytCnX(1, P);
  NumericVector Cny(N);
  NumericVector XcolMeans(P);
  
  // Compute mean of each column of design matrix
  for( i = 0 ; i < P ; i++ )
  {
    XcolMeans(i) = mean( X( _, i ) );
  }
  
  // Compute centered matrices
  for( i = 0 ; i < N ; i++ )
  {
      Cny(i) = y(i) - ybar;
      sumSq += Cny(i) * Cny(i);
      for( j = 0 ; j < P ; j++ )
      {
        CnX(i,j) = X(i,j) - XcolMeans(j);  
      }
  }
  
  MatrixXd XtCnX(MatrixXd(P,P).setZero().selfadjointView<Lower>().rankUpdate( MatrixXd((as<Map<MatrixXd> >(CnX))).transpose()));
  
  // Construct prior cov matrix for continuous covariates from X
  if(incCont){
    for( i = 0 ; i < incCont ; i++ ){
      for( j = 0 ; j <= i ; j++ ){
        priorX(i,j) = sum( CnX( _ , i) * CnX( _ , j) ) / N;
        priorX(j,i) = priorX(i,j);
      }
    }
    logDetPriorX = log_determinant_pos_def(MatrixXd((as<Map<MatrixXd> >(priorX))));
  }
  
  // Compute t(Cy) %*% CX
  CnytCnX = wrap(MatrixXd((as<Map<MatrixXd> >(Cny))).transpose() * MatrixXd((as<Map<MatrixXd> >(CnX))));
  
  if( which == 0 ){ // mc sampler
    jzs_mc_sampler(&logsamples, iterations, sumSq, N, wrap(XtCnX), CnytCnX, rscale, gMap, gMapCounts, priorX, logDetPriorX, incCont, progress, callback, callbackInterval);
  }else if( which == 1 ){ // importance sampler
    jzs_importance_sampler(&logsamples, iterations, importanceMu, importanceSig, sumSq, N, wrap(XtCnX), CnytCnX, rscale, gMap, gMapCounts, priorX, logDetPriorX, incCont, progress, callback, callbackInterval);  
  }else{
    Rcpp::stop("Invalid sampler specified.");
  }
  return logsamples;
}


