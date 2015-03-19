#include "progress.h"
#include "bfcommon.h"
#include <RcppEigen.h>
#include <time.h>

using namespace Rcpp;
using Eigen::MatrixXd;
using Eigen::Map;
using Eigen::Lower;

// [[Rcpp::depends(RcppEigen)]]

// [[Rcpp::export]]
NumericMatrix jzs_Gibbs(const int iterations, const NumericVector y, const NumericMatrix X, const NumericVector rscale, const double sig2start, const IntegerVector gMap, const NumericVector gMapCounts, const int incCont, bool nullModel, const IntegerVector ignoreCols, const int thin, const int progress, const Function callback, const double callbackInterval)
{
  RNGScope scope;

  // setting last_cb to the beginning of the epoch 
  // ensures that the callback is called once, first
  time_t last_cb = static_cast<time_t>(int(0));    

  const int N = X.nrow();
  const int P = X.ncol() - 1;
  
  if(P != gMap.size())
    Rcpp::stop("Design matrix size / gMap mismatch; ncol(X) should be length(gMap)+1");
  
  const int nGs = gMapCounts.size(); 
  const int nPars = P + 2 + nGs;
  const int nOutputPars = nPars - sum(ignoreCols);
  const int effectiveIterations = iterations / thin;
  const double ybar = mean(y);
  int i = 0, j = 0, k = 0, colCounter = 0, rowCounter = 0;

  const Map<MatrixXd> Xm(as<Map<MatrixXd> >(X));
  const Map<MatrixXd> ym(as<Map<MatrixXd> >(y));
  const MatrixXd XtX( MatrixXd( P + 1, P + 1 ).setZero().selfadjointView<Lower>().rankUpdate( Xm.transpose() ) );
  NumericMatrix chains( effectiveIterations, nOutputPars );

  NumericVector g(nGs, 1.0);
  NumericVector SSqG(nGs, 0.0);
  double sig2 = sig2start, SSq = 0;
  MatrixXd beta( MatrixXd( P + 1, 1 ).setZero() );
  MatrixXd yResid( MatrixXd( N, 1 ).setZero() );  
  MatrixXd Sigma( XtX.selfadjointView<Lower>() );
  MatrixXd Xty( Xm.transpose() * ym );  
  MatrixXd mu( MatrixXd( P + 1, 1 ).setZero() );

  // create progress bar
  class Progress p(iterations, (bool) progress);
  
  for( i = 0 ; i < iterations ; i++ ){ //Gibbs sampler
    
    // Check interrupt
    if (Progress::check_abort() )
      Rcpp::stop("Operation cancelled by interrupt.");
      
    p.increment(); // update progress
      
    // Check callback
    if( RcppCallback( &last_cb, callback, ( 1000.0 * ( i + 1 ) ) / iterations, callbackInterval) )
        Rcpp::stop("Operation cancelled by callback function.");
    
    
    // sample beta
    if(nullModel){ // all beta except grand mean are 0
      beta( 0, 0 ) = Rf_rnorm( ybar, sqrt( sig2 / N ) );
      
      for( j = 1 ; j < N ; j++ )
        yResid( j, 0 ) = ym( j, 0 ) - beta( 0, 0 );
    }else{ // null is false
      
      Sigma = XtX.selfadjointView<Lower>();
      // This only creates the lower triangle.
      for( j = 0; j < P ; j++ ){
        if(j < incCont){
          Sigma( j + 1, j + 1 ) += XtX( j + 1, j + 1 ) / N / g( gMap(j) ) ;
          for( k = 0; k < j; k++ ){
            Sigma( j + 1, k + 1 ) += XtX( j + 1, k + 1 ) / N / g( gMap(j) );
          }
        }else{
          Sigma( j + 1, j + 1 ) += 1/g( gMap(j) );
        }
      }
      
      Sigma = sig2 * Sigma.selfadjointView<Lower>().llt().solve( MatrixXd::Identity( P + 1, P + 1 ) );
      mu = Sigma * Xty / sig2;
      beta = random_multivariate_normal( mu, Sigma );

      yResid = ym - Xm * beta;
      SSqG = SSqG * 0.0;
      if(incCont){
        SSqG( gMap(0) ) += ( beta.block(1,0,incCont,1).transpose() * XtX.block(1, 1, incCont, incCont) * beta.block(1, 0, incCont, 1) )(0,0) / N;
       }
      for( j = incCont ; j < P ; j++ ){
        SSqG( gMap(j) ) += beta( j + 1, 0 ) * beta( j + 1, 0 );
      }
    }
    
    // sample sig2
    SSq = ( yResid.transpose() * yResid )(0,0);
    if(!nullModel) SSq += sum(SSqG / g);
    sig2 = 1 / Rf_rgamma( 0.5 * ( N + P*(!nullModel) ), 2 / SSq );

    // sample g
    for( j = 0 ; j < nGs ; j++ ){
      if(nullModel){
        g(j) = NA_REAL; 
      }else{
        g(j) = 1 / Rf_rgamma( 0.5 * ( gMapCounts(j)*(!nullModel) + 1 ) , 2 / ( SSqG(j)/sig2 + rscale(j)*rscale(j) ) );
      }
    }
    
    // copy to chain
    if(!(i % thin)){
      colCounter = 0;
      chains( rowCounter, 0 ) = beta( 0, 0 );
      for( j = 0 ; j < P ; j++ ){ // beta parameters
        if( !ignoreCols( j ) ){ // ignore filtered parameters
          chains( rowCounter, ++colCounter ) = beta( j + 1 , 0 );
        }
      }
      chains( rowCounter, ++colCounter ) = sig2;
      for( j = 0 ; j < nGs ; j++ ){ // copy g parameters
        chains( rowCounter, ++colCounter ) = g(j);
      }
      rowCounter++;
    }
  
  } // end Gibbs sampler

  return chains;

}
