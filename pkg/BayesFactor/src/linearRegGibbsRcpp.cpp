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
NumericMatrix GibbsLinearRegRcpp(const int iterations, const NumericVector y, const NumericMatrix X, const double r, const double sig2start, const bool nullModel, const int progress, const Function callback, const double callbackInterval)
{
    RNGScope scope;
    // setting last_cb to the beginning of the epoch 
    // ensures that the callback is called once, first
    time_t last_cb = static_cast<time_t>(int(0));

    const Map<MatrixXd> Xm(as<Map<MatrixXd> >(X));
    
    const int P = X.ncol();
    const int N = y.size();
    int i = 0, j = 0;
    NumericMatrix XcolMeans(1, P);
    MatrixXd Sigma(MatrixXd(P, P).setZero());
    MatrixXd beta(MatrixXd(P, 1).setZero());
    MatrixXd XtCny(MatrixXd(P, 1).setZero());
    MatrixXd mu(MatrixXd(P, 1).setZero());
    MatrixXd Cny(MatrixXd(N, 1).setZero());
    MatrixXd yResid(MatrixXd(N, 1).setZero());
    MatrixXd CnX(MatrixXd(N, P).setZero());
    const MatrixXd XtXoN(MatrixXd(P, P).setZero().
      selfadjointView<Lower>().rankUpdate(Xm.adjoint(), 1 / (1.0*N) ));

    double g=1, sig2 = sig2start, SSq=0, meanY = mean(y), SSq_temp = 0;
    const double sumy2 = sum(y * y);
    
    for( i = 0 ; i < X.ncol() ; i++ ){
      XcolMeans(0, i) = sum( X( _, i) ) / ( N * 1.0 );
    }

    // compute centered values
    for( i = 0 ; i < N ; i++ ){
      Cny(i) = y(i) - meanY; 
      for( j = 0 ; j < P ; j++){
        CnX(i, j) = X(i, j) - XcolMeans(0, j);
      }
    }
    
    XtCny = Xm.transpose() * Cny;
    
    const MatrixXd XtCnX(MatrixXd(P, P).setZero().
      selfadjointView<Lower>().rankUpdate(CnX.adjoint()));

    // create progress bar
    class Progress p(iterations, (bool) progress);

    // Create matrix for chains
    NumericMatrix chains(iterations, P + 2);

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
   
      // Sample beta
      if(!nullModel){
        Sigma = ( XtCnX + ( XtXoN / g ) ) / sig2;
        Sigma = Sigma.llt().solve(MatrixXd::Identity(P,P));
        mu = ( Sigma / sig2 ) * XtCny;
        beta = random_multivariate_normal( mu, Sigma );
        
        SSq_temp = ( beta.transpose() * XtXoN * beta )(0,0);
        SSq = SSq_temp / g;
        yResid = Cny - Xm * beta;
        SSq += ( yResid.transpose() * yResid )(0,0);
      }else{
        SSq_temp = 0;
        SSq = sumy2;
      }

      // Sample sig2
      sig2 = 1 / Rf_rgamma( ( N + P*(!nullModel) ) / 2, 2 / SSq );

      // Sample g
		  SSq =  SSq_temp / sig2 + r*r;
      g  = 1/Rf_rgamma( 0.5 * ( P*(!nullModel) + 1 ) , 2 / SSq );

      // Copy samples to chain
      for( j = 0 ; j < P ; j++ )
        chains(i, j) = beta(j, 0);
		  chains(i, P) = sig2;
		  chains(i, P + 1) = g;
    }
    
    return chains;
}


