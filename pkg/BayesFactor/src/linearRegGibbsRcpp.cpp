#include "progress.h"
#include "bfcommon.h"
#include <RcppEigen.h>
#include <time.h>

using namespace Rcpp;

// [[Rcpp::depends(RcppEigen)]]

// void GibbsLinearReg(double *chains, int iters, double *Cny, double *X, double *XtX, double *XtCnX, double *XtCny, int N, int P, double r, double sig2start, int progress, SEXP pBar, SEXP callback, SEXP rho)
// [[Rcpp::export]]
NumericMatrix GibbsLinearRegRcpp(int iterations, NumericVector y, NumericMatrix X, double r, double sig2start, bool nullModel, int progress, Function callback, double callbackInterval)
{
    RNGScope scope;  
    // setting last_cb to the beginning of the epoch 
    // ensures that the callback is called once, first
    time_t last_cb = static_cast<time_t>(int(0));

    const Eigen::Map<Eigen::MatrixXd> Xm(as<Eigen::Map<Eigen::MatrixXd> >(X));
    
    int P = X.ncol(), N = y.size(), i = 0, j = 0;
    NumericMatrix XcolMeans(1, P);
    Eigen::MatrixXd Sigma(Eigen::MatrixXd(P, P).setZero());
    Eigen::MatrixXd beta(Eigen::MatrixXd(P, 1).setZero());
    Eigen::MatrixXd XtCny(Eigen::MatrixXd(P, 1).setZero());
    Eigen::MatrixXd mu(Eigen::MatrixXd(P, 1).setZero());
    Eigen::MatrixXd Cny(Eigen::MatrixXd(N, 1).setZero());
    Eigen::MatrixXd yResid(Eigen::MatrixXd(N, 1).setZero());
    Eigen::MatrixXd CnX(Eigen::MatrixXd(N, P).setZero());
    Eigen::MatrixXd XtXoN(Eigen::MatrixXd(P, P).setZero().
      selfadjointView<Eigen::Lower>().rankUpdate(Xm.adjoint(), 1 / (1.0*N) ));

    double g=1, sig2 = sig2start, SSq=0, meanY = mean(y), SSq_temp = 0;

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
    
    Eigen::MatrixXd XtCnX(Eigen::MatrixXd(P, P).setZero().
      selfadjointView<Eigen::Lower>().rankUpdate(CnX.adjoint()));

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
        Sigma = Sigma.llt().solve(Eigen::MatrixXd::Identity(P,P));
        mu = ( Sigma / sig2 ) * XtCny;
        beta = random_multivariate_normal( mu, Sigma );
      }

      SSq_temp = ( beta.transpose() * XtXoN * beta )(0,0);
      
      // Sample sig2
      SSq =  SSq_temp / g;
      yResid = Xm * beta - Cny;
      for( j = 0; j < N; j++ )
        SSq += yResid( j, 0 ) * yResid( j, 0 );
      sig2 = 1 / Rf_rgamma( N / 2, 2 / SSq );
	      
      // Sample g
		  SSq =  SSq_temp / sig2 + r*r;
      g  = 1/Rf_rgamma( 1, 2 / SSq );

      // Copy samples to chain
      for( j = 0 ; j < P ; j++ )
        chains(i, j) = beta(j, 0);	
		  chains(i, P) = sig2;	
		  chains(i, P + 1) = g;	
    }
    
    return chains;
}


