#include "bfcommon.h"
#include <RcppEigen.h>

using namespace Rcpp;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using Eigen::Map;
using Eigen::Lower;

// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::export]]
Rcpp::List jzs_log_marginal_posterior_logg(const NumericVector q, const double sumSq, const int N, const NumericMatrix XtCnX0, const NumericMatrix CnytCnX0, const NumericVector rscale, const IntegerVector gMap, const NumericVector gMapCounts, const NumericMatrix priorX, const int incCont, const bool limit, const NumericVector limits, const int which)
{
 
  const int P = XtCnX0.ncol();
  const int nGs = q.size();
  double d0g = NA_REAL;
  int i = 0, j = 0;
  
  const NumericVector g = exp(q);
  
  if( nGs != rscale.size() ){
    Rcpp::stop("length mismatch: q and r");
  }
  if( nGs != gMapCounts.size() ){
    Rcpp::stop("length mismatch: q and gMapCount");
  }
  
  NumericVector d1g(nGs, NA_REAL);
  NumericVector d2g(nGs, NA_REAL);
  NumericVector ddensg(nGs, 0.0);
  NumericVector tempVec1(nGs, 0.0);
  NumericVector tempVec2(nGs, 0.0);
  NumericVector tempVec3(nGs, 0.0);
  NumericVector tempVec4(nGs, 0.0); 
  MatrixXd VgInv(MatrixXd(P, P).setZero());
  const MatrixXd XtCnX(as<Map<MatrixXd> >(XtCnX0));
  const MatrixXd CnytCnX(as<Map<MatrixXd> >(CnytCnX0));
  MatrixXd gInv( XtCnX );
  MatrixXd CnytCnXVg(MatrixXd(1, P).setZero());
  MatrixXd CnytCnXVg2(MatrixXd(1, P).setZero());
  MatrixXd VgInv2(MatrixXd(P, P));

  double logDetVg, yXVXy, sumLogg = 0, sumInvGammaDens = 0;
  
  if( ( P == 1 ) && incCont ){
    Rcpp::stop("Inappropriate use of Gaussian approximation with single column, continuous X.");
  }
  
  if( limit ){
    for( i = 0; i < nGs ; i++ ){
      if( ( q(i) < limits(0) ) | ( q(i) > limits(1) ) ){
        return Rcpp::List::create(Rcpp::Named("d0g") = -INFINITY,
                          Rcpp::Named("d1g") = wrap(d1g),
                          Rcpp::Named("d2g") = wrap(d2g));;
      }
    }
  }
  
  // Build g matrix
  for( i = incCont ; i < P ; i++ ){
    sumLogg += q( gMap(i) );
    gInv(i,i) +=  1 / g( gMap(i) );
  }
  if(incCont){ // Continuous covariates included
    if( priorX.nrow() != incCont )
      Rcpp::stop("priorX matrix size does not match argument incCont.");
    if( priorX.nrow() != priorX.ncol() )
      Rcpp::stop("priorX matrix must be square.");
    for( i = 0; i < incCont ; i++ ){
      sumLogg += q( gMap(i) );
      for( j = 0 ; j <= i ; j++ ){
        gInv(i,j) += priorX(i,j) / g(gMap(i));
      }
    }
  }

  VgInv = gInv.selfadjointView<Lower>().llt().solve(MatrixXd::Identity(P, P));
  CnytCnXVg =  CnytCnX * VgInv;
  yXVXy = ( CnytCnXVg * CnytCnX.transpose() )(0,0);
  
  if(which == 0 || which == -1){ // (log) marginal posterior
    logDetVg = -log_determinant_pos_def( VgInv );

    for( i = 0 ; i < nGs ; i++ ){
      sumInvGammaDens += dinvgamma1_Rcpp(g(i), 0.5, rscale(i) * rscale(i) / 2);
    }
  
    d0g = -0.5 * sumLogg - 0.5 * logDetVg - 0.5*(N-1)*log1p( -yXVXy/sumSq ) + sumInvGammaDens + sum(q);
  }

  if(which > 0 || which == -1){ // first derivative of log
    
    for( i = 0 ; i < nGs ; i++ ){
      ddensg(i) = g(i) * ddinvgamma1_Rcpp(g(i), 0.5, rscale(i) * rscale(i) / 2);
    }
    
    for( i = 0 ; i < P ; i++ ){
      tempVec1( gMap(i) ) += VgInv.diagonal()(i) / g( gMap(i) );
      tempVec2( gMap(i) ) += CnytCnXVg(0, i) * CnytCnXVg(0, i) / g( gMap(i) ) / (sumSq - yXVXy);
    }
    
    d1g = -0.5 * gMapCounts + 0.5 * tempVec1 + 0.5*(N-1) * tempVec2  + ddensg + 1.0;
  }
  
  if(which > 1 || which == -1){ // second derivative
    
    CnytCnXVg2 = CnytCnXVg * VgInv;
    VgInv2 = MatrixXd(P, P).setZero().selfadjointView<Lower>().rankUpdate(VgInv);
    
    for( i = 0 ; i < nGs ; i++ ){
      ddensg(i) = g(i) * g(i) * d2dinvgamma1_Rcpp(g(i), 0.5, rscale(i) * rscale(i) / 2);
    }

    for( i = 0 ; i < P ; i++ ){
      tempVec3( gMap(i) ) += VgInv2.diagonal()(i) / ( g( gMap(i) ) * g( gMap(i) ) );
      tempVec4( gMap(i) ) += CnytCnXVg(0, i) * CnytCnXVg2(0, i) / ( g( gMap(i) ) * g( gMap(i) ) * (sumSq - yXVXy) );
    }

    d2g = d1g + 0.5 * gMapCounts - 1.0 + ddensg + 0.5 * (tempVec3 - 2*tempVec1) + 0.5*(N-1) * ( tempVec2 * tempVec2 - 2*tempVec2  + 2*tempVec4 );

  }
  
  return Rcpp::List::create(Rcpp::Named("d0g") = d0g,
                          Rcpp::Named("d1g") = d1g,
                          Rcpp::Named("d2g") = d2g);
}
