#include "bfcommon.h"
#include <RcppEigen.h>

using namespace Rcpp;
using Eigen::MatrixXd;
using Eigen::Map;

// [[Rcpp::depends(RcppEigen)]]
double jzs_mc_marg_like(const NumericVector g, const double sumSq, const NumericVector Cny, const NumericMatrix CnX, const NumericMatrix XtCnX, const NumericMatrix CnytCnX, const NumericVector rscale, const IntegerVector gMap, const NumericVector gMapCounts, const NumericMatrix priorX, const double logDetPriorX, const int incCont)
{
  double ans = 0, sumInvGammaDens = 0;
  const NumericVector q = log(g);
  const NumericVector limits(2);
  const int nGs = g.size();
  int i = 0;


  // The jzs_log_marginal_posterior_logg function adds
  // this in, but we don't need it so substract it back out
  for( i = 0 ; i < nGs ; i++ ){
    sumInvGammaDens += dinvgamma1_Rcpp(g(i), 0.5, rscale(i) * rscale(i) / 2);
  }
  
  // Warning: this function already has the null likelihood built in.
  ans = ( jzs_log_marginal_posterior_logg(q, sumSq, Cny, CnX, XtCnX, CnytCnX, rscale, gMap, gMapCounts, priorX, incCont, false, limits, 0) )["d0g"]; 
  
  // substract sum(q) for the log transformation
  return ans - sum(q) - sumInvGammaDens + .5*logDetPriorX;
}

double jzs_mc_marg_like2(const NumericVector g, const double N, const double sumSq, const NumericMatrix XtCnX0, const NumericMatrix XtCny0, const IntegerVector gMap, const NumericMatrix priorX, const double logDetPriorX, const int incCont)
{
  int P = XtCnX0.ncol();
  const MatrixXd XtCnX(as<Map<MatrixXd> >(XtCnX0));  
  MatrixXd W( XtCnX );
  MatrixXd WInvXtCny( MatrixXd(P, 1).setZero() );
  const MatrixXd XtCny( as<Map<MatrixXd> >(XtCny0) );

  double ldetS=0,ldetW,top,bottom1,bottom2,q=0;
	int i,j;
  	
  ldetS += -logDetPriorX;
  for(i=0;i<incCont;i++){
    ldetS += log( g( gMap(i) ) );
    for(j=0;j<=i;j++){
      W( i, j ) += priorX( i, j ) / g( gMap( i ) ); 
    }
  }
  
  for(i=incCont;i<P;i++)
	{
		ldetS += log(g(gMap(i)));
		W(i,i) += 1/g(gMap(i));
	}

  ldetW = -log_determinant_pos_def(W);
  WInvXtCny = W.llt().solve(XtCny);
	for( i = 0 ; i < P ; i++ )
    q += XtCny( i, 0 ) * WInvXtCny( i, 0 );
	
	top = 0.5*ldetW;
	bottom1 = 0.5*(N-1) * log1p( -q / sumSq );
	bottom2 = 0.5*ldetS;
	
	return(top-bottom1-bottom2);	
}



double jzs_importance_marg_like(const NumericVector q, const NumericVector mu, const NumericVector sig, const double sumSq, const NumericVector Cny, const NumericMatrix CnX, const NumericMatrix XtCnX, const NumericMatrix CnytCnX, const NumericVector rscale, const IntegerVector gMap, const NumericVector gMapCounts, const NumericMatrix priorX, const double logDetPriorX, const int incCont)
{
  double ans = 0, sumNormDens = 0;
  const NumericVector limits(2);
  const int nGs = q.size();
  int i = 0;

  
  // The jzs_log_marginal_posterior_logg function adds
  // this in, but we don't need it so substract it back out
  for( i = 0 ; i < nGs ; i++ ){
    sumNormDens += Rf_dnorm4(q(i), mu(i), sig(i), 1);
  }
  
  // Warning: this function already has the null likelihood built in.
  ans = ( jzs_log_marginal_posterior_logg(q, sumSq, Cny, CnX, XtCnX, CnytCnX, rscale, gMap, gMapCounts, priorX, incCont, false, limits, 0) )["d0g"]; 
  
  // substract sum(q) for the log transformation
  return ans - sumNormDens + .5*logDetPriorX;
}




