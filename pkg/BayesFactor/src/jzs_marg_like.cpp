#include "bfcommon.h"

using namespace Rcpp;

// [[Rcpp::depends(Rcpp)]]
// [[Rcpp::export]]
double jzs_mc_marg_like(const NumericVector g, const double sumSq, const NumericVector Cny, const NumericMatrix CnX, const NumericMatrix CnytCnX, const NumericVector rscale, const IntegerVector gMap, const NumericVector gMapCounts, const NumericMatrix priorX, const double logDetPriorX, const int incCont)
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
  ans = ( jzs_log_marginal_posterior_logg(q, sumSq, Cny, CnX, CnytCnX, rscale, gMap, gMapCounts, priorX, incCont, false, limits, 0) )["d0g"]; 
  
  // substract sum(q) for the log transformation
  return ans - sum(q) - sumInvGammaDens + .5*logDetPriorX;
}

double jzs_importance_marg_like(const NumericVector q, const NumericVector mu, const NumericVector sig, const double sumSq, const NumericVector Cny, const NumericMatrix CnX, const NumericMatrix CnytCnX, const NumericVector rscale, const IntegerVector gMap, const NumericVector gMapCounts, const NumericMatrix priorX, const double logDetPriorX, const int incCont)
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
  ans = ( jzs_log_marginal_posterior_logg(q, sumSq, Cny, CnX, CnytCnX, rscale, gMap, gMapCounts, priorX, incCont, false, limits, 0) )["d0g"]; 
  
  // substract sum(q) for the log transformation
  return ans - sumNormDens + .5*logDetPriorX;
}




