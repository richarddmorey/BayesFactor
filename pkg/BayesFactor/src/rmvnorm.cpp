#include <RcppEigen.h>
#include "bfcommon.h"

using namespace Rcpp;

// [[Rcpp::depends(RcppEigen)]]
Eigen::MatrixXd random_multivariate_normal(const Eigen::MatrixXd mu, const Eigen::MatrixXd Sigma) 
{

  int P = mu.rows(), i = 0;
  Eigen::MatrixXd y(Eigen::MatrixXd(P, 1).setZero());
  Eigen::MatrixXd z(Eigen::MatrixXd(P, 1).setZero());
  
  for( i = 0 ; i < P ; i++ )
    z(i, 0) = Rf_rnorm( 0, 1 );
  
  y = mu + Sigma.llt().matrixL() * z;

  return y;
}

