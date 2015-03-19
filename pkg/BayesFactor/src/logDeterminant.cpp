#include <RcppEigen.h>

using namespace Rcpp;
using Eigen::MatrixXd;
using Eigen::VectorXd;

// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::export]]
double log_determinant_pos_def(Eigen::MatrixXd A)
{
  const VectorXd Dvec(A.ldlt().vectorD());
  return Dvec.array().log().sum();
}
