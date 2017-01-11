#include "bfcommon.h"

using namespace Rcpp;

// return log(1 + exp(x)), preventing cancellation and overflow */
// From http://www.codeproject.com/KB/recipes/avoiding_overflow.aspx
// [[Rcpp::export]]
double log1pExp(double x)
{
    const double LOG_DBL_EPSILON = log(DBL_EPSILON);
    const double LOG_ONE_QUARTER = log(0.25);

    if (x > -LOG_DBL_EPSILON)
    {
        // log(exp(x) + 1) == x to machine precision
        return x;
    }
    else if (x > LOG_ONE_QUARTER)
    {
        return log( 1.0 + exp(x) );
    }
    else
    {
        // Prevent loss of precision that would result from adding small argument to 1.
        return log1p( exp(x) );
    }
}

// Compute log(exp(x) + exp(y))
// [[Rcpp::export]]
double logExpXplusExpY( const double x, const double y )
{
  return x + log1pExp( y - x );
}

// Compute log(exp(x) - exp(y))
// [[Rcpp::export]]
double logExpXminusExpY( const double x, const double y )
{
  return x + Rf_pexp( x - y, 1, true, true );
}


