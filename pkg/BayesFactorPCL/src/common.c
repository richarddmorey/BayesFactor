#include "BFPCL.h"

SEXP rmvGaussianR(SEXP mu, SEXP Sigma)
{
    SEXP ans;
    int *dims = INTEGER(getAttrib(Sigma, R_DimSymbol)), psqr, p;
    double *scCp, *ansp;//, one = 1, zero = 0;

    if (!isMatrix(Sigma) || !isReal(Sigma) || dims[0] != dims[1])
	error("Sigma must be a square, real matrix");
    
    p = dims[0];
    psqr=p*p;

    PROTECT(ans = allocVector(REALSXP, p));
    ansp = REAL(ans);
    Memcpy(ansp, REAL(mu), p);

    scCp = Memcpy(Calloc(psqr,double), REAL(Sigma), psqr);

    
    rmvGaussianC(ansp, scCp, p);
    
    Free(scCp);

    UNPROTECT(1);
    return ans;
}

void rmvGaussianC(double *mu, double *Sigma, int p)
{
  double ans[p];
  int info, psqr,j=0, intOne=1;
  double *scCp, one = 1; //zero = 0;
  
  psqr = p * p;
  scCp = Memcpy(Calloc(psqr,double), Sigma, psqr);

  F77_NAME(dpotrf)("L", &p, scCp, &p, &info);
  if (info){
	error("Nonzero info from dpotrf: Sigma matrix is not positive-definite");
  }
  //GetRNGstate();
  for(j=0;j<p;j++)
    {
      ans[j] = rnorm(0,1);
    }
  F77_NAME(dtrmv)("L","N","N", &p, scCp, &p, ans, &intOne);
  F77_NAME(daxpy)(&p, &one, ans, &intOne, mu, &intOne);
  //PutRNGstate();
  Free(scCp);
}



/**
 * Allocate a 3-dimensional array
 *
 * @param mode The R mode (e.g. INTSXP)
 * @param nrow number of rows
 * @param ncol number of columns
 * @param nface number of faces
 *
 * @return A 3-dimensional array of the indicated dimensions and mode
 */
SEXP alloc3Darray(SEXPTYPE mode, int nrow, int ncol, int nface)
{
    SEXP s, t;
    int n;

    if (nrow < 0 || ncol < 0 || nface < 0)
	error("negative extents to 3D array");
    if ((double)nrow * (double)ncol * (double)nface > INT_MAX)
	error("alloc3Darray: too many elements specified");
    n = nrow * ncol * nface;
    PROTECT(s = allocVector(mode, n));
    PROTECT(t = allocVector(INTSXP, 3));
    INTEGER(t)[0] = nrow;
    INTEGER(t)[1] = ncol;
    INTEGER(t)[2] = nface;
    setAttrib(s, R_DimSymbol, t);
    UNPROTECT(2);
    return s;
}

/*
double LogOnePlusX(double x)
{
    if (x <= -1.0)
    {
		error("Error: attempt to compute log(1+X) where X is %f\n.",x);
	}

    if (fabs(x) > 1/(10000.0))
    {
        // x is large enough that the obvious evaluation is OK
        return( log(1.0 + x) );
    }

    // Use Taylor approx. log(1 + x) = x - x^2/2 with error roughly x^3/3
    // Since |x| < 10^-4, |x|^3 < 10^-12, relative error less than 10^-8

    return( (-0.5*x + 1.0)*x );
}
*/

// Calculate log(1 + x), preventing loss of precision for small values of x.
// The input x must be larger than -1 so that log(1 + x) is real.
// Taken from http://www.codeproject.com/KB/recipes/avoiding_overflow.aspx
double LogOnePlusX(double x)
{
    if (x <= -1.0)
    {
        
        error("Attempt to compute log(1+x) on x<=-1!");
    }

	if (fabs(x) > 0.375)
    {
        // x is sufficiently large that the obvious evaluation is OK
        return log(1.0 + x);
    }

	// For smaller arguments we use a rational approximation
	// to the function log(1+x) to avoid the loss of precision
	// that would occur if we simply added 1 to x then took the log.

    const double p1 =  -0.129418923021993e+01;
    const double p2 =   0.405303492862024e+00;
    const double p3 =  -0.178874546012214e-01;
    const double q1 =  -0.162752256355323e+01;
    const double q2 =   0.747811014037616e+00;
    const double q3 =  -0.845104217945565e-01;
    double t, t2, w;

    t = x/(x + 2.0);
    t2 = t*t;
    w = (((p3*t2 + p2)*t2 + p1)*t2 + 1.0)/(((q3*t2 + q2)*t2 + q1)*t2 + 1.0);
    return 2.0*t*w;
}

// return log(1 + exp(x)), preventing cancellation and overflow */
// From http://www.codeproject.com/KB/recipes/avoiding_overflow.aspx
double LogOnePlusExpX(double x)
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
        return LogOnePlusX( exp(x) );
    }
}
