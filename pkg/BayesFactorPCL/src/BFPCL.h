#include <R.h>
#include <Rmath.h>  
#include <Rdefines.h>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <R_ext/Utils.h>
#include <Rversion.h>
#include <Rconfig.h>
#include <R_ext/Constants.h>
#include <R_ext/Random.h>
#include <R_ext/Lapack.h>
#include <R_ext/BLAS.h>


void gibbsOneSample(double *y, int N, double rscale, int iterations, double *chains, int doInterval, double *interval, int progress, SEXP pBar, SEXP rho);

void gibbsOneWayAnova(double *y, int *N, int J, int sumN, int *whichJ, double rscale, int iterations, double *chains, double *CMDE, SEXP debug, int progress, SEXP pBar, SEXP rho);

double sampleSig2EqVar(double sig2, double *mu, double tau, double *yBar, double *SS, int *N, int sumN, int J, double sdMetrop, double *acc);

SEXP RjeffSamplerNwayAov(SEXP Riters, SEXP RXtCX, SEXP RXtCy, SEXP RytCy, SEXP RN, SEXP RP, SEXP RnGs, SEXP RgMap, SEXP Ra, SEXP Rb, SEXP progressR, SEXP pBar, SEXP rho);
double jeffSamplerNwayAov(double *samples, double *gsamples, int iters, double *XtCX, double *XtCy, double ytCy, int N, int P,int nGs, int *gMap, double *a, double *b, int progress, SEXP pBar, SEXP rho);

SEXP RimportanceSamplerNwayAov(SEXP Riters, SEXP RXtCX, SEXP RXtCy, SEXP RytCy, SEXP RN, SEXP RP, SEXP RnGs, SEXP RgMap, SEXP Ra, SEXP Rb, SEXP Rmu, SEXP Rsig, SEXP progressR, SEXP pBar, SEXP rho);
double importanceSamplerNwayAov(double *samples, double *gsamples, int iters, double *XtCX, double *XtCy, double ytCy, int N, int P,int nGs, int *gMap, double *a, double *b, double *mu, double *sig, int progress, SEXP pBar, SEXP rho);


double jeffmlikeNWayAov(double *XtCX, double *XtCy, double ytCy, int N, int P, double *g);
SEXP RjeffmlikeNWayAov(SEXP XtCXR, SEXP XtCyR, SEXP ytCyR, SEXP NR, SEXP PR, SEXP gR);

SEXP RLogMeanExpLogs(SEXP Rv, SEXP Rlen);

SEXP RLogCumMeanExpLogs(SEXP Rv, SEXP Rlen);
void logCumMeanExpLogs(double *v, int len, double *ret);


double logSumExpLogs(double *v, int len);

double logMeanExpLogs(double *v, int len);

double LogOnePlusExpX(double x);

double LogOnePlusX(double x);

double quadform(double *x, double *A, int N, int incx, int LDA);

double quadform2(double *x, double *A, int N, int incx, int LDA);

SEXP rmvGaussianR(SEXP mu, SEXP Sigma);

void rmvGaussianC(double *mu, double *Sigma, int p);

int InvMatrixUpper(double *A, int p);

double matrixDet(double *A, int N, int LDA, int doLog);

void debugPrintMatrix(double *X, int rows, int cols);

void debugPrintVector(double *x, int len);

void GibbsNwayAov(double *chains, int iters, double *y, double *X, double *XtX, double *Xty, int N, int P, int nGs, int *gMap, double *r, int progress, SEXP pBar, SEXP rho);

SEXP RGibbsNwayAov(SEXP Riters, SEXP Ry, SEXP RX, SEXP RXtX, SEXP RXty, SEXP RN, SEXP RP, SEXP RnGs, SEXP RgMap, SEXP Rr, SEXP progressR, SEXP pBar, SEXP rho);

SEXP RGibbsLinearReg(SEXP Riters, SEXP RCny, SEXP RX, SEXP RXtX, SEXP RXtCnX, SEXP RXtCny, SEXP RN, SEXP RP, SEXP Rr, SEXP sig2start, SEXP progressR, SEXP pBar, SEXP rho);

void GibbsLinearReg(double *chains, int iters, double *Cny, double *X, double *XtX, double *XtCnX, double *XtCny, int N, int P, double r, double sig2start, int progress, SEXP pBar, SEXP rho);


#define AZERO(x, n) {int _I_, _SZ_ = (n); for(_I_ = 0; _I_ < _SZ_; _I_++) (x)[_I_] = 0;}

/**
 * Symmetrize a matrix by copying the strict upper triangle into the
 * lower triangle.
 *
 * @param a pointer to a matrix in Fortran storage mode
 * @param nc number of columns (and rows and leading dimension) in the matrix
 *
 * @return a, symmetrized
 */
static R_INLINE double*
internal_symmetrize(double *a, int nc)
{
    int i,j;
    for (i = 1; i < nc; i++)
	for (j = 0; j < i; j++)
	    a[i + j*nc] = a[j + i*nc];
    return a;
}

