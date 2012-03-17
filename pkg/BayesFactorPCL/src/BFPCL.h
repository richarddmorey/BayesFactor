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

void gibbsEqVariance(double *y, int *N, int J, int I, double lambda, int iterations, double *chains, double sdMetropSig2, double sdMetropTau, int progress, SEXP pBar, SEXP rho);

void gibbsOneWayAnova(double *y, int *N, int J, int sumN, int *whichJ, double rscale, int iterations, double *chains, double *CMDE, SEXP debug, int progress, SEXP pBar, SEXP rho);

double sampleSig2EqVar(double sig2, double *mu, double tau, double *yBar, double *SS, int *N, int sumN, int J, double sdMetrop, double *acc);

double logFullCondTauEqVar(double tau, double *mu, double sig2, double *yBar, double *SS, int *N, int J, double lambda);

double logFullCondSig2EqVar(double sig2, double *mu, double tau, double *yBar, double *SS, int *N, int sumN, int J);

double sampleTauEqVar(double tau, double *mu, double sig2, double *yBar, double *SS, int *N, int J, double lambda, double sdMetrop, double *acc);

void gibbsEqVarianceM2(double *y, int *N, int J, int I, double alpha, double beta, int iterations, double *chains, double *debug, double sdMetropg, double sdDecorr, int newtonSteps, int progress, SEXP pBar, SEXP rho);

void decorrStepEqVarM2(double sumg, double sig2g, int J, double *g, double *sig2, double sdDecorr, double *acc);

double newtonMethodEqVar2(double xt, double c1, double c2, double c3, int iterations, int max);

double samplegEqVarM2(double g, double sig2, double mu, double sig2g, double yBar, double sumy2, int N, double sdMetrop, double *acc);

double fullCondgEqVarM2(double g, double sig2, double mu, double sig2g, double yBar, double sumy2, int N);

double thetaLogLikeAR(double theta, double mu, double delta, double sig2, double g, double *y, int N, double *t, double alphaTheta, double betaTheta);

double sampThetaAR(double theta, double mu, double delta, double sig2, double g, double *y, int N, double *t, double alphaTheta, double betaTheta , double sdmet);

void gibbsTwoSampleAR(double *y, int N, double *t, double rscale, double alphaTheta, double betaTheta, double loArea, double upArea, int iterations, double sdmet, double *chains, int progress, SEXP pBar, SEXP rho);

SEXP RgibbsTwoSampleAR(SEXP yR, SEXP NR, SEXP tR, SEXP rscaleR, SEXP alphaThetaR, SEXP betaThetaR, SEXP loArea, SEXP upArea, SEXP iterationsR, SEXP sdmet, SEXP progressR, SEXP pBar, SEXP rho);

double MCTwoSampleAR(double *y, int N, double *t, double rscale, double alphaTheta, double betaTheta, int iterations);

SEXP RMCTwoSampleAR(SEXP yR, SEXP NR, SEXP tR, SEXP rscaleR, SEXP alphaThetaR, SEXP betaThetaR, SEXP iterationsR);


SEXP RjeffSamplerNwayAov(SEXP Riters, SEXP RXtCX, SEXP RXtCy, SEXP RytCy, SEXP RN, SEXP RP, SEXP RnGs, SEXP RgMap, SEXP Ra, SEXP Rb, SEXP progressR, SEXP pBar, SEXP rho);

double jeffSamplerNwayAov(double *samples, int iters, double *XtCX, double *XtCy, double ytCy, int N, int P,int nGs, int *gMap, double *a, double *b, int progress, SEXP pBar, SEXP rho);

double jeffmlikeNWayAov(double *XtCX, double *XtCy, double ytCy, int N, int P, double *g);

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

SEXP RgibbsTwoSampleARmulti(SEXP yR, SEXP NR, SEXP NobsR, SEXP tR, SEXP rscaleR, SEXP betaThetaR, SEXP iterationsR, SEXP sdmetR, SEXP progressR, SEXP pBar, SEXP rho);

void gibbsTwoSampleARmulti(double *y, int N, int *Nobs, double *t, double rscale, double betaTheta, int iterations, double sdmet, double *chains, int progress, SEXP pBar, SEXP rho);

double sampThetaARmulti(double theta, double *mu0, double *beta, double *sig2, double g, double *y, int N, double *t, int *Nobs, int *cumobs, int maxobs, double betaTheta , double sdmet);

double thetaLogLikeARmulti(double theta, double *mu0, double *beta, double *sig2, double g, double *y, int N, double *t, int *Nobs, int *cumobs, int maxobs, double betaTheta);

SEXP RthetaLogLikeARmulti(SEXP thetaR, SEXP mu0R, SEXP betaR, SEXP sig2R, SEXP gR, SEXP yR, SEXP NR, SEXP tR, SEXP NobsR, SEXP cumobsR, SEXP maxobsR, SEXP betaThetaR);





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

