// Analysis of Visual Array Data
// Richard D. Morey (richarddmorey@gmail.com)

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
#include <R_ext/Lapack.h>
#include <R_ext/Random.h>
#include <R_ext/BLAS.h>

void R_CheckUserInterrupt(void);



SEXP doModel(SEXP designMatrix, SEXP Xdim, SEXP dataVec, SEXP NsParGrps, SEXP Nsims, SEXP burnin, SEXP progress, SEXP pBar, SEXP rho, SEXP MCMCCI, SEXP rscale);
void WM2_rmvGaussianC(double *mu, double *Sigma, int p);
double WM2_quadform(double *x, double *A, int N, int incx);


void WM2_rmvGaussianC(double *mu, double *Sigma, int p)
{
  double ans[p];
  int info, psqr,j=0, intOne=1;
  double *scCp, one = 1; //zero = 0;
  
  psqr = p * p;
  scCp = Memcpy(Calloc(psqr,double), Sigma, psqr);

  F77_NAME(dpotrf)("U", &p, scCp, &p, &info);
  if (info)
    error("Sigma matrix is not positive-definite");
  
  GetRNGstate();
  for(j=0;j<p;j++)
    {
      ans[j] = rnorm(0,1);
    }
  F77_NAME(dtrmv)("U","N","N", &p, scCp, &p, ans, &intOne);
  F77_NAME(daxpy)(&p, &one, ans, &intOne, mu, &intOne);
  PutRNGstate();
  Free(scCp);
}



SEXP doModel(SEXP designMatrix, SEXP Xdim, SEXP dataVec, SEXP NsParGrps, SEXP Nsims, SEXP burnin, SEXP progress, SEXP pBar, SEXP rho, SEXP MCMCCI, SEXP rscale)
{
  int N = INTEGER_POINTER(Xdim)[1];
  int p = INTEGER_POINTER(Xdim)[2]-1;
  double *Nks = INTEGER_POINTER(NsParGrps);
  double *y   = REAL(dataVec);
  int M = INTEGER_VALUE(Nsims);
  int burn = INTEGER_VALUE(burnin);
  double dRscale = REAL(rscale)[0];
  int m=0, k=0, K=length(NsParGrps);
  double *X = REAL(designMatrix);

  SEXP R_fcall, sampCounter, betaChain, gChain, sig2Chain, dens;
  int *pSampCounter = INTEGER_POINTER(sampCounter);
  int iProgress = INTEGER_VALUE(progress);
  double *B = Calloc((p+1)*(p+1),double);
  
  double *XtX = Calloc((p+1)*(p+1),double);

  PROTECT(R_fcall = lang2(pBar, R_NilValue));
  PROTECT(sampCounter = NEW_INTEGER(1));

  PROTECT(betaChain = allocMatrix(REALSXP, p+1, M));
  PROTECT(gChain = allocMatrix(REALSXP, K, M));
  PROTECT(dens = allocMatrix(REALSXP, K, M-1));
  PROTECT(sig2Chain = allocVector(VECSXP, M));
  
  
   for(m=1;m<M;m++)
     {
  
       //printf("Iteration %d\n",i);
       //Check to see if we need to update the progress bar
       if(iProgress && !((i+1)%iProgress)){
	 pSampCounter[0]=i+1;
	 SETCADR(R_fcall, sampCounter);
	 eval(R_fcall, rho); //Update the progress bar
       }

       for(k=0;k<K;k++)
	 {
	   
	   rgamma((double)(Nks)[k]/2+.5,);
	 }

       

     }

   Free(B);
   Free(XtX);

   UNPROTECT(6);

}

