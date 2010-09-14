
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
#include <R_ext/Random.h>
#include <R_ext/BLAS.h>


void gibbsOneSample(double *y, int N, double rscale, int iterations, double *chains, int progress, SEXP pBar, SEXP rho);
void gibbsEqVariance(double *y, int *N, int J, int I, double lambda, int iterations, double *chains, double sdMetropSig2, double sdMetropTau, int progress, SEXP pBar, SEXP rho);
void gibbsOneWayAnova(double *y, int *N, int J, int sumN, int *whichJ, double rscale, int iterations, double *chains, int progress, SEXP pBar, SEXP rho);
double sampleSig2EqVar(double sig2, double *mu, double tau, double *yBar, double *SS, int *N, int sumN, int J, double sdMetrop, double *acc);
double logFullCondTauEqVar(double tau, double *mu, double sig2, double *yBar, double *SS, int *N, int J, double lambda);
double logFullCondSig2EqVar(double sig2, double *mu, double tau, double *yBar, double *SS, int *N, int sumN, int J);
double sampleTauEqVar(double tau, double *mu, double sig2, double *yBar, double *SS, int *N, int J, double lambda, double sdMetrop, double *acc);

double quadform(double *x, double *A, int N, int incx);
SEXP rmvGaussianR(SEXP mu, SEXP Sigma);
void rmvGaussianC(double *mu, double *Sigma, int p);
int InvMatrixUpper(double *A, int p);
double matrixDet(double *A, int N, int log);

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


SEXP RgibbsOneSample(SEXP yR, SEXP NR, SEXP rscaleR, SEXP iterationsR, SEXP progressR, SEXP pBar, SEXP rho)
{
	
	int npars = 5,iterations = INTEGER_VALUE(iterationsR);
	int N = INTEGER_VALUE(NR), progress = INTEGER_VALUE(progressR);
	double rscale = NUMERIC_VALUE(rscaleR);
	double *y = REAL(yR);
	
	
	SEXP chainsR;
	PROTECT(chainsR = allocMatrix(REALSXP, npars, iterations));

	gibbsOneSample(y, N, rscale, iterations, REAL(chainsR), progress, pBar, rho);
	
	UNPROTECT(1);
	
	return(chainsR);
}

void gibbsOneSample(double *y, int N, double rscale, int iterations, double *chains, int progress, SEXP pBar, SEXP rho)
{
	int i=0;
	double yBar=0,sumy2=0,rscaleSq=pow(rscale,2),densDelta=0,meanDelta=0,varDelta=1;
	double mu=0,sig2=1,g=1;
	double meanMu=0, varMu=0, shapeSig2=(1.0*N)/2+0.5,scaleSig2=1,scaleg=1;
	int npars = 5;
	
	// progress stuff
	SEXP sampCounter, R_fcall;
	int *pSampCounter;
    PROTECT(R_fcall = lang2(pBar, R_NilValue));
	PROTECT(sampCounter = NEW_INTEGER(1));
	pSampCounter = INTEGER_POINTER(sampCounter);
	
	GetRNGstate();
	
	for(i=0;i<N;i++)
	{
		yBar += y[i]/((float)(N));
		sumy2 += pow(y[i],2);
	}
	
	// start MCMC
	for(i=0; i<iterations; i++)
	{
		R_CheckUserInterrupt();
	
		
		//Check progress
		if(progress && !((i+1)%progress)){
			pSampCounter[0]=i+1;
			SETCADR(R_fcall, sampCounter);
			eval(R_fcall, rho); //Update the progress bar
		}
	
	
		// sample mu
		meanMu = yBar * N/(1.0*N+1/g);
		meanDelta = meanMu/sqrt(sig2);
		varDelta = 1/(1.0*N + 1/g);
		varMu  = sig2*varDelta;
		mu = rnorm(meanMu,sqrt(varMu));
		densDelta = dnorm(0,meanDelta,sqrt(varDelta),0);
		
		// sample sig2
		scaleSig2 = 0.5*(sumy2 - 2.0 * N * yBar * mu + (N + 1/g)*pow(mu,2));
		sig2 = 1/rgamma(shapeSig2,1/scaleSig2);
	
		// sample g
		scaleg = 0.5 * (pow(mu,2)/sig2+rscaleSq);
		g = 1/rgamma(1,1/scaleg);
		
		
		chains[npars*i + 0] = mu;
		chains[npars*i + 1] = sig2;
		chains[npars*i + 2] = g;
		chains[npars*i + 3] = mu/sqrt(sig2);
		chains[npars*i + 4] = densDelta;
		
	}

	UNPROTECT(2);
	PutRNGstate();
	
}


SEXP RgibbsEqVariance(SEXP yR, SEXP NR, SEXP JR, SEXP IR, SEXP lambdaR, SEXP iterationsR, SEXP sdMetropSig2R, SEXP sdMetropTauR, SEXP progressR, SEXP pBar, SEXP rho)
{
	int iterations = INTEGER_VALUE(iterationsR);
	int *N = INTEGER_POINTER(NR), progress = INTEGER_VALUE(progressR);
	double lambda = NUMERIC_VALUE(lambdaR);
	double sdMetropSig2 = NUMERIC_VALUE(sdMetropSig2R);
	double sdMetropTau = NUMERIC_VALUE(sdMetropTauR);
	double *y = REAL(yR);
	int J = INTEGER_VALUE(JR),I = INTEGER_VALUE(IR);
	
	int npars = J + 5;
	
	SEXP chainsR;
	PROTECT(chainsR = allocMatrix(REALSXP, npars, iterations));

	gibbsEqVariance(y, N, J, I, lambda, iterations, REAL(chainsR), sdMetropSig2, sdMetropTau, progress, pBar, rho);
	
	UNPROTECT(1);
	
	return(chainsR);
}

void gibbsEqVariance(double *y, int *N, int J, int I, double lambda, int iterations, double *chains, double sdMetropSig2, double sdMetropTau, int progress, SEXP pBar, SEXP rho)
{
	int i=0,j=0,m=0, sumN=0;
	double yBar[J],sumy2[J],logdensg;
	double mu[J],sig2=1,tau=1,SS[J];
	double scaleMu=0,dfMu=0,alpha=0,beta=0;
	int npars = J + 5;
	
	// progress stuff
	SEXP sampCounter, R_fcall;
	int *pSampCounter;
    PROTECT(R_fcall = lang2(pBar, R_NilValue));
	PROTECT(sampCounter = NEW_INTEGER(1));
	pSampCounter = INTEGER_POINTER(sampCounter);
	
	GetRNGstate();
	
	for(j=0;j<J;j++)
	{
		yBar[j]=0;
		sumy2[j]=0;
		sumN+=N[j];
		for(i=0;i<N[j];i++)
		{
			yBar[j] += y[j*I+i]/((double)(N[j]));
			sumy2[j] += pow(y[j*I+i],2);
		}
		SS[j] = sumy2[j] - N[j]*pow(yBar[j],2);
	}
	
	// start MCMC
	for(m=0; m<iterations; m++)
	{
		R_CheckUserInterrupt();
	
		
		//Check progress
		if(progress && !((m+1)%progress)){
			pSampCounter[0]=m+1;
			SETCADR(R_fcall, sampCounter);
			eval(R_fcall, rho); //Update the progress bar
		}


		logdensg = 0;
		// sample mu
		for(j=0;j<J;j++)
		{
			dfMu = N[j] + 2*tau - 1;
			scaleMu  = sqrt((2.0*tau*sig2 + SS[j])/(N[j]*dfMu));
			mu[j] = scaleMu * rt((double)(dfMu)) + yBar[j];
			chains[npars*m + j] = mu[j];
			alpha = N[j]*0.5 + tau;
			beta = (sumy2[j] - 2.0*N[j]*yBar[j]*mu[j] + N[j]*pow(mu[j],2))/(2*sig2) + tau;
			logdensg += alpha*log(beta) - lgammafn(alpha) - beta;
		}		
		chains[npars*m + J] = exp(logdensg);
		

		// sample sig2
		sig2 = sampleSig2EqVar(sig2, mu, tau, yBar, SS, N, sumN, J, sdMetropSig2, &chains[npars*m + J + 2]);
		chains[npars*m + J + 1] = sig2;
	
		// sample tau
		tau = sampleTauEqVar(tau, mu, sig2, yBar, SS, N, J, lambda, sdMetropTau, &chains[npars*m + J + 4]);
		chains[npars*m + J + 3] = tau;	
		
	}

	UNPROTECT(2);
	PutRNGstate();
	
}

double sampleTauEqVar(double tau, double *mu, double sig2, double *yBar, double *SS, int *N, int J, double lambda, double sdMetrop, double *acc)
{
	// we're going to do log(tau) instead.
	double z, b, logratio, logTau = log(tau), cand, expCand;
	
	z = rnorm(0,sdMetrop);
	cand = logTau + z;
	expCand = exp(cand);
	
	logratio = logFullCondTauEqVar(expCand, mu, sig2, yBar, SS, N, J, lambda) - logFullCondTauEqVar(tau, mu, sig2, yBar, SS, N, J, lambda)
			+ (cand - logTau);  
	
	b = log(runif(0,1));
	if(b > logratio){
		acc[0]=0;
		return(tau);
	}else{
		acc[0]=1;
		return(expCand);
	}
}


double sampleSig2EqVar(double sig2, double *mu, double tau, double *yBar, double *SS, int *N, int sumN, int J, double sdMetrop, double *acc)
{
	// we're going to do log(sig2) instead.
	double z, b, logratio, logSig2 = log(sig2), cand, expCand;
	
	z = rnorm(0,sdMetrop);
	cand = logSig2 + z;
	expCand = exp(cand);
	
	logratio = logFullCondSig2EqVar(expCand, mu, tau, yBar, SS, N, sumN, J) - logFullCondSig2EqVar(sig2, mu, tau, yBar, SS, N, sumN, J)
			+ (cand - logSig2);  
	
	
	b = log(runif(0,1));
	if(b > logratio){
		acc[0]=0;
		return(sig2);
	}else{
		acc[0]=1;
		return(expCand);
	}
}


double logFullCondTauEqVar(double tau, double *mu, double sig2, double *yBar, double *SS, int *N, int J, double lambda)
{
	if(tau<=0){ return(-DBL_MAX);}
	double jpart=0, logDens = 0;
	int j=0;
	
	for(j=0;j<J;j++)
	{
		//jpart += lgammafn(N[j]*0.5 + tau) - (N[j]*0.5 + tau)*(log(tau + (SS[j]/sig2)/2) + log(1+pow(mu[j]-yBar[j],2)/(2*tau*sig2/(1.0*N[j]) + SS[j]/N[j])));
		jpart += lgammafn(N[j]*0.5 + tau) - (N[j]*0.5 + tau)*log((SS[j] + N[j]*pow(mu[j]-yBar[j],2))/(2*sig2) + tau);
	}
	logDens = jpart - lambda*tau + J*tau*log(tau) - J*lgammafn(tau);
	
	return(logDens);
}

double logFullCondSig2EqVar(double sig2, double *mu, double tau, double *yBar, double *SS, int *N, int sumN, int J)
{
	if(sig2<=0){ return(-DBL_MAX);}
	double jpart=0, logDens = 0;
	int j=0;
	
	
	for(j=0;j<J;j++)
	{
		//jpart += -(N[j]*0.5 + tau)*(log(tau + (SS[j]/sig2)/2) + log(1+pow(mu[j]-yBar[j],2)/(2*tau*sig2/(1.0*N[j]) + SS[j]/N[j])));
		jpart += - (N[j]*0.5 + tau)*log((SS[j] + N[j]*pow(mu[j]-yBar[j],2))/(2*sig2) + tau);
	}
	logDens = jpart + -(0.5*sumN+1)*log(sig2);
	
	
	return(logDens);
}



SEXP RgibbsOneWayAnova(SEXP yR, SEXP XR, SEXP NR, SEXP JR, SEXP IR, SEXP rscaleR, SEXP iterationsR, SEXP progressR, SEXP pBar, SEXP rho)
{
	int iterations = INTEGER_VALUE(iterationsR);
	int *N = INTEGER_POINTER(NR), progress = INTEGER_VALUE(progressR);
	double rscale = NUMERIC_VALUE(rscaleR);
	double *y = REAL(yR);
	int J = INTEGER_VALUE(JR),I = INTEGER_VALUE(IR);
	int j=0,sumN=0,counter=0;

	npars = J+5;
	
	for(j=0;j<J;j++){
		sumN += N[j];
	}

	double X[sumN * (J+1)]; 
	double yVec[sumN];
	
	int whichJ[sumN]; //which j the element belongs to
	
	for(j=0;j<J;j++)
	{
		X[counter]=1;
		for(i=0;i<N[j];i++){
			X[(j+1)*sumN + counter] = 1;
			whichJ[counter] = j;
			yVec[counter]=y[j*I + i];
			counter++;
		}
	}
	
	SEXP chainsR;
	PROTECT(chainsR = allocMatrix(REALSXP, npars, iterations));

	gibbsOneWayAnova(y, X, N, J, sumN, whichJ, rscale, iterations, REAL(chainsR), progress, pBar, rho);
	
	UNPROTECT(1);
	
	return(chainsR);
}

void gibbsOneWayAnova(double *y, int *X, int *N, int J, int sumN, int *whichJ, double rscale, int iterations, double *chains, int progress, SEXP pBar, SEXP rho)
{
	int i=0,j=0,m=0,Jp1sq = Jp1sq*Jp1sq,Jsq=J*J,Jp1=J+1;
	double ySum[J],sumy2[J],densDelta=0;
	double sig2=1,g=1;
	double B[Jp1sq], double XtX[Jp1sq], B2[Jsq], ZtZ[Jsq];
	double Btemp[Jp1sq],B2temp[Jsq];
	double muTemp[J],oneOverSig2temp=0;
	double beta[J+1],grandSum=0,grandSumSq=0;
	double shapeSig2 = (sumN+J)/2, shapeg = (J+1)/2;
	double rateSig2=0, rateg=0;
	double Xty[J+1];
	double logDet=0;
	double rscaleSq=rscale*rscale;
	
	int iOne=1;
	double dOne=1,dZero=0;
	
	// progress stuff
	SEXP sampCounter, R_fcall;
	int *pSampCounter;
    PROTECT(R_fcall = lang2(pBar, R_NilValue));
	PROTECT(sampCounter = NEW_INTEGER(1));
	pSampCounter = INTEGER_POINTER(sampCounter);
	
	npars=J+5;
	
	GetRNGstate();
	AZERO(B,pow(J+1,2));
	AZERO(B2,pow(J,2));
	
	XtX[0]=sumN;

	for(i=0;i<sumN;i++)
	{
		j = whichJ[i];
		ySum[j] += y[i];
		sumy2[j] += pow(y[i],2);
	}
	
	
	for(j=0;j<J;j++)
	{
		ySum[j]=0;
		sumy2[j]=0;
		grandSum += ySum[j];
		grandSumSq += sumy2[j];
		XtX[j+1]=N[j];
		XtX[(J+1)*(j+1)]=N[j];
		for(i=0;i<J;i++)
		{
			if(i==j)
			{
				XtX[(i+1)*J+(j+1)] = N[j];  
				ZtZ[i*J + j] = N[j];
			}else{
				XtX[(i+1)*J+(j+1)] = 0;
				XtX[i*J + j] = 0;
			}
			
		}
	}
	
	Xty[0] = grandSum;
	Memcpy(&Xty[1],ySum,J+1);
	beta[0]=0;
	AZER0(&beta[1],J);
	
	
	// start MCMC
	for(m=0; m<iterations; m++)
	{
		R_CheckUserInterrupt();
	
		//Check progress
		if(progress && !((m+1)%progress)){
			pSampCounter[0]=m+1;
			SETCADR(R_fcall, sampCounter);
			eval(R_fcall, rho); //Update the progress bar
		}

		// sample beta
		Memcpy(Btemp,XtX,Jp1sq);
		
		for(j=0;j<J;j++){
			Btemp[(j+1)*(J+1)+(j+1)] += 1/g;
		}
		InvMatrixUpper(Btemp, J+1);
		internal_symmetrize(Btemp,J+1);
		for(j=0;j<Jp1sq;j++)
			Btemp[j] *= sig2;
		oneOverSig2temp = 1/sig2;
		F77_CALL(dsymv)("U", &(Jp1), &oneOverSig2temp, Btemp, &(Jp1), Xty, &iOne, &dZero, beta, &iOne);
		rmvGaussianC(beta, Btemp, J+1);
		Memcpy(&chains[npars*m],beta,J+1);	
		
		
		// calculate density (double Standardized)
		Memcpy(B2temp,ZtZ,Jsq);
		densDelta = -J*0.5*log(2*M_PI);
		for(j=0;j<J;j++)
		{
			B2temp[j*J+j] += 1/g;
			muTemp[j] = ySum[j]-N[j]*beta[0];
		}
		for(j=0;j<Jsq;j++)
			B2temp[j] *= g;
		InvMatrixUpper(B2temp, J);
		internal_symmetrize(B2temp,J);
		densDelta += -0.5*matrixDet(B2temp,J,1);
		densDelta += quadform(muTemp, B2temp, J, 1);
		chains[npars*m + (J+1) + 0] = exp(densDelta);
		
		// calculate density (single Standardized)
		densDelta = -J*0.5*log(2*M_PI);
		//chains[npars*m + (J+1) + 1] = densDelta;
		
		tempBetaY = 0;
		tempBetaSq= 0;
		// sample sig2
		for(j=0;j<J;j++)
		{
			tempBetaY += beta[j+1]*ySum[j];
			tempBetaSq += pow(beta[j+1],2);
		}
		rateSig2 = 0.5*(sumy2 - 2*tempBetaY + (1+1/g)*tempBetaSq);
		sig2 = rgamma(shapeSig2,1/rateSig2);
		chains[npars*m + (J+1) + 2] = sig2;
	
		// sample g
		rateg = 0.5*(tempBetaSq/sig2) + rscaleSq/2;
		g = 1/rgamma(shapeg,1/rateg);
		chains[npars*m + (J+1) + 3] = g;
	}

	UNPROTECT(2);
	PutRNGstate();
	
}


// Compute determinant of an N by N matrix A
double matrixDet(double *A, int N, int log)
{
//SUBROUTINE DPOTRF( UPLO, N, A, LDA, INFO )
	int i=0, info=0;
	double B[N*N], logDet=0;
	
	Memcpy(B,A,N*N);
	
	F77_CALL(dpotrf)("U", &N, B, &N, &info);
	if(info){
		Rprintf("Cholesky decomposition in matrixDet() returned nonzero info %d.\n",info);
	}
	
	for(i=0;i<N;i++)
	{
		logDet += 2 * log(B[i*N+i]);
	}
	
	if(log){
		return(logDet);
	}else{
		return(exp(logDet));
	}
}

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


double quadform(double *x, double *A, int N, int incx)
{
  
  int Nsqr = N*N,info,i=0;
  double *B = Calloc(Nsqr,double);
  //double one=1;
  //double zero=0;
  double sumSq=0;
  double y[N];
  int iOne=1;

  for(i=0;i<N;i++){
    y[i] = x[i*incx];
    //printf("%d %f\n",i,y[i]);
  }
  Memcpy(B,A,Nsqr);
  
  F77_NAME(dpotrf)("U", &N, B, &N, &info);
  F77_NAME(dtrmv)("U","N","N", &N, B, &N, y, &iOne);
  
  for(i=0;i<N;i++){
    sumSq += y[i]*y[i];
  }
  
  Free(B);
  
  return(sumSq);
}

int InvMatrixUpper(double *A, int p)
{
      int info1, info2;
      F77_NAME(dpotrf)("U", &p, A, &p, &info1);
      F77_NAME(dpotri)("U", &p, A, &p, &info2);      
      //make sure you make it symmetric later...
      return(info1);
}


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

