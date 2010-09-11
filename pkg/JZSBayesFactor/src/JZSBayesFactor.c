
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

void gibbsOneSample(double *y, int N, double rscale, int iterations, double *chains, int progress, SEXP pBar, SEXP rho);
void gibbsEqVariance(double *y, int *N, int J, int I, double lambda, int iterations, double *chains, double sdMetrop, int progress, SEXP pBar, SEXP rho);


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
		scaleSig2 = 0.5*(sumy2 - 2.0 * N * yBar * mu + 1.0*N*pow(mu,2));
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


SEXP RgibbsEqVariance(SEXP yR, SEXP NR, SEXP JR, SEXP IR, SEXP lambdaR, SEXP iterationsR, SEXP sdMetropR, SEXP progressR, SEXP pBar, SEXP rho)
{
	
	int iterations = INTEGER_VALUE(iterationsR);
	int *N = INTEGER_POINTER(NR), progress = INTEGER_VALUE(progressR);
	double lambda = NUMERIC_VALUE(rscaleR);
	double sdMetrop = NUMBERIC_VALUE(sdMetropR);
	double *y = REAL(yR);
	int J = INTEGER_VALUE(JR),I = INTEGER_VALUE(IR);
	
	int npars = 2*J + 3;
	
	SEXP chainsR;
	PROTECT(chainsR = allocMatrix(REALSXP, npars, iterations));

	gibbsEqVariance(y, N, J, I, lambda, iterations, REAL(chainsR), sdMetrop, progress, pBar, rho);
	
	UNPROTECT(1);
	
	return(chainsR);
}

void gibbsEqVariance(double *y, int *N, int J, int I, double lambda, int iterations, double *chains, double sdMetrop, int progress, SEXP pBar, SEXP rho)
{
	int i=0,j=0,m=0;
	double yBar[J],sumy2[J],densg[J];
	double mu[J],sig2=1,g[J],tau=1;
	double varMu=0, shapeSig2=(1.0*N)/2+0.5,scaleSig2=1,scaleg=1;
	int npars = 2*J + 3;
	
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
		for(i=0;i<N[j];i++)
		{
			yBar[j] += y[j*I+i]/((float)(N[j]));
			sumy2[j] += pow(y[j*I+i],2);
		}
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


		// sample mu
		for(j=0;j<J;j++)
		{
			varMu  = sig2/(N[j]*g[j]);
			mu[j] = rnorm(yBar,sqrt(varMu));
			chains[npars*m + j] = mu[j];
		}
				
		
		logdensg = 0;
		// sample g
		for(j=0;j<J;j++)
		{
			//scaleg = 0.5 * (pow(mu,2)/sig2+rscaleSq);
			//shapeg = 
			g[j] = rgamma(shapeg,scaleg);
			chains[npars*m + J + j] = g[j];
			//logdensg +=  
		}
		chains[npars*m + J + J] = exp(logdensg);
		
		
		// sample sig2
		//scaleSig2 = 0.5*(sumy2 - 2.0 * N * yBar * mu + 1.0*N*pow(mu,2));
		//sig2 = 1/rgamma(shapeSig2,1/scaleSig2);
		chains[npars*m + J + J + 1] = sig2;
	
		// sample tau
		// tau = //MUST USE METROPOLIS-HASTINGS
		chains[npars*m + J + J + 2] = tau;	
		
	}

	UNPROTECT(2);
	PutRNGstate();
	
}
