#include "BFPCL.h"



SEXP RgibbsOneSample(SEXP yR, SEXP NR, SEXP rscaleR, SEXP iterationsR, SEXP doIntervalR, SEXP intervalR, SEXP progressR, SEXP pBar, SEXP callback, SEXP rho)
{
	
	int npars = 6,iterations = INTEGER_VALUE(iterationsR);
	int N = INTEGER_VALUE(NR), progress = INTEGER_VALUE(progressR);
	double rscale = NUMERIC_VALUE(rscaleR);
	double *y = REAL(yR);
	double *interval = REAL(intervalR);
	int doInterval = INTEGER_VALUE(doIntervalR);
	
	SEXP chainsR;
	PROTECT(chainsR = allocMatrix(REALSXP, npars, iterations));

	gibbsOneSample(y, N, rscale, iterations, REAL(chainsR), doInterval, interval, progress, pBar, callback, rho);
	
	UNPROTECT(1);
	
	return(chainsR);
}

void gibbsOneSample(double *y, int N, double rscale, int iterations, double *chains, int doInterval, double *interval, int progress, SEXP pBar, SEXP callback, SEXP rho)
{
	int i=0;
	double yBar=0,sumy2=0,rscaleSq=pow(rscale,2),densDelta=0,meanDelta=0,varDelta=1;
	double mu=0,sig2=1,g=1;
	double meanMu=0, varMu=0, shapeSig2=(1.0*N)/2+0.5,scaleSig2=1,scaleg=1;
	int npars = 6;
	double deltaArea;
	
	// progress stuff
	SEXP sampCounter, R_fcall;
	int *pSampCounter;
    PROTECT(R_fcall = lang2(pBar, R_NilValue));
	PROTECT(sampCounter = NEW_INTEGER(1));
	pSampCounter = INTEGER_POINTER(sampCounter);
	
  //callback stuff
  SEXP callbackCounter, R_fcall2, callbackResult;
	int *pCallbackCounter;
    PROTECT(R_fcall2 = lang2(callback, R_NilValue));
	PROTECT(callbackCounter = NEW_INTEGER(1));
  PROTECT(callbackResult = NEW_INTEGER(1));
	pCallbackCounter = INTEGER_POINTER(callbackCounter);
  
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
	  //Do callback
    pCallbackCounter[0] = ( 1000 * (i+1) ) / iterations;
	  SETCADR(R_fcall2, callbackCounter);
    callbackResult = eval(R_fcall2, rho);
    if(INTEGER(callbackResult)[0]){
        error("Operation cancelled: code %d",INTEGER(callbackResult)[0]);
    }
    
		// sample mu
		meanMu = yBar * N/(1.0*N+1/g);
		meanDelta = meanMu/sqrt(sig2);
		varDelta = 1/(1.0*N + 1/g);
		varMu  = sig2*varDelta;
		mu = rnorm(meanMu,sqrt(varMu));
		densDelta = dnorm(0,meanDelta,sqrt(varDelta),0);
		
		if(doInterval)
		{
			deltaArea = pnorm(interval[1],meanDelta,sqrt(varDelta),1,0) - pnorm(interval[0],meanDelta,sqrt(varDelta),1,0);
		}else{
			deltaArea = 1;
		}
		
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
		chains[npars*i + 5] = deltaArea;
		
	}

	UNPROTECT(5);
	PutRNGstate();
	
}

