
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



SEXP RgibbsOneSample(SEXP yR, SEXP NR, SEXP rscaleR, SEXP iterationsR)
{
	
	int npars = 3,iterations = INTEGER_VALUE(iterationsR);
	int N = INTEGER_VALUE(NR);
	double rscale = NUMERIC_VALUE(rscaleR);
	double *y = REAL(yR);
	
	SEXP chainsR;
	PROTECT(chainsR = allocMatrix(REALSXP, npars, iterations));

	gibbsOneSample(y, N, rscale, iterations, REAL(chainsR));
	
	return(chainsR);
}

void gibbsOneSample(double *y, int N, double rscale, int iterations, double *chains)
{
	int i=0;
	double yBar=0,sumy2=0,rscaleSq=pow(rscale,2);
	double mu=0,sig2=1,g=1;
	double meanMu=0, varMu=0, shapeSig2=(1.0*N)/2+0.5,scaleSig2=1,scaleg=1;
	int npars = 3;
	
	
	for(i=0;i<N;i++)
	{
		yBar += y[i]/((float)(N));
		sumy2 += pow(y[i],2);
	}
	
	// start MCMC
	for(i=0; i<iterations; i++)
	{
		
		// sample mu
		meanMu = yBar * N/(1.0*N+1/g);
		varMu  = sig2/(1.0*N + 1/g);
		mu = rnorm(meanMu,sqrt(varMu));
		
		// sample sig2
		scaleSig2 = 0.5*(sumy2 - 2.0 * N * yBar + 1.0*N*pow(mu,2));
		sig2 = 1/rgamma(shapeSig2,1/scaleSig2);
		
		// sample g
		scaleg = 0.5 * (pow(mu,2)/sig2+rscaleSq);
		g = 1/rgamma(1,1/scaleg);
		
		chains[npars*i + 0] = mu;
		chains[npars*i + 1] = sig2;
		chains[npars*i + 2] = g;
	}

}
