
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
void gibbsEqVariance(double *y, int *N, int J, int I, double lambda, int iterations, double *chains, double sdMetropSig2, double sdMetropTau, int progress, SEXP pBar, SEXP rho);
double sampleSig2EqVar(double sig2, double *mu, double tau, double *yBar, double *SS, int *N, int sumN, int J, double sdMetrop, double *acc);
double logFullCondTauEqVar(double tau, double *mu, double sig2, double *yBar, double *SS, int *N, int J, double lambda);
double logFullCondSig2EqVar(double sig2, double *mu, double tau, double *yBar, double *SS, int *N, int sumN, int J);
double sampleTauEqVar(double tau, double *mu, double sig2, double *yBar, double *SS, int *N, int J, double lambda, double sdMetrop, double *acc);

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
