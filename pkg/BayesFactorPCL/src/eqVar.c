#include "BFPCL.h"


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



SEXP RgibbsEqVarianceM2(SEXP yR, SEXP NR, SEXP JR, SEXP IR, SEXP alphaR, SEXP betaR, SEXP iterationsR, SEXP sdMetropgR, SEXP sdDecorrR, SEXP newtonStepsR, SEXP progressR, SEXP pBar, SEXP rho)
{
	int iterations = INTEGER_VALUE(iterationsR);
	int *N = INTEGER_POINTER(NR), progress = INTEGER_VALUE(progressR);
	double alpha = NUMERIC_VALUE(alphaR);
	double beta = NUMERIC_VALUE(betaR);
	double sdMetropg = NUMERIC_VALUE(sdMetropgR);
	double sdDecorr = NUMERIC_VALUE(sdDecorrR);
	double *y = REAL(yR);
	int J = INTEGER_VALUE(JR),I = INTEGER_VALUE(IR);
	
	int newtonSteps = INTEGER_VALUE(newtonStepsR);	
	int npars = 3*J + 4;
	
	SEXP chainsR, debugR, returnList;
	PROTECT(chainsR = allocMatrix(REALSXP, npars, iterations));
	PROTECT(debugR = allocMatrix(REALSXP, J*5, iterations));
	PROTECT(returnList = allocVector(VECSXP, 2));
	
	gibbsEqVarianceM2(y, N, J, I, alpha, beta, iterations, REAL(chainsR), REAL(debugR), sdMetropg, sdDecorr, newtonSteps, progress, pBar, rho);
	
	SET_VECTOR_ELT(returnList, 0, chainsR);
    SET_VECTOR_ELT(returnList, 1, debugR);
	
	UNPROTECT(3);
	
	return(returnList);
}

void gibbsEqVarianceM2(double *y, int *N, int J, int I, double alpha, double beta, int iterations, double *chains, double *debug, double sdMetropg, double sdDecorr, int newtonSteps, int progress, SEXP pBar, SEXP rho)
{
	int i=0,j=0,m=0, sumN=0;
	double yBar[J],sumy2[J],IWMDE=0;
	double g[J],mu[J],sig2=1,sig2g=1,SS[J];
	double modeIWMDE = 0, sdIWMDE = 0;
	double sumSS=0, logVar[J], sumVar=0;
	int npars = 3*J + 4;
	
	double shapeSig2 = 0, shapeSig2g = 0;
	double scaleSig2 = 0, scaleSig2g = 0;
	double sumg2=0, sumg=0;
	double c1=0,c2=0,c3=0;
	
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
		g[j]=0;
		for(i=0;i<N[j];i++)
		{
			yBar[j] += y[j*I+i]/((double)(N[j]));
			sumy2[j] += pow(y[j*I+i],2);
		}
		SS[j] = sumy2[j] - N[j]*pow(yBar[j],2);
		mu[j] = yBar[j];
		sumSS += SS[j];
		logVar[j] = log(SS[j]/N[j]);
		sumVar += logVar[j];
	}
	sig2 = sumSS/(sumN-J);
	shapeSig2 = 0.5*sumN;
	shapeSig2g = 0.5*J + alpha;
	
	
	for(j=0;j<J;j++)
		logVar[j] = logVar[j] - sumVar/J;
	
	
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


		// sample mu and g
		IWMDE=0;
		sumg2=0;
		sumg=0;
		scaleSig2=0;
		for(j=0;j<J;j++)
		{
			//sigma
			//mu[j] = rnorm(yBar[j],exp(g[j])*sqrt(sig2/N[j]));
			//sigma2
			mu[j] = rnorm(yBar[j],exp(0.5*g[j])*sqrt(sig2/N[j]));
			g[j] = samplegEqVarM2(g[j], sig2, mu[j], sig2g, yBar[j], sumy2[j], N[j], sdMetropg, &chains[npars*m + 2*J + 4 + j]);
			sumg2 += g[j]*g[j];	
			sumg += g[j];
		
			//sigma
			//scaleSig2 += (sumy2[j] - 2*mu[j]*N[j]*yBar[j] + N[j]*mu[j]*mu[j])/(2*exp(2*g[j]));
			//sigma2
			scaleSig2 += (sumy2[j] - 2*mu[j]*N[j]*yBar[j] + N[j]*mu[j]*mu[j])/(2*exp(g[j]));
			
			// IWMDE
			c1 = -1/(2*sig2g);
			c2 = 0.5*N[j]*sig2g;
			c3 = -(sumy2[j] - 2*mu[j]*N[j]*yBar[j] + N[j]*mu[j]*mu[j])/(2*sig2);
			
			//modeIWMDE = newtonMethodEqVar2(logVar[j],c1,c2,c3,1,newtonSteps);
			modeIWMDE = newtonMethodEqVar2(g[j],c1,c2,c3,1,newtonSteps);
			sdIWMDE = 1/sqrt(-2*c1 - c3*exp(-modeIWMDE));
			IWMDE += dnorm(g[j],modeIWMDE,sdIWMDE,1) + -c3*(exp(-g[j]) - 1) + 0.5*N[j]*g[j] + g[j]*g[j]/(2*sig2g);
			debug[m*(2*J) + 2*j + 0] = modeIWMDE;
			debug[m*(2*J) + 2*j + 1] = sdIWMDE;
			debug[m*(2*J) + 2*j + 2] = c1;
			debug[m*(2*J) + 2*j + 3] = c2;
			debug[m*(2*J) + 2*j + 4] = c3;
		}		
		
		// sample sig2
		sig2 = 1/rgamma(shapeSig2,1/scaleSig2);
	
		// sample sig2g
		scaleSig2g = sumg2/2 + beta;
		sig2g = 1/rgamma(shapeSig2g,1/scaleSig2g);
	
		// Decorelating step
		decorrStepEqVarM2(sumg, sig2g, J, &g[0], &sig2, sdDecorr, &chains[npars*m + 2*J + 3]);
		
		// copy to chains
		Memcpy(chains + npars*m, mu, J);
		Memcpy(chains + npars*m + J, g, J);
		chains[npars*m + 2*J + 0] = IWMDE; 
		chains[npars*m + 2*J + 1] = sig2; 
		chains[npars*m + 2*J + 2] = sig2g; 
	}

	UNPROTECT(2);
	PutRNGstate();
	
}


void decorrStepEqVarM2(double sumg, double sig2g, int J, double *g, double *sig2, double sdDecorr, double *acc)
{
	double z, logratio,b=0;
	int j=0;
	
	z = rnorm(0, sdDecorr);
	logratio = -(J*z*z - 2*z*sumg)/(2*sig2g);
	
	
	b = log(runif(0,1));
	if(b > logratio){
		acc[0]=0;
	}else{
		acc[0]=1;
		for(j=0;j<J;j++)
			g[j] = g[j]-z;
			sig2[0] = sig2[0]*exp(z);
	}
}

double newtonMethodEqVar2(double xt, double c1, double c2, double c3, int iterations, int max)
{
	
	xt =	xt - (2*c1*(xt + c2) - c3*exp(-xt))/(2*c1 + c3*exp(-xt));
	if(iterations==max){
		return(xt);
	}else{
		return(newtonMethodEqVar2(xt,c1,c2,c3,iterations+1,max));
	}
}


double samplegEqVarM2(double g, double sig2, double mu, double sig2g, double yBar, double sumy2, int N, double sdMetrop, double *acc)
{
	double z, b, logratio, cand;
	
	z = rnorm(0,sdMetrop);
	cand = g + z;
	
	logratio = fullCondgEqVarM2(cand, sig2, mu, sig2g, yBar, sumy2, N) - fullCondgEqVarM2(g, sig2, mu, sig2g, yBar, sumy2, N); 
	
	b = log(runif(0,1));
	if(b > logratio){
		acc[0]=0;
		return(g);
	}else{
		acc[0]=1;
		return(cand);
	}
}

double fullCondgEqVarM2(double g, double sig2, double mu, double sig2g, double yBar, double sumy2, int N)
{
	//sigma
	//return(-0.5/sig2g * pow(g + N*sig2g,2) - exp(-2*g)*(sumy2 - 2*mu*N*yBar + N*mu*mu)/(2*sig2));
	//sigma2
	return(-0.5/sig2g * pow(g + 0.5*N*sig2g,2) - exp(-g)*(sumy2 - 2*mu*N*yBar + N*mu*mu)/(2*sig2));
}

