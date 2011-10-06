#include "BFPCL.h"

SEXP RgibbsTwoSampleARmulti(SEXP yR, SEXP NR, SEXP NobsR, SEXP tR, SEXP rscaleR, SEXP betaThetaR, SEXP iterationsR, SEXP sdmetR, SEXP progressR, SEXP pBar, SEXP rho)
{
	int iterations = INTEGER_VALUE(iterationsR);
	int N = INTEGER_VALUE(NR), progress = INTEGER_VALUE(progressR);
	double rscale = NUMERIC_VALUE(rscaleR);
	double betaTheta = NUMERIC_VALUE(betaThetaR);
	double sdmet = NUMERIC_VALUE(sdmetR);
	double *y = REAL(yR);
	double *t = REAL(tR);
	int *Nobs = INTEGER(NobsR);
	
	int npars = 3*N + 4;
	
	SEXP chainsR;
	PROTECT(chainsR = allocMatrix(REALSXP, npars, iterations));

	gibbsTwoSampleARmulti(y, N, Nobs, t, rscale, betaTheta, iterations, sdmet, REAL(chainsR), progress, pBar, rho);
	
	UNPROTECT(1);
	
	return(chainsR);
}


void gibbsTwoSampleARmulti(double *y, int N, int *Nobs, double *t, double rscale,  double betaTheta, int iterations, double sdmet, double *chains, int progress, SEXP pBar, SEXP rho)
{
	int npars = 3*N + 4;
	int i=0, j=0, k=0, m=0,maxobs=0,iOne=1,cumobs[N];
	double SSq,varMu, meanMu, varBeta, meanBeta, aSig2, bSig2, ag, bg, rscsq=rscale*rscale;
	double mu0[N], sig2[N], beta[N], theta = 0, g = 1;
	double muBeta=0,sig2Beta=1;
	double tOnePsiOne,dZero=0,dOne=1;
	double ttPsit,loglikeTheta;
	
	for(i=0;i<N;i++)
	{
		sig2[i] = 1;
		beta[i] = 0;
		if(Nobs[i]>maxobs) maxobs = Nobs[i];
		if(i==0)
		{
			cumobs[i] = Nobs[i];
		}else{
			cumobs[i] = cumobs[i-1] + Nobs[i];
		}
		
	}

	double tempV[maxobs];
	double ones[maxobs];
	double psiOne[maxobs], psit[maxobs];
	double invPsi[maxobs*maxobs];
	AZERO(invPsi,maxobs*maxobs);
	
	for(i=0;i<maxobs;i++)
	{
		ones[i] = 1;		
	}
	

	// progress stuff
	SEXP sampCounter, R_fcall;
	int *pSampCounter;
    PROTECT(R_fcall = lang2(pBar, R_NilValue));
	PROTECT(sampCounter = NEW_INTEGER(1));
	pSampCounter = INTEGER_POINTER(sampCounter);

	GetRNGstate();

	for(m=0;m<iterations;m++)
	{

		R_CheckUserInterrupt();
	
		//Check progress
		if(progress && !((m+1)%progress)){
			pSampCounter[0]=m+1;
			SETCADR(R_fcall, sampCounter);
			eval(R_fcall, rho); //Update the progress bar
		}

		for(i=0;i<maxobs;i++)
		{
			invPsi[i + maxobs*i] = 1/(1-theta*theta);		
			for(j=0;j<i;j++)
			{
				invPsi[i + maxobs*j] = invPsi[i + maxobs*i] * pow(theta,abs(i-j));
				invPsi[j + maxobs*i] = invPsi[i + maxobs*j];
			}
		}
		
		InvMatrixUpper(invPsi, maxobs);
		internal_symmetrize(invPsi, maxobs);
	

		//mu0
		for(i=0;i<N;i++)
		{
			tOnePsiOne = quadform(ones, invPsi, Nobs[i], 1, maxobs);
			F77_NAME(dsymv)("U", &Nobs[i], &dOne, invPsi, &maxobs, ones, &iOne, &dZero, psiOne, &iOne);

			meanMu = 0;
			varMu = sig2[i]/tOnePsiOne;
			for(j=0;j<Nobs[i];j++)
			{
				k = j + cumobs[i]-Nobs[i];
				meanMu += (y[k] - beta[i]*t[k])*psiOne[j];
			}
			meanMu = meanMu/tOnePsiOne;
			//Rprintf("%d %d: meanMu:%f varMu:%f\n",m,i,meanMu,varMu);
			mu0[i] = rnorm(meanMu,sqrt(varMu));
		}
		
		//beta
		for(i=0;i<N;i++)
		{
			k = cumobs[i]-Nobs[i];
			ttPsit = quadform(t + k, invPsi, Nobs[i], 1, maxobs);
			F77_NAME(dsymv)("U", &Nobs[i], &dOne, invPsi, &maxobs, t + k, &iOne, &dZero, psit, &iOne);
			

			meanBeta = 0;
			varBeta = 1/(ttPsit/sig2[i] + 1/sig2Beta);
			for(j=0;j<Nobs[i];j++)
			{
				k = j + cumobs[i]-Nobs[i];
				meanBeta += (y[k] - mu0[i])*psit[j]/sig2[i];
			}
			meanBeta = (meanBeta + muBeta/sig2Beta)/(ttPsit/sig2[i] + 1/sig2Beta);
			beta[i] = rnorm(meanBeta, sqrt(varBeta));
		}
			
		
		//muBeta
		// reuse meanBeta and varBeta, but for sums here
		meanBeta=0;
		varBeta=0;
		for(i=0;i<N;i++)
		{
			meanBeta += beta[i];
			varBeta += beta[i]*beta[i];
		}
		muBeta = rnorm(meanBeta/(N+1/g), sqrt(sig2Beta/(N+1/g)));
	
		//sig2Beta
		SSq = varBeta - 2*muBeta*meanBeta + N*muBeta*muBeta + muBeta*muBeta/g;
		sig2Beta = 1/rgamma((N+1)*0.5, 2/SSq);
		
		//sig2
		for(i=0;i<N;i++)
		{
			aSig2 = 0.5*Nobs[i];
			for(j=0;j<Nobs[i];j++)
			{
				k = j + cumobs[i]-Nobs[i];
				tempV[j] = (y[k] - mu0[i] - beta[i]*t[k]);
			}
			bSig2 = 0.5*quadform(tempV,invPsi,Nobs[i], 1, maxobs);
			sig2[i] = 1/rgamma(aSig2,1/bSig2);
		}
		
		//g
		ag = 1;
		bg = 0.5*(muBeta*muBeta/sig2Beta + rscsq);
		g = 1/rgamma(ag,1/bg);
		
		//theta
		theta = sampThetaARmulti(theta, mu0, beta, sig2, g, y, N, t, Nobs, cumobs, maxobs, betaTheta, sdmet);
	
		// write chain
		
		Memcpy(&chains[npars*m + 0], mu0, N);
		Memcpy(&chains[npars*m + N], beta, N);
		Memcpy(&chains[npars*m + N + N], &muBeta, 1);
		Memcpy(&chains[npars*m + N + N + 1], &sig2Beta, 1);
		Memcpy(&chains[npars*m + N + N + 2], sig2, N);
		Memcpy(&chains[npars*m + N + N + 2 + N], &g, 1);
		Memcpy(&chains[npars*m + N + N + 2 + N + 1], &theta, 1);
	}

	UNPROTECT(2);
	PutRNGstate();
}


double sampThetaARmulti(double theta, double *mu0, double *beta, double *sig2, double g, double *y, int N, double *t, int *Nobs, int *cumobs, int maxobs, double betaTheta , double sdmet)
{
	// sample theta with Metropolis-Hastings
	double cand, likeRat, b;
	
	//cand = rnorm(theta,sdmet);
	cand = theta + sdmet*norm_rand();
	//Rprintf("prev %f cand %f\n",theta,cand);
	
	if(cand<0 || cand>1)
	{
		return(theta);
	}
	//b = log(runif(0,1));
	b = log(unif_rand());
	likeRat = thetaLogLikeARmulti(cand, mu0, beta, sig2, g, y, N, t, Nobs, cumobs, maxobs, betaTheta) - thetaLogLikeARmulti(theta, mu0, beta, sig2, g, y, N, t, Nobs, cumobs, maxobs, betaTheta);
	
	if(b>likeRat){
		return(theta);
	}else{
		return(cand);
	}
}

double thetaLogLikeARmulti(double theta, double *mu0, double *beta, double *sig2, double g, double *y, int N, double *t, int *Nobs, int *cumobs, int maxobs, double betaTheta)
{
	int i,j,k;
	double loglike=0,tempV[maxobs],invPsi[maxobs*maxobs];
	
	for(i=0;i<maxobs;i++)
	{
		invPsi[i + maxobs*i] = 1/(1-pow(theta,2));		
		
		for(j=0;j<i;j++)
		{
			invPsi[i + maxobs*j] = invPsi[i + maxobs*i] * pow(theta,abs(i-j));
			invPsi[j + maxobs*i] = invPsi[i + maxobs*j];
		}
	}

	InvMatrixUpper(invPsi, maxobs);
	internal_symmetrize(invPsi, maxobs);
	
	for(i=0;i<N;i++)
	{
		for(j=0;j<Nobs[i];j++)
		{
			k = j + cumobs[i]-Nobs[i];
			tempV[j] = y[k] - mu0[i] - beta[i]*t[k];
		}
		loglike +=  0.5 * matrixDet(invPsi, Nobs[i], maxobs, 1) - 
		0.5/sig2[i] * quadform(tempV,invPsi,Nobs[i], 1, maxobs);
	}
	
	loglike +=  (betaTheta-1)*log(1-theta);
	
	//Rprintf("theta: %f, loglike: %f\n",theta, loglike);

	return(loglike);
}


SEXP RthetaLogLikeARmulti(SEXP thetaR, SEXP mu0R, SEXP betaR, SEXP sig2R, SEXP gR, SEXP yR, SEXP NR, SEXP tR, SEXP NobsR, SEXP cumobsR, SEXP maxobsR, SEXP betaThetaR)
{
	int N = INTEGER_VALUE(NR), maxobs = INTEGER_VALUE(maxobsR);
	int *Nobs = INTEGER(NobsR), *cumobs=INTEGER(cumobsR);
	double theta = NUMERIC_VALUE(thetaR);
	double g = NUMERIC_VALUE(gR);
	double *mu0 = REAL(mu0R);
	double *beta = REAL(betaR);
	double *sig2 = REAL(sig2R);
	double *y = REAL(yR);
	double *t = REAL(tR);
	double betaTheta = NUMERIC_VALUE(betaThetaR);
	
	double loglike = thetaLogLikeARmulti(theta, mu0, beta, sig2, g, y, N, t, Nobs, cumobs, maxobs, betaTheta);
	
	
	SEXP ret;
	PROTECT(ret = allocVector(REALSXP, 1));
	REAL(ret)[0] = loglike;
	
	UNPROTECT(1);
	return(ret);
	
}
