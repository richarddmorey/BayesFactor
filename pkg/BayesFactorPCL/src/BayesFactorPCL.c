
#include "BFPCL.h"


void debugPrintMatrix(double *X, int rows, int cols)
{
	int i=0,j=0;
	
	for(i=0;i<rows;i++){
		for(j=0;j<cols;j++){
			Rprintf("%f ",X[j*rows+i]);
		}
		Rprintf("\n");
	}
}

void debugPrintVector(double *x, int len)
{
	int i=0;
	
	for(i=0;i<len;i++){
		Rprintf("%f ",x[i]);
	}
	Rprintf("\n");
}

/**
 * Allocate a 3-dimensional array
 *
 * @param mode The R mode (e.g. INTSXP)
 * @param nrow number of rows
 * @param ncol number of columns
 * @param nface number of faces
 *
 * @return A 3-dimensional array of the indicated dimensions and mode
 */
SEXP alloc3Darray(SEXPTYPE mode, int nrow, int ncol, int nface)
{
    SEXP s, t;
    int n;

    if (nrow < 0 || ncol < 0 || nface < 0)
	error("negative extents to 3D array");
    if ((double)nrow * (double)ncol * (double)nface > INT_MAX)
	error("alloc3Darray: too many elements specified");
    n = nrow * ncol * nface;
    PROTECT(s = allocVector(mode, n));
    PROTECT(t = allocVector(INTSXP, 3));
    INTEGER(t)[0] = nrow;
    INTEGER(t)[1] = ncol;
    INTEGER(t)[2] = nface;
    setAttrib(s, R_DimSymbol, t);
    UNPROTECT(2);
    return s;
}

SEXP RMCTwoSampleAR(SEXP yR, SEXP NR, SEXP tR, SEXP rscaleR, SEXP alphaThetaR, SEXP betaThetaR, SEXP iterationsR)
{
	int iterations = INTEGER_VALUE(iterationsR);
	int N = INTEGER_VALUE(NR);
	double rscale = NUMERIC_VALUE(rscaleR);
	double alphaTheta = NUMERIC_VALUE(alphaThetaR);
	double betaTheta = NUMERIC_VALUE(betaThetaR);
	double *y = REAL(yR);
	double *t = REAL(tR);
	
	SEXP logBFR;
	PROTECT(logBFR = allocVector(REALSXP, 1));
	double *logBF = REAL(logBFR);

	logBF[0] = MCTwoSampleAR(y, N, t, rscale, alphaTheta, betaTheta, iterations);
	
	UNPROTECT(1);
	
	return(logBFR);
}

double MCTwoSampleAR(double *y, int N, double *t, double rscale, double alphaTheta, double betaTheta, int iterations)
{
	int i=0, j=0, m=0, Nsqr = N*N, iOne=1;
	double theta,g,rscsq=rscale*rscale,mLike;
	double ones[N],dOne=1,dZero=0;
	double invPsi[Nsqr],invPsi0[Nsqr],invPsi1[Nsqr];
	double tempV[N];
	double tOnePsiOne,ttPsi0t;
	double tempS, devs, logmlike[iterations];
	
	AZERO(invPsi,Nsqr);
	
	
	for(i=0;i<N;i++)
	{
		ones[i] = 1;
	}

	GetRNGstate();
	
	for(m=0;m<iterations;m++)
	{
		g = 1/rgamma(0.5,2/rscsq);
		theta = rbeta(alphaTheta,betaTheta);
		
		for(i=0;i<N;i++)
		{
			invPsi[i + N*i] = 1/(1-theta*theta);		
			for(j=0;j<i;j++)
			{
				invPsi[i + N*j] = invPsi[i + N*i] * pow(theta,abs(i-j));
				invPsi[j + N*i] = invPsi[i + N*j];
			}
		}
		
		
		InvMatrixUpper(invPsi, N);
		internal_symmetrize(invPsi, N);

		Memcpy(invPsi0,invPsi,Nsqr);
		
		tOnePsiOne = quadform(ones, invPsi, N, 1, N);
		tempS = -1/tOnePsiOne;
		
		F77_NAME(dsymv)("U", &N, &dOne, invPsi, &N, ones, &iOne, &dZero, tempV, &iOne);
		F77_NAME(dsyr)("U", &N, &tempS, tempV, &iOne, invPsi0, &N);
		
		Memcpy(invPsi1,invPsi0,Nsqr);

		ttPsi0t = quadform(t, invPsi0, N, 1, N) + 1/g;
		tempS = -1/ttPsi0t;
		
		F77_NAME(dsymv)("U", &N, &dOne, invPsi0, &N, t, &iOne, &dZero, tempV, &iOne);
		F77_NAME(dsyr)("U", &N, &tempS, tempV, &iOne, invPsi1, &N);
		
		devs = quadform(y,invPsi1,N,1,N);
		
		logmlike[m] = -0.5*(1.0*N-1)*log(devs) - 0.5*log(tOnePsiOne) - 0.5*log(ttPsi0t) +
					0.5*matrixDet(invPsi, N, N, 1) - 0.5*log(g);
		
	}

	PutRNGstate();
	
	return(logMeanExpLogs(logmlike,iterations));
	
}


SEXP RgibbsTwoSampleAR(SEXP yR, SEXP NR, SEXP tR, SEXP rscaleR, SEXP alphaThetaR, SEXP betaThetaR, SEXP loAreaR, SEXP upAreaR, SEXP iterationsR, SEXP sdmetR, SEXP progressR, SEXP pBar, SEXP rho)
{
	int npars = 7,iterations = INTEGER_VALUE(iterationsR);
	int N = INTEGER_VALUE(NR), progress = INTEGER_VALUE(progressR);
	double rscale = NUMERIC_VALUE(rscaleR);
	double alphaTheta = NUMERIC_VALUE(alphaThetaR);
	double betaTheta = NUMERIC_VALUE(betaThetaR);
	double sdmet = NUMERIC_VALUE(sdmetR);
	double *y = REAL(yR);
	double *t = REAL(tR);
	double loArea = REAL(loAreaR)[0];
	double upArea = REAL(upAreaR)[0];
	
	SEXP chainsR;
	PROTECT(chainsR = allocMatrix(REALSXP, npars, iterations));

	gibbsTwoSampleAR(y, N, t, rscale, alphaTheta, betaTheta, loArea, upArea, iterations, sdmet, REAL(chainsR), progress, pBar, rho);
	
	UNPROTECT(1);
	
	return(chainsR);
}


void gibbsTwoSampleAR(double *y, int N, double *t, double rscale, double alphaTheta, double betaTheta, double loArea, double upArea, int iterations, double sdmet, double *chains, int progress, SEXP pBar, SEXP rho)
{
	int i=0, j=0, m=0,Nsqr=N*N,iOne=1;
	double varMu, meanMu, varDelta, meanDelta, aSig2, bSig2, ag, bg, rscsq=rscale*rscale;
	double mu = 0, sig2=0, delta = 0, theta = 0, g = 1, ldensDelta;
	double tOnePsiOne,psiOne[N],dZero=0,dOne=1;
	double ttPsit,psit[N],loglikeTheta;
	double nullArea;
	
	double tempV[N];
	double ones[N];
	double invPsi[Nsqr];
	AZERO(invPsi,Nsqr);

	for(i=0;i<N;i++)
	{
		ones[i]=1;
		mu += y[i]/(N*1.0);
		sig2 += y[i]*y[i];
	}
	sig2 = (sig2 - N*mu*mu)/(N*1.0-1);
	

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

		// Build invPsi matrix
		invPsi[0] = 1;
		invPsi[N*N-1] = 1;
		invPsi[1] = -theta;
		invPsi[N] = -theta;
		
		for(i=1;i<(N-1);i++)
		{
			invPsi[i + N*i] = (1 + theta*theta);
			invPsi[i + N*(i+1)] = -theta; 
			invPsi[(i+1) + N*i] = -theta; 			
		}
	

		//mu
		tOnePsiOne = quadform(ones, invPsi, N, 1,N);
		F77_NAME(dsymv)("U", &N, &dOne, invPsi, &N, ones, &iOne, &dZero, psiOne, &iOne);

		meanMu = 0;
		varMu = sig2/tOnePsiOne;
		for(i=0;i<N;i++)
		{
			meanMu += (y[i] - delta*t[i])*psiOne[i];
		}
		meanMu = meanMu/tOnePsiOne;
		mu = rnorm(meanMu,sqrt(varMu));
		
		//delta
		ttPsit = quadform(t, invPsi, N, 1,N);
		F77_NAME(dsymv)("U", &N, &dOne, invPsi, &N, t, &iOne, &dZero, psit, &iOne);

		meanDelta = 0;
		varDelta = sig2/(ttPsit + 1/g);
		for(i=0;i<N;i++)
		{
			meanDelta += (y[i] - mu)*psit[i];
		}
		meanDelta = meanDelta/(ttPsit + 1/g);
		delta = rnorm(meanDelta, sqrt(varDelta));

		//deltaDens
		varDelta = 1/(ttPsit + 1/g);
		meanDelta = meanDelta/sqrt(sig2);
		ldensDelta = dnorm(0,meanDelta, sqrt(varDelta), 1);		
		
		//area
		nullArea = pnorm(upArea,meanDelta, sqrt(varDelta), 1, 0) -
		pnorm(loArea,meanDelta, sqrt(varDelta), 1, 0);	
		
		//sig2
		aSig2 = 0.5*(N+1);
		for(i=0;i<N;i++)
		{
			tempV[i] = (y[i] - mu - delta*t[i]);
		}
		bSig2 = 0.5*(quadform(tempV,invPsi,N,1,N) + delta*delta/g);
		sig2 = 1/rgamma(aSig2,1/bSig2);
		
		//g
		ag = 1;
		bg = 0.5*(delta*delta/sig2 + rscsq);
		g = 1/rgamma(ag,1/bg);
		
		//theta
		theta = sampThetaAR(theta, mu, delta, sig2, g, y, N, t, alphaTheta, betaTheta, sdmet);
	
		// write chain
		chains[0*iterations + m] = mu;
		chains[1*iterations + m] = delta;
		chains[2*iterations + m] = ldensDelta;
		chains[3*iterations + m] = sig2;
		chains[4*iterations + m] = g;
		chains[5*iterations + m] = theta;
		chains[6*iterations + m] = nullArea;
		
	
	}

	UNPROTECT(2);
	PutRNGstate();
}


double sampThetaAR(double theta, double mu, double delta, double sig2, double g, double *y, int N, double *t, double alphaTheta, double betaTheta , double sdmet)
{
	// sample theta with Metropolis-Hastings
	double cand, likeRat, b;
	
	cand = theta + rnorm(0,sdmet);
	
	if(cand<0 || cand>1)
	{
		return(theta);
	}
	b = log(runif(0,1));
	likeRat = thetaLogLikeAR(cand, mu, delta, sig2, g, y, N, t, alphaTheta, betaTheta) - thetaLogLikeAR(theta, mu, delta, sig2, g, y, N, t, alphaTheta, betaTheta);
	
	if(b>likeRat){
		return(theta);
	}else{
		return(cand);
	}
}

double thetaLogLikeAR(double theta, double mu, double delta, double sig2, double g, double *y, int N, double *t, double alphaTheta, double betaTheta)
{
	int i;
	double loglike=0,tempV[N],invPsi[N*N],det;
	
	AZERO(invPsi,N*N);
	
	// Build invPsi matrix	
	invPsi[0] = 1;
	invPsi[N*N-1] = 1;
	invPsi[1] = -theta;
	invPsi[N] = -theta;
	for(i=0;i<N;i++)
	{
		tempV[i] = y[i] - mu - delta*t[i];
		
		if(i>0 && i<(N-1)){
			invPsi[i + N*i] = (1 + theta*theta);
			invPsi[i + N*(i+1)] = -theta; 
			invPsi[(i+1) + N*i] = -theta; 			
		}
	}
	det = log(1-theta*theta);

	loglike = 0.5 * det - 0.5/sig2 * quadform(tempV,invPsi,N,1,N) + (alphaTheta-1)*log(theta) + (betaTheta-1)*log(1-theta);
	
	return(loglike);
}




SEXP RgibbsOneSample(SEXP yR, SEXP NR, SEXP rscaleR, SEXP iterationsR, SEXP doIntervalR, SEXP intervalR, SEXP progressR, SEXP pBar, SEXP rho)
{
	
	int npars = 6,iterations = INTEGER_VALUE(iterationsR);
	int N = INTEGER_VALUE(NR), progress = INTEGER_VALUE(progressR);
	double rscale = NUMERIC_VALUE(rscaleR);
	double *y = REAL(yR);
	double *interval = REAL(intervalR);
	int doInterval = INTEGER_VALUE(doIntervalR);
	
	SEXP chainsR;
	PROTECT(chainsR = allocMatrix(REALSXP, npars, iterations));

	gibbsOneSample(y, N, rscale, iterations, REAL(chainsR), doInterval, interval, progress, pBar, rho);
	
	UNPROTECT(1);
	
	return(chainsR);
}

void gibbsOneSample(double *y, int N, double rscale, int iterations, double *chains, int doInterval, double *interval, int progress, SEXP pBar, SEXP rho)
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



SEXP RgibbsOneWayAnova(SEXP yR, SEXP NR, SEXP JR, SEXP IR, SEXP rscaleR, SEXP iterationsR, SEXP progressR, SEXP pBar, SEXP rho)
{
	int iterations = INTEGER_VALUE(iterationsR);
	int *N = INTEGER_POINTER(NR), progress = INTEGER_VALUE(progressR);
	double rscale = NUMERIC_VALUE(rscaleR);
	double *y = REAL(yR);
	int J = INTEGER_VALUE(JR),I = INTEGER_VALUE(IR);
	int j=0,i=0,sumN=0,counter=0,npars=0;
	
	npars = J+5;
	
	for(j=0;j<J;j++){
		sumN += N[j];
	}

	double yVec[sumN];
	
	int whichJ[sumN]; //which j the element belongs to
	
	for(j=0;j<J;j++)
	{
		for(i=0;i<N[j];i++){
			whichJ[counter] = j;
			yVec[counter]=y[j*I + i];
			//Rprintf("i %d, j %d, y %f\n",i,j,yVec[counter]);
			counter++;
		}
	}
	
	//We're going to add another element to returnList for debugging.
	
	SEXP chainsR,returnListR,CMDER,debug;
	PROTECT(chainsR = allocMatrix(REALSXP, npars, iterations));
	PROTECT(returnListR = allocVector(VECSXP,3));
	PROTECT(CMDER = allocVector(REALSXP,4));
	PROTECT(debug = allocVector(VECSXP,2));
	
	gibbsOneWayAnova(yVec, N, J, sumN, whichJ, rscale, iterations, REAL(chainsR), REAL(CMDER), debug, progress, pBar, rho);
	
	SET_VECTOR_ELT(returnListR, 0, chainsR);
    SET_VECTOR_ELT(returnListR, 1, CMDER);
	SET_VECTOR_ELT(returnListR, 2, debug);
	
	UNPROTECT(4);
	
	return(returnListR);
}


double jeffmlikeNWayAov(double *XtCX, double *XtCy, double ytCy, int N, int P, double *g)
{
	double *W,ldetS=0,ldetW,top,bottom1,bottom2,q;
	int i = 0, Psqr=P*P;
	
	W = Memcpy(Calloc(Psqr,double),XtCX,Psqr);
	
	for(i=0;i<P;i++)
	{
		ldetS += log(g[i]);
		W[i + P*i] = W[i + P*i] + 1/g[i];
	}
	InvMatrixUpper(W, P);
	internal_symmetrize(W, P);
	ldetW = matrixDet(W, P, P, 1);
	
	q = quadform(XtCy, W, P, 1, P);
	
	//top = lgamma((N-1)*0.5) + 0.5*ldetW;
	top = 0.5*ldetW;
	bottom1 = 0.5*(N-1) * log(ytCy - q);
	//bottom2 = 0.5*(N-1) * log(M_PI) + 0.5*ldetS;
	bottom2 = 0.5*ldetS;
	
	Free(W);
	
	return(top-bottom1-bottom2);	
}

double jeffSamplerNwayAov(double *samples, int iters, double *XtCX, double *XtCy, double ytCy, int N, int P,int nGs, int *gMap, double *a, double *b, int progress, SEXP pBar, SEXP rho)
{
	int i=0,j=0;
	double avg = 0, g1[nGs],g2[P];
	
	// progress stuff
	SEXP sampCounter, R_fcall;
	int *pSampCounter;
    PROTECT(R_fcall = lang2(pBar, R_NilValue));
	PROTECT(sampCounter = NEW_INTEGER(1));
	pSampCounter = INTEGER_POINTER(sampCounter);

	
	for(i=0;i<iters;i++)
	{
		R_CheckUserInterrupt();
	
		
		//Check progress
		if(progress && !((i+1)%progress)){
			pSampCounter[0]=i+1;
			SETCADR(R_fcall, sampCounter);
			eval(R_fcall, rho); //Update the progress bar
		}


	
		for(j=0;j<nGs;j++)
		{
			g1[j] = 1/rgamma(a[j],1/b[j]);
		}
	
	
		for(j=0;j<P;j++)
		{
			g2[j] = g1[gMap[j]];
		}
	
		samples[i] = jeffmlikeNWayAov(XtCX, XtCy, ytCy, N, P, g2);
		if(i==0)
		{
			avg = samples[i];
		}else{
			avg = LogOnePlusExpX(samples[i]-avg)+avg;
			//avg = LogOnePlusX(exp(avg-samples[i]))+samples[i];
		}
	}
	
	UNPROTECT(2);
	
	return(avg-log(iters));
	
}

SEXP RLogMeanExpLogs(SEXP Rv, SEXP Rlen)
{
	double *v;
	int len;
	SEXP ret;
	double *retp;
	
	len = INTEGER_VALUE(Rlen);
	v = REAL(Rv);
	PROTECT(ret = allocVector(REALSXP,1));
	retp = REAL(ret);
	
	retp[0] = logMeanExpLogs(v,len);
	
	UNPROTECT(1);
	return(ret);
}

SEXP RLogCumMeanExpLogs(SEXP Rv, SEXP Rlen)
{
	double *v;
	int len;
	SEXP ret;
	double *retp;
	
	len = INTEGER_VALUE(Rlen);
	v = REAL(Rv);
	PROTECT(ret = allocVector(REALSXP,len));
	retp = REAL(ret);
	
	logCumMeanExpLogs(v, len, retp);
	
	UNPROTECT(1);
	return(ret);
}



double logSumExpLogs(double *v, int len)
{
	int i=0;
	double sum=v[0];
	
	for(i=1;i<len;i++)
	{
		sum = LogOnePlusExpX(v[i]-sum)+sum;
	}
	
	return(sum);
}

double logMeanExpLogs(double *v, int len)
{	
	double sum=0;
	sum = logSumExpLogs(v,len);
	return(sum - log(len));
}

void logCumMeanExpLogs(double *v, int len, double *ret)
{	
	int i=0;
	double sums[len];
	
	ret[0]=v[0];
	sums[0]=v[0];
	
	for(i=1;i<len;i++)
	{
		sums[i] = LogOnePlusExpX(v[i]-sums[i-1])+sums[i-1];
		ret[i] = sums[i] - log(i+1); 
	}
}

SEXP RjeffSamplerNwayAov(SEXP Riters, SEXP RXtCX, SEXP RXtCy, SEXP RytCy, SEXP RN, SEXP RP, SEXP RnGs, SEXP RgMap, SEXP Ra, SEXP Rb, SEXP progressR, SEXP pBar, SEXP rho)
{
	double *a,*b,*XtCy,*XtCX,ytCy,*samples,*avg;
	int iters,nGs,*gMap,N,P,progress;

	SEXP Rsamples,Ravg,returnList;
	
	iters = INTEGER_VALUE(Riters);
	XtCX = REAL(RXtCX);
	XtCy = REAL(RXtCy);
	ytCy = REAL(RytCy)[0];
	nGs = INTEGER_VALUE(RnGs);
	gMap = INTEGER_POINTER(RgMap);
	a = REAL(Ra);
	b = REAL(Rb);
	N = INTEGER_VALUE(RN);
	P = INTEGER_VALUE(RP);
	progress = INTEGER_VALUE(progressR);
	
	PROTECT(Rsamples = allocVector(REALSXP,iters));
	PROTECT(Ravg = allocVector(REALSXP,1));
	PROTECT(returnList = allocVector(VECSXP,2));
	
	samples = REAL(Rsamples);
	avg = REAL(Ravg);
	
	GetRNGstate();
	avg[0] = jeffSamplerNwayAov(samples, iters, XtCX, XtCy, ytCy, N, P, nGs, gMap, a, b, progress, pBar, rho);
	PutRNGstate();
	
	SET_VECTOR_ELT(returnList, 0, Ravg);
    SET_VECTOR_ELT(returnList, 1, Rsamples);
	
	UNPROTECT(3);
	
	return(returnList);
}


void gibbsOneWayAnova(double *y, int *N, int J, int sumN, int *whichJ, double rscale, int iterations, double *chains, double *CMDE, SEXP debug, int progress, SEXP pBar, SEXP rho)
{
	int i=0,j=0,m=0,Jp1sq = (J+1)*(J+1),Jsq=J*J,Jp1=J+1,npars=0;
	double ySum[J],yBar[J],sumy2[J],densDelta=0;
	double sig2=1,g=1;
	double XtX[Jp1sq], ZtZ[Jsq];
	double Btemp[Jp1sq],B2temp[Jsq],tempBetaSq=0;
	double muTemp[J],oneOverSig2temp=0;
	double beta[J+1],grandSum=0,grandSumSq=0;
	double shapeSig2 = (sumN+J*1.0)/2, shapeg = (J+1.0)/2;
	double scaleSig2=0, scaleg=0;
	double Xty[J+1],Zty[J];
	double logDet=0;
	double rscaleSq=rscale*rscale;
	
	double logSumSingle=0,logSumDouble=0;

	// for Kahan sum
	double kahanSumSingle=0, kahanSumDouble=0;
	double kahanCSingle=0,kahanCDouble=0;
	double kahanTempT=0, kahanTempY=0;
	
	int iOne=1;
	double dZero=0;
	
	/* debug stuff 
	SEXP debugSigmas, debugMus;
	PROTECT(debugMus = allocMatrix(REALSXP, J+1, iterations));
	PROTECT(debugSigmas = alloc3Darray(REALSXP,J+1,J+1,iterations));
	double *cDebugSigmas = REAL(debugSigmas);
	double *cDebugMus = REAL(debugMus);
	
	
	SEXP R_mvtnormCall,SigmaR,MuR,mvtArgs,mvtnormSamp;
	PROTECT(R_mvtnormCall = lang2(mvtnorm, R_NilValue));
	PROTECT(SigmaR = allocMatrix(REALSXP, J+1, J+1));
	PROTECT(MuR = allocMatrix(REALSXP, J+1, 1));
	PROTECT(mvtArgs = allocVector(VECSXP, 2));
	PROTECT(mvtnormSamp = allocVector(VECSXP,1));
	 end debug stuff */
	
	
	
	// progress stuff
	SEXP sampCounter, R_fcall;
	int *pSampCounter;
    PROTECT(R_fcall = lang2(pBar, R_NilValue));
	PROTECT(sampCounter = NEW_INTEGER(1));
	pSampCounter = INTEGER_POINTER(sampCounter);
	
	npars=J+5;
	
	GetRNGstate();

	// Initialize to 0
	AZERO(XtX,Jp1sq);
	AZERO(ZtZ,Jsq);
	AZERO(beta,Jp1);
	AZERO(ySum,J);
	AZERO(sumy2,J);
	
	// Create vectors
	for(i=0;i<sumN;i++)
	{
		j = whichJ[i];
		ySum[j] += y[i];
		//Rprintf("%d %d %f\n",j,whichJ[i],y[i]);
		sumy2[j] += y[i]*y[i];
		grandSum += y[i];
		grandSumSq += y[i]*y[i];
	}
	
	
	// create design matrices
	XtX[0]=sumN;	
	for(j=0;j<J;j++)
	{
		XtX[j+1]=N[j];
		XtX[(J+1)*(j+1)]=N[j];
		XtX[(j+1)*(J+1) + (j+1)] = N[j];
		ZtZ[j*J + j] = N[j];
		yBar[j] = ySum[j]/(1.0*N[j]);
	}
	
	Xty[0] = grandSum;	
	Memcpy(Xty+1,ySum,J);
	Memcpy(Zty,ySum,J);
	
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
		F77_CALL(dsymv)("U", &Jp1, &oneOverSig2temp, Btemp, &Jp1, Xty, &iOne, &dZero, beta, &iOne);
		
		//AZERO(beta,Jp1);
		/* for debugging 
		Memcpy(cDebugMus + m*Jp1,beta,Jp1);
		Memcpy(cDebugSigmas + m*Jp1sq,Btemp,Jp1sq);
		
	
		Memcpy(REAL(SigmaR),Btemp,Jp1sq);
		Memcpy(REAL(MuR),beta,Jp1);
		
		SET_VECTOR_ELT(mvtArgs, 0, MuR);
		SET_VECTOR_ELT(mvtArgs, 1, SigmaR);

		SETCADR(R_mvtnormCall, mvtArgs);
		SET_VECTOR_ELT(mvtnormSamp, 0, eval(R_mvtnormCall, rho)); // Get sample
		Memcpy(beta,REAL(VECTOR_ELT(mvtnormSamp,0)),Jp1);
		 end debugging */
		rmvGaussianC(beta, Btemp, J+1);
		Memcpy(&chains[npars*m],beta,J+1);	
		
		
		// calculate density (Single Standardized)
		
		Memcpy(B2temp,ZtZ,Jsq);
		densDelta = -J*0.5*log(2*M_PI);
		for(j=0;j<J;j++)
		{
			B2temp[j*J+j] += 1/g;
			muTemp[j] = (ySum[j]-N[j]*beta[0])/sqrt(sig2);
		}
		InvMatrixUpper(B2temp, J);
		internal_symmetrize(B2temp,J);
		logDet = matrixDet(B2temp,J,J,1);
		densDelta += -0.5*quadform(muTemp, B2temp, J, 1, J);
		densDelta += -0.5*logDet;
		if(m==0){
			logSumSingle = densDelta;
			kahanSumSingle = exp(densDelta);
		}else{
			logSumSingle =  logSumSingle + LogOnePlusX(exp(densDelta-logSumSingle));
			kahanTempY = exp(densDelta) - kahanCSingle;
			kahanTempT = kahanSumSingle + kahanTempY;
			kahanCSingle = (kahanTempT - kahanSumSingle) - kahanTempY;
			kahanSumSingle = kahanTempT;
		}
		chains[npars*m + (J+1) + 0] = densDelta;
		
		
		// calculate density (Double Standardized)
		densDelta += 0.5*J*log(g);
		if(m==0){
			logSumDouble = densDelta;
			kahanSumDouble = exp(densDelta);
		}else{
			logSumDouble =  logSumDouble + LogOnePlusX(exp(densDelta-logSumDouble));
			kahanTempY = exp(densDelta) - kahanCDouble;
			kahanTempT = kahanSumDouble + kahanTempY;
			kahanCDouble = (kahanTempT - kahanSumDouble) - kahanTempY;
			kahanSumDouble = kahanTempT;
		}
		chains[npars*m + (J+1) + 1] = densDelta;
		
		
		
		// sample sig2
		tempBetaSq = 0;
		scaleSig2 = grandSumSq - 2*beta[0]*grandSum + beta[0]*beta[0]*sumN;
		for(j=0;j<J;j++)
		{
			scaleSig2 += -2.0*(yBar[j]-beta[0])*N[j]*beta[j+1] + (N[j]+1/g)*beta[j+1]*beta[j+1];
			tempBetaSq += beta[j+1]*beta[j+1];
		}
		scaleSig2 *= 0.5;
		sig2 = 1/rgamma(shapeSig2,1/scaleSig2);
		chains[npars*m + (J+1) + 2] = sig2;
	
		// sample g
		scaleg = 0.5*(tempBetaSq/sig2 + rscaleSq);
		g = 1/rgamma(shapeg,1/scaleg);
		chains[npars*m + (J+1) + 3] = g;

	}
	
	CMDE[0] = logSumSingle - log(iterations);
	CMDE[1] = logSumDouble - log(iterations);
	CMDE[2] = log(kahanSumSingle) - log(iterations);
	CMDE[3] = log(kahanSumDouble) - log(iterations);
	
	/* for debugging 
	SET_VECTOR_ELT(debug, 0, debugMus);
    SET_VECTOR_ELT(debug, 1, debugSigmas);
	 end debugging */
	
	UNPROTECT(2);
	PutRNGstate();
	
}


/*
double LogOnePlusX(double x)
{
    if (x <= -1.0)
    {
		error("Error: attempt to compute log(1+X) where X is %f\n.",x);
	}

    if (fabs(x) > 1/(10000.0))
    {
        // x is large enough that the obvious evaluation is OK
        return( log(1.0 + x) );
    }

    // Use Taylor approx. log(1 + x) = x - x^2/2 with error roughly x^3/3
    // Since |x| < 10^-4, |x|^3 < 10^-12, relative error less than 10^-8

    return( (-0.5*x + 1.0)*x );
}
*/

// Calculate log(1 + x), preventing loss of precision for small values of x.
// The input x must be larger than -1 so that log(1 + x) is real.
// Taken from http://www.codeproject.com/KB/recipes/avoiding_overflow.aspx
double LogOnePlusX(double x)
{
    if (x <= -1.0)
    {
        
        error("Attempt to compute log(1+x) on x<=-1!");
    }

	if (fabs(x) > 0.375)
    {
        // x is sufficiently large that the obvious evaluation is OK
        return log(1.0 + x);
    }

	// For smaller arguments we use a rational approximation
	// to the function log(1+x) to avoid the loss of precision
	// that would occur if we simply added 1 to x then took the log.

    const double p1 =  -0.129418923021993e+01;
    const double p2 =   0.405303492862024e+00;
    const double p3 =  -0.178874546012214e-01;
    const double q1 =  -0.162752256355323e+01;
    const double q2 =   0.747811014037616e+00;
    const double q3 =  -0.845104217945565e-01;
    double t, t2, w;

    t = x/(x + 2.0);
    t2 = t*t;
    w = (((p3*t2 + p2)*t2 + p1)*t2 + 1.0)/(((q3*t2 + q2)*t2 + q1)*t2 + 1.0);
    return 2.0*t*w;
}

// return log(1 + exp(x)), preventing cancellation and overflow */
// From http://www.codeproject.com/KB/recipes/avoiding_overflow.aspx
double LogOnePlusExpX(double x)
{
    const double LOG_DBL_EPSILON = log(DBL_EPSILON);
    const double LOG_ONE_QUARTER = log(0.25);

    if (x > -LOG_DBL_EPSILON)
    {
        // log(exp(x) + 1) == x to machine precision
        return x;
    }
    else if (x > LOG_ONE_QUARTER)
    {
        return log( 1.0 + exp(x) );
    }
    else
    {
        // Prevent loss of precision that would result from adding small argument to 1.
        return LogOnePlusX( exp(x) );
    }
}


// Compute determinant of an N by N matrix A
double matrixDet(double *A, int N, int LDA, int doLog)
{
//SUBROUTINE DPOTRF( UPLO, N, A, LDA, INFO )
	int i=0, info=0;
	double B[N*N], logDet=0;
	
	//Memcpy(B,A,N*N);
	for(i=0;i<N;i++){
		Memcpy(&B[i*N],&A[i*LDA],N);
	}
	
	F77_CALL(dpotrf)("U", &N, B, &N, &info);
	if(info){
		Rprintf("Cholesky decomposition in matrixDet() returned nonzero info %d.\n",info);
	}
	
	for(i=0;i<N;i++)
	{
		logDet += 2 * log(B[i*N+i]);
	}
	
	if(doLog){
		return(logDet);
	}else{
		return(exp(logDet));
	}
}

double quadform(double *x, double *A, int N, int incx, int LDA)
{
  
  int Nsqr = N*N,info,i=0,j=0;
  double *B = Calloc(Nsqr,double);
  //double one=1;
  //double zero=0;
  double sumSq=0;
  double y[N];
  int iOne=1;

  for(i=0;i<N;i++){
    y[i] = x[i*incx];
  }
  for(i=0;i<N;i++){
	Memcpy(&B[i*N],&A[i*LDA],N);
  }

  F77_NAME(dpotrf)("U", &N, B, &N, &info);
  F77_NAME(dtrmv)("U","N","N", &N, B, &N, y, &iOne);
  
  for(i=0;i<N;i++){
    sumSq += y[i]*y[i];
  }
  
  Free(B);
  
  return(sumSq);
}

//version without cholesky decomp for non positive definite matrix A
double quadform2(double *x, double *A, int N, int incx, int LDA)
{
  
  int i=0;
  double dOne=1;
  double dZero=0;
  double sumSq=0;
  double y[N];
  int iOne=1;
  
  F77_NAME(dgemv)("N", &N, &N, &dOne, A, &LDA, x, &incx, &dZero, y, &iOne);
 
  for(i=0;i<N;i++){
    sumSq += y[i]*x[i];
  }
    
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

  F77_NAME(dpotrf)("L", &p, scCp, &p, &info);
  if (info){
	error("Nonzero info from dpotrf: Sigma matrix is not positive-definite");
  }
  //GetRNGstate();
  for(j=0;j<p;j++)
    {
      ans[j] = rnorm(0,1);
    }
  F77_NAME(dtrmv)("L","N","N", &p, scCp, &p, ans, &intOne);
  F77_NAME(daxpy)(&p, &one, ans, &intOne, mu, &intOne);
  //PutRNGstate();
  Free(scCp);
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



