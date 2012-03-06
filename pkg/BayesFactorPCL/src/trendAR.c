#include "BFPCL.h"
#include "trendAR.h"


SEXP RgibbsTwoSampleAR_trend(SEXP yR, SEXP NR, SEXP XR, SEXP pR, SEXP rscaleIntR, SEXP rscaleSlpR, SEXP alphaThetaR, SEXP betaThetaR, SEXP loIntR, SEXP upIntR, SEXP loSlpR, SEXP upSlpR, SEXP iterationsR, SEXP sdmetR, SEXP progressR, SEXP pBar, SEXP rho)
{
	int iterations = INTEGER_VALUE(iterationsR);
	int N = INTEGER_VALUE(NR), progress = INTEGER_VALUE(progressR);
	double rscaleInt = NUMERIC_VALUE(rscaleIntR);
	double rscaleSlp = NUMERIC_VALUE(rscaleSlpR);
	double alphaTheta = NUMERIC_VALUE(alphaThetaR);
	double betaTheta = NUMERIC_VALUE(betaThetaR);
	double sdmet = NUMERIC_VALUE(sdmetR);
	double *y = REAL(yR);
	double *X = REAL(XR);
	int p = INTEGER_VALUE(pR);
	int npars = p + 4;
	
	double loInt = REAL(loIntR)[0];
	double upInt = REAL(upIntR)[0];
	double loSlp = REAL(loSlpR)[0];
	double upSlp = REAL(upSlpR)[0];
	
	
	SEXP chainsR;
	PROTECT(chainsR = allocMatrix(REALSXP, npars, iterations));

	// means
	SEXP postdensR;
	PROTECT(postdensR = allocMatrix(REALSXP, 1, 5));
	SEXP returnList;
	PROTECT(returnList = allocVector(VECSXP, 2));
	

	gibbsTwoSampleAR_trend(y, N, X, p, rscaleInt, rscaleSlp, alphaTheta, betaTheta, loInt, upInt, loSlp, upSlp, iterations, sdmet, REAL(chainsR), REAL(postdensR), progress, pBar, rho);
	
	SET_VECTOR_ELT(returnList, 0, chainsR);
    SET_VECTOR_ELT(returnList, 1, postdensR);
	
	UNPROTECT(3);
	
	return(returnList);
}


void gibbsTwoSampleAR_trend(double *y, int N, double *X, int p, double rscaleInt, double rscaleSlp, double alphaTheta, double betaTheta, double loInt, double upInt, double loSlp, double upSlp, int iterations, double sdmet, double *chains, double *postdens, int progress, SEXP pBar, SEXP rho)
{
	int i=0, j=0, m=0,Nsqr=N*N,iOne=1,iTwo=2,iThree=3;
	double aSig2, bSig2, ag, bg, rscIntsq=rscaleInt*rscaleInt, rscSlpsq=rscaleSlp*rscaleSlp;
	double sig2=0, theta = 0, g1 = 1, g2=1;
	double dZero=0,dOne=1,dNegOne=-1;
	double dOneOverSig2;
	
	double ldensFullRestrict,ldensSlpRestrict,ldensIntRestrict,lAreaSlp,lAreaInt;
	
	double tempV[N];
	double tempV2[N];
	double invPsi[Nsqr];
	AZERO(invPsi,Nsqr);
	double beta[p];
	AZERO(beta,p);
	double tempNbyP[N*p];
	double Sigma[p*p];
	double tempNby2[2*N];
	double temp2by2[2*2];
	double temp2by2_2[2*2];
	double temp2by1[2];
	double temp2by1_2[2];

	
	double Xg[2*N], X2[2*N], X1[N], X3[3*N], beta2[2], beta3[3], betag[2];
	double betaVar;
	double betaMean;

	
	Memcpy(Xg, X, N);
	Memcpy(Xg + N, X + 2*N, N);
	Memcpy(X2, X + N, N);
	Memcpy(X2 + N, X + 3*N, N);
	
	int npars = p + 4;
	
	for(i=0;i<N;i++)
	{
		beta[0] += y[i]/(N*1.0);
		sig2 += y[i]*y[i];
	}
	sig2 = (sig2 - N*beta[0]*beta[0])/(N*1.0-1);
	

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
						
		dOneOverSig2 = 1/sig2;
		F77_NAME(dgemm)("N", "N", &N, &p, &N, &dOneOverSig2, invPsi, &N, X, &N, &dZero, tempNbyP, &N);
				
		F77_NAME(dgemm)("T", "N", &p, &p, &N, &dOne, X, &N, tempNbyP, &N, &dZero, Sigma, &p);

		Sigma[1 + p] += 1/(sig2 * g1);
		Sigma[3 + p*3] += 1/(sig2 * g2);
		
		InvMatrixUpper(Sigma, p);
		internal_symmetrize(Sigma, p);			
		
		F77_NAME(dgemv)("T", &N, &p, &dOne, tempNbyP, &N, y, &iOne, &dZero, tempV, &iOne);
		F77_NAME(dgemv)("N", &p, &p, &dOne, Sigma, &p, tempV, &iOne, &dZero, tempV2, &iOne);
		
		rmvGaussianC(tempV2, Sigma, p);
		Memcpy(beta, tempV2, p);
	
		//densities
		
		//slope restricted
		Memcpy(X1, X + 3*N, N);
		Memcpy(X3, Xg, 2*N);
		Memcpy(X3 + 2*N, X + N, N);
		beta3[0] = beta[0];
		beta3[1] = beta[2];
		beta3[2] = beta[1];
		Memcpy(tempV,y,N);
		F77_NAME(dgemv)("N", &N, &iThree, &dNegOne, X3, &N, beta3, &iOne, &dOne, tempV, &iOne);
		betaVar = 0;
		betaMean = 0;
		betaVar = quadform(X1,invPsi,N,1,N);
		F77_NAME(dgemv)("N", &N, &N, &dOne, invPsi, &N, X1, &iOne, &dZero, tempV2, &iOne);
		
		for(i=0;i<N;i++)
			betaMean += tempV2[i] * tempV[i];

		betaVar = 1 / (betaVar + 1/g2);
		betaMean = betaVar * betaMean;

		ldensSlpRestrict = dnorm(0,betaMean/sqrt(sig2),sqrt(betaVar),1);		
		lAreaSlp = pnorm(upSlp,betaMean/sqrt(sig2),sqrt(betaVar),1,0) - pnorm(loSlp,betaMean/sqrt(sig2),sqrt(betaVar),1,0);
		
		if(m==0)
		{
			postdens[1] = ldensSlpRestrict;
			postdens[4] = lAreaSlp;
		}else{
			postdens[1] = LogOnePlusExpX(ldensSlpRestrict-postdens[1])+postdens[1];
			postdens[4] += lAreaSlp;
		}
		
		
		
		//intercept restricted
		Memcpy(X1, X + N, N);
		Memcpy(X3 + 2*N, X + 3*N, N);
		beta3[0] = beta[0];
		beta3[1] = beta[2];
		beta3[2] = beta[3];
		Memcpy(tempV,y,N);
		F77_NAME(dgemv)("N", &N, &iThree, &dNegOne, X3, &N, beta3, &iOne, &dOne, tempV, &iOne);
		betaVar = 0;
		betaMean = 0;
		betaVar = quadform(X1,invPsi,N,1,N);
		F77_NAME(dgemv)("N", &N, &N, &dOne, invPsi, &N, X1, &iOne, &dZero, tempV2, &iOne);
		
		for(i=0;i<N;i++)
			betaMean += tempV2[i] * tempV[i];

		betaVar = 1 / (betaVar + 1/g1);
		betaMean = betaVar * betaMean;
		ldensIntRestrict = dnorm(0,betaMean/sqrt(sig2),sqrt(betaVar),1);	
		lAreaInt = pnorm(upInt,betaMean/sqrt(sig2),sqrt(betaVar),1,0) - pnorm(loInt,betaMean/sqrt(sig2),sqrt(betaVar),1,0);
		
		if(m==0)
		{
			postdens[2] = ldensIntRestrict;
			postdens[3] = lAreaInt;
		}else{
			postdens[2] = LogOnePlusExpX(ldensIntRestrict-postdens[2])+postdens[2];
			postdens[3] += lAreaInt;
		}
		
		
		//Both restricted
		betag[0] = beta[0];
		betag[1] = beta[2];
		Memcpy(tempV,y,N);
		
		//(y - Xg%*%betag)
		F77_NAME(dgemv)("N", &N, &iTwo, &dNegOne, Xg, &N, betag, &iOne, &dOne, tempV, &iOne);
				
		// invPsi%*%X2		
		F77_NAME(dgemm)("N", "N", &N, &iTwo, &N, &dOne, invPsi, &N, X2, &N, &dZero, tempNby2, &N);	
		//t(X2)%*%invPsi%*%(y - Xg%*%betag)
		F77_NAME(dgemv)("T", &N, &iTwo, &dOne, tempNby2, &N, tempV, &iOne, &dZero, temp2by1, &iOne);
		
		// t(X2)%*%invPsi%*%X2		
		F77_NAME(dgemm)("T", "N", &iTwo, &iTwo, &N, &dOne, X2, &N, tempNby2, &N, &dZero, temp2by2, &iTwo);
		
		temp2by2[0] += 1/g1;
		temp2by2[3] += 1/g2;
		
		Memcpy(temp2by2_2,temp2by2,4);
		
		InvMatrixUpper(temp2by2, 2);
		internal_symmetrize(temp2by2, 2);	
		
		dOneOverSig2 = 1/sqrt(sig2);
		
		//Sigma%*%t(X2)%*%invPsi%*%X2
		F77_NAME(dgemv)("N", &iTwo, &iTwo, &dOneOverSig2, temp2by2, &iTwo, temp2by1, &iOne, &dZero, temp2by1_2, &iOne);
		
		ldensFullRestrict = -log(2 * M_PI) - 0.5*matrixDet(temp2by2,2,2,1) - 0.5*quadform(temp2by1_2,temp2by2_2,2,1,2);
		
		if(m==0)
		{
			postdens[0] = ldensFullRestrict;
		}else{
			postdens[0] = LogOnePlusExpX(ldensFullRestrict-postdens[0])+postdens[0];
		}
		
		//sig2
		Memcpy(tempV,y,N);
		F77_NAME(dgemv)("N", &N, &p, &dNegOne, X, &N, beta, &iOne, &dOne, tempV, &iOne);
		aSig2 = 0.5*(N+2);
		bSig2 = 0.5*(quadform(tempV,invPsi,N,1,N) + beta[1]*beta[1]/g1 + beta[3]*beta[3]/g2);

		sig2 = 1/rgamma(aSig2,1/bSig2);
		
		//g1
		ag = 1;
		bg = 0.5*(beta[1]*beta[1]/sig2 + rscIntsq);
		g1 = 1/rgamma(ag,1/bg);

		//g2
		ag = 1;
		bg = 0.5*(beta[3]*beta[3]/sig2 + rscSlpsq);
		g2 = 1/rgamma(ag,1/bg);
		
		//theta
		theta = sampThetaAR_trend(theta, beta, X, sig2, y, N, p, alphaTheta, betaTheta, sdmet);
	
		// write chain
		Memcpy(chains + m*npars, beta, p);
		chains[p + m*npars] = sig2;
		chains[1 + p + m*npars] = g1;
		chains[2 + p + m*npars] = g2;
		chains[3 + p + m*npars] = theta;
		//chains[4 + p + m*npars] = ldensFullRestrict;
		//chains[5 + p + m*npars] = ldensSlpRestrict;
		//chains[6 + p + m*npars] = ldensIntRestrict;
		//chains[7 + p + m*npars] = lAreaInt;
		//chains[8 + p + m*npars] = lAreaSlp;
	
	}

	UNPROTECT(2);
	PutRNGstate();
}


double sampThetaAR_trend(double theta, double *beta, double *X, double sig2, double *y, int N, int p, double alphaTheta, double betaTheta, double sdmet)
{
	// sample theta with Metropolis-Hastings
	double cand, likeRat, b;
	
	cand = theta + rnorm(0,sdmet);
	
	if(cand<0 || cand>1)
	{
		return(theta);
	}
	b = log(runif(0,1));
	likeRat = thetaLogLikeAR_trend(cand, beta, X, sig2, y, N, p, alphaTheta, betaTheta) - thetaLogLikeAR_trend(theta, beta, X, sig2, y, N, p, alphaTheta, betaTheta);
	
	if(b>likeRat){
		return(theta);
	}else{
		return(cand);
	}
}

double thetaLogLikeAR_trend(double theta, double *beta, double *X, double sig2, double *y, int N, int p, double alphaTheta, double betaTheta)
{
	int i,iOne=1;
	double loglike=0,tempV[N],invPsi[N*N],det,dNegOne=-1,dOne=1;
	
	Memcpy(tempV,y,N);
	F77_NAME(dgemv)("N", &N, &p, &dNegOne, X, &N, beta, &iOne, &dOne, tempV, &iOne);

	AZERO(invPsi,N*N);
	
	// Build invPsi matrix	
	invPsi[0] = 1;
	invPsi[N*N-1] = 1;
	invPsi[1] = -theta;
	invPsi[N] = -theta;
	for(i=0;i<N;i++)
	{
		
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
