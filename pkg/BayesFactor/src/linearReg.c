#include "BFPCL.h"


SEXP RGibbsLinearReg(SEXP Riters, SEXP RCny, SEXP RX, SEXP RXtX, SEXP RXtCnX, SEXP RXtCny, SEXP RN, SEXP RP, SEXP Rr, SEXP Rsig2start, SEXP progressR, SEXP pBar, SEXP callback, SEXP rho)
{
	double *XtX,*XtCnX,*XtCny,*Cny,*samples,*X, r, sig2start;
	int iters,N,P,progress;

	SEXP Rsamples;
	
	iters = INTEGER_VALUE(Riters);
	X = REAL(RX);
	Cny = REAL(RCny);
	XtX = REAL(RXtX);
	XtCnX = REAL(RXtCnX);
	XtCny = REAL(RXtCny);
	r = REAL(Rr)[0];
  sig2start = REAL(Rsig2start)[0];
  N = INTEGER_VALUE(RN);
	P = INTEGER_VALUE(RP);
	progress = INTEGER_VALUE(progressR);
	
	PROTECT(Rsamples = allocMatrix(REALSXP, 2 + P, iters));

	samples = REAL(Rsamples);
	
	GetRNGstate();
	GibbsLinearReg(samples, iters, Cny, X, XtX, XtCnX, XtCny, N, P, r, sig2start, progress, pBar, callback, rho);
	PutRNGstate();
	
	
	UNPROTECT(1);
	
	return(Rsamples);
}

void GibbsLinearReg(double *chains, int iters, double *Cny, double *X, double *XtX, double *XtCnX, double *XtCny, int N, int P, double r, double sig2start, int progress, SEXP pBar, SEXP callback, SEXP rho)
{
	int i=0,j=0,k=0, nPars=P+2, PSq=P*P, iOne=1;
	double g=1, Sigma[PSq], SSq, oneOverSig2,dZero=0,dOne=1,dnegOne=-1;
	
	double beta[P], sig2=sig2start, yTemp[N],bTemp, XtXoN[PSq];
	
	// progress stuff
	SEXP sampCounter, R_fcall;
	int *pSampCounter;
    PROTECT(R_fcall = lang2(pBar, R_NilValue));
	PROTECT(sampCounter = NEW_INTEGER(1));
	pSampCounter = INTEGER_POINTER(sampCounter);

	Memcpy(XtXoN,XtX,PSq);
	for(i=0;i<PSq;i++){
		XtXoN[i] = XtXoN[i] / (N*1.0); 
	}
	
	for(i=0;i<iters;i++)
	{
		R_CheckUserInterrupt();
			
		//Check progress
		if(progress && !((i+1)%progress)){
			pSampCounter[0]=i+1;
			SETCADR(R_fcall, sampCounter);
			eval(R_fcall, rho); //Update the progress bar
		}
		
		// Sample beta
		
    Memcpy(Sigma,XtX,PSq);
		for(j=0;j<PSq;j++){
			Sigma[j] = ( XtCnX[j] + XtXoN[j] / g ) / sig2; 
		}
	
		InvMatrixUpper(Sigma,P);
		internal_symmetrize(Sigma,P);
		oneOverSig2 = 1/sig2;
		F77_CALL(dsymv)("U", &P, &oneOverSig2, Sigma, &P, XtCny, &iOne, &dZero, beta, &iOne);	
		rmvGaussianC(beta, Sigma, P);					
		
    
		// Sample Sig2
		SSq = 0;
		Memcpy(yTemp,Cny,N);
		
		F77_NAME(dgemv)("N", &N, &P, &dOne, X, &N, beta, &iOne, &dnegOne, yTemp, &iOne);
		for(j=0;j<N;j++)
		{
			SSq += yTemp[j] * yTemp[j];
		}
		SSq += quadform(beta,XtX,P,1,P) / (g*N);
    sig2 = 1/rgamma( N/2.0, 2.0/SSq );
		
		
		// Sample g
		SSq = 0;
		SSq += quadform(beta,XtX,P,1,P) / (sig2*N) + r*r;
    g  = 1/rgamma( 1, 2.0/SSq);
		
	
		Memcpy(chains + nPars*i,beta,P);	
		chains[nPars*i + P] = sig2;	
		chains[nPars*i + P + 1] = g;	
	}
	
	UNPROTECT(2);
		
}


