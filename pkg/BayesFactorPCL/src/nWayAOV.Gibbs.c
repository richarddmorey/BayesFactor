#include "BFPCL.h"


SEXP RGibbsNwayAov(SEXP Riters, SEXP Ry, SEXP RX, SEXP RXtX, SEXP RXty, SEXP RN, SEXP RP, SEXP RnGs, SEXP RgMap, SEXP Rr, SEXP progressR, SEXP pBar, SEXP rho)
{
	double *a,*b,*Xty,*XtX,*samples,*X, *y, *r;
	int iters,nGs,*gMap,N,P,progress;

	SEXP Rsamples;
	
	iters = INTEGER_VALUE(Riters);
	X = REAL(RX);
	y = REAL(Ry);
	XtX = REAL(RXtX);
	Xty = REAL(RXty);
	nGs = INTEGER_VALUE(RnGs);
	gMap = INTEGER_POINTER(RgMap);
	r = REAL(Rr);
	N = INTEGER_VALUE(RN);
	P = INTEGER_VALUE(RP);
	progress = INTEGER_VALUE(progressR);
	
	PROTECT(Rsamples = allocMatrix(REALSXP,iters, 2 + P + nGs));

	samples = REAL(Rsamples);
	
	GetRNGstate();
	GibbsNwayAov(samples, iters, y, X, XtX, Xty, N, P, nGs, gMap, r, progress, pBar, rho);
	PutRNGstate();
	
	
	UNPROTECT(1);
	
	return(Rsamples);
}

void GibbsNwayAov(double *chains, int iters, double *y, double *X, double *XtX, double *Xty, int N, int P, int nGs, int *gMap, double *r, int progress, SEXP pBar, SEXP rho)
{
	int i=0,j=0,k=0, nPars=2+P+nGs, P1Sq=(P+1)*(P+1), P1=P+1, iOne=1;
	double *g1, Sigma[P1Sq], SSq, oneOverSig2,dZero=0,dOne=1,dnegOne=-1;
	
	double g[nGs], beta[P+1], sig2, yTemp[N], SSqG[nGs],bTemp,nParG[nGs];
	
	// progress stuff
	SEXP sampCounter, R_fcall;
	int *pSampCounter;
    PROTECT(R_fcall = lang2(pBar, R_NilValue));
	PROTECT(sampCounter = NEW_INTEGER(1));
	pSampCounter = INTEGER_POINTER(sampCounter);

	AZERO(nParG,nGs);
	sig2 = 1;
	for(i=0;i<P;i++)
	{
		nParG[gMap[i]]++;
	}
	for(i=0;i<nGs;i++)
	{
		g[i] = 1;
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
		Memcpy(Sigma,XtX,P1Sq);
		for(j=0;j<(P+1);j++)
		{
			if(j>0){
				Sigma[ j + (P+1)*j] = (Sigma[ j + (P+1)*j] + 1/g[gMap[j-1]])/sig2;
			}else{
				Sigma[ j + (P+1)*j] = Sigma[ j + (P+1)*j]/sig2;
			}
			for(k=j+1;k<(P+1);k++)
			{
				Sigma[j + (P+1)*k] = Sigma[j + (P+1)*k]/sig2;
			}
		}
		InvMatrixUpper(Sigma,P+1);
		internal_symmetrize(Sigma,P+1);
		oneOverSig2 = 1/sig2;
		F77_CALL(dsymv)("U", &P1, &oneOverSig2, Sigma, &P1, Xty, &iOne, &dZero, beta, &iOne);	
		rmvGaussianC(beta, Sigma, P+1);
				
		// Sample Sig2
		SSq = 0;
		Memcpy(yTemp,y,N);
		
		F77_NAME(dgemv)("N", &N, &P1, &dOne, X, &N, beta, &iOne, &dnegOne, yTemp, &iOne);
		for(j=0;j<N;j++)
		{
			SSq += yTemp[j] * yTemp[j];
		}
		for(j=0;j<P;j++)
		{
			SSq += beta[j+1] * beta[j+1] / g[gMap[j]];
		}
		sig2 = 1/rgamma( (N+P)/2.0, 2.0/SSq );
		
		
		// Sample g
		AZERO(SSqG,nGs);
		for(j=0;j<P;j++)
		{
			SSqG[gMap[j]] += beta[j+1]*beta[j+1];
		}
		
		for(j=0;j<nGs;j++)
		{
			bTemp = SSqG[j]/sig2 + r[j]*r[j];
			g[j]  = 1/rgamma( (nParG[j]+1)/2.0, 2.0/bTemp);
		}
	
		Memcpy(chains + nPars*i,beta,P+1);	
		chains[nPars*i+P+1] = sig2;	
		Memcpy(chains + nPars*i + P + 2,g,nGs);	
	
		
	}
	
	UNPROTECT(2);
		
}

