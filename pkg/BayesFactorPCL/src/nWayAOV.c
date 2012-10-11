#include "BFPCL.h"


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

SEXP RjeffmlikeNWayAov(SEXP XtCXR, SEXP XtCyR, SEXP ytCyR, SEXP NR, SEXP PR, SEXP gR)
{
	int N = INTEGER_VALUE(NR);
	int P = INTEGER_VALUE(PR);
	double *XtCX = REAL(XtCXR);
	double *XtCy = REAL(XtCyR);
	double ytCy = REAL(ytCyR)[0];
	double *g = REAL(gR);
	
	SEXP ans;
	PROTECT(ans = allocVector(REALSXP,1));
	
	REAL(ans)[0] = jeffmlikeNWayAov(XtCX, XtCy, ytCy, N, P, g);
	
	UNPROTECT(1);
	return(ans);
}

double jeffSamplerNwayAov(double *samples, double *gsamples, int iters, double *XtCX, double *XtCy, double ytCy, int N, int P,int nGs, int *gMap, double *a, double *b, int progress, SEXP pBar, SEXP rho)
{
	int i=0,j=0;
	double avg = 0, *g1, g2[P];
	
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

		g1 = gsamples + i*nGs;
	
		for(j=0;j<nGs;j++)
		{
			g1[j]  = 1/rgamma(a[j],1/b[j]);
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



SEXP RjeffSamplerNwayAov(SEXP Riters, SEXP RXtCX, SEXP RXtCy, SEXP RytCy, SEXP RN, SEXP RP, SEXP RnGs, SEXP RgMap, SEXP Ra, SEXP Rb, SEXP progressR, SEXP pBar, SEXP rho)
{
	double *a,*b,*XtCy,*XtCX,ytCy,*samples,*avg;
	int iters,nGs,*gMap,N,P,progress;

	SEXP Rsamples,Rgsamples,Ravg,returnList;
	
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
	PROTECT(returnList = allocVector(VECSXP,3));
	PROTECT(Rgsamples = allocMatrix(REALSXP,nGs,iters));	

	samples = REAL(Rsamples);
	avg = REAL(Ravg);
	
	GetRNGstate();
	avg[0] = jeffSamplerNwayAov(samples, REAL(Rgsamples), iters, XtCX, XtCy, ytCy, N, P, nGs, gMap, a, b, progress, pBar, rho);
	PutRNGstate();
	
	SET_VECTOR_ELT(returnList, 0, Ravg);
    SET_VECTOR_ELT(returnList, 1, Rsamples);
    SET_VECTOR_ELT(returnList, 2, Rgsamples);
	
	UNPROTECT(4);
	
	return(returnList);
}

SEXP RGibbsNwayAov(SEXP Riters, SEXP Ry, SEXP RX, SEXP RXtX, SEXP RXty, SEXP RN, SEXP RP, SEXP RnGs, SEXP RgMap, SEXP Rr, SEXP progressR, SEXP pBar, SEXP rho)
{
	double *Xty,*XtX,*samples,*X, *y, *r;
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
	
	PROTECT(Rsamples = allocMatrix(REALSXP,iters, 3 + P + nGs));

	samples = REAL(Rsamples);
	
	GetRNGstate();
	GibbsNwayAov(samples, iters, y, X, XtX, Xty, N, P, nGs, gMap, r, progress, pBar, rho);
	PutRNGstate();
	
	
	UNPROTECT(1);
	
	return(Rsamples);
}

void GibbsNwayAov(double *chains, int iters, double *y, double *X, double *XtX, double *Xty, int N, int P, int nGs, int *gMap, double *r, int progress, SEXP pBar, SEXP rho)
{
	int i=0,j=0,k=0, nPars=3+P+nGs, P1Sq=(P+1)*(P+1), P1=P+1, iOne=1;
	double Sigma[P1Sq], SSq, oneOverSig2,dZero=0,dOne=1,dnegOne=-1;
	
	double g[nGs], beta[P+1], sig2, yTemp[N], SSqG[nGs],bTemp,nParG[nGs];
	//double betaSub[P], SigSub[P*P], logCMDE;
	
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
		
		//CMDE density
		// Removed for now...
    /*
    for(j=0;j<P;j++){
			betaSub[j] = beta[j+1];
			SigSub[j*P+j] = Sigma[(j+1)*(P+1) + (j+1)] - Sigma[j+1] * Sigma[j+1] / Sigma[0];
			for(k=0;k<j;k++){
				SigSub[j*P+k] = Sigma[(j+1)*(P+1) + (k+1)] - Sigma[j+1] * Sigma[k+1] / Sigma[0];
				SigSub[k*P+j] = SigSub[j*P+k];
			}
		}
		logCMDE = -P/2 * log(2 * M_PI) - 0.5 * matrixDet(SigSub, P, P, 1) - 0.5 * quadform(betaSub,SigSub,P,1,P);
		*/
    
    
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
	  //chains[nPars*i + P + nGs + 3] = logCMDE;
		
	}
	
	UNPROTECT(2);
		
}


