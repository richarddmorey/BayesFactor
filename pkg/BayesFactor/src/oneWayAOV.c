#include "BFPCL.h"



SEXP RgibbsOneWayAnova(SEXP yR, SEXP NR, SEXP JR, SEXP IR, SEXP rscaleR, SEXP iterationsR, SEXP progressR, SEXP pBar, SEXP callback, SEXP rho)
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
	
	gibbsOneWayAnova(yVec, N, J, sumN, whichJ, rscale, iterations, REAL(chainsR), REAL(CMDER), debug, progress, pBar, callback, rho);
	
	SET_VECTOR_ELT(returnListR, 0, chainsR);
    SET_VECTOR_ELT(returnListR, 1, CMDER);
	SET_VECTOR_ELT(returnListR, 2, debug);
	
	UNPROTECT(4);
	
	return(returnListR);
}



void gibbsOneWayAnova(double *y, int *N, int J, int sumN, int *whichJ, double rscale, int iterations, double *chains, double *CMDE, SEXP debug, int progress, SEXP pBar, SEXP callback, SEXP rho)
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
	
	int iOne=1, info;
	double dZero=0;
		

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
		logDet = matrixDet(B2temp,J,J,1, &info);
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
	
	UNPROTECT(2);
	PutRNGstate();
	
}
