#include "BFPCL.h"


double jeffmlikeNWayAov(double *XtCX, double *priorX, double *XtCy, double ytCy, int N, int P, double *g, int incCont, double logDetPrX, int *colInfo)
{
	double *W,ldetS=0,ldetW,top,bottom1,bottom2,q;
	int i,j, Psqr=P*P, info=0;
  
	W = Memcpy(Calloc(Psqr,double),XtCX,Psqr);
	
  ldetS += -logDetPrX;
  for(i=0;i<incCont;i++){
    ldetS += log(g[i]);
    for(j=0;j<=i;j++){
      W[j + P*i] = W[j + P*i] + priorX[i+j*incCont]/g[i]; 
    }
  }
  
  for(i=incCont;i<P;i++)
	{
		ldetS += log(g[i]);
		W[i + P*i] = W[i + P*i] + 1/g[i];
	}

  ldetW = -matrixDet(W, P, P, 1, &info);
  InvMatrixUpper(W, P);
	internal_symmetrize(W, P);
	
	q = quadform(XtCy, W, P, 1, P);
	
	top = 0.5*ldetW;
	bottom1 = 0.5*(N-1) * log(ytCy - q);
	bottom2 = 0.5*ldetS;
	
	Free(W);
	
	return(top-bottom1-bottom2);	
}

SEXP RjeffmlikeNWayAov(SEXP XtCXR, SEXP priorXR, SEXP XtCyR, SEXP ytCyR, SEXP NR, SEXP PR, SEXP gR, SEXP incContR, SEXP logDetPrXR)
{
	int N = INTEGER_VALUE(NR);
	int P = INTEGER_VALUE(PR);
	double *XtCX = REAL(XtCXR);
  double *priorX = REAL(priorXR);
  double *XtCy = REAL(XtCyR);
	double ytCy = REAL(ytCyR)[0];
	double *g = REAL(gR);
  int incCont  = INTEGER_VALUE(incContR);
  double logDetPrX = REAL(logDetPrXR)[0];
  int cholInfo=0;
  
	SEXP ans;
	PROTECT(ans = allocVector(REALSXP,1));
	

	REAL(ans)[0] = jeffmlikeNWayAov(XtCX, priorX, XtCy, ytCy, N, P, g, incCont, logDetPrX, &cholInfo);
	
	UNPROTECT(1);
  
	return(ans);
}

double jeffSamplerNwayAov(double *samples, double *gsamples, int iters, double *XtCX, double *priorX, double *XtCy, double ytCy, int N, int P,int nGs, int *gMap, double *a, double *b, int incCont, int *badSamples, int progress, SEXP pBar, SEXP rho)
{
  int i=0,j=0, info=0, cholInfo=0,isFirst=1;
	double avg = 0, *g1, g2[P], *W;
  double logDetPrX=0;
  
  *badSamples = 0;
  
	// progress stuff
	SEXP sampCounter, R_fcall;
	int *pSampCounter;
    PROTECT(R_fcall = lang2(pBar, R_NilValue));
	PROTECT(sampCounter = NEW_INTEGER(1));
	pSampCounter = INTEGER_POINTER(sampCounter);

  if(incCont){
    logDetPrX = matrixDet(priorX, incCont, incCont, 1, &info);
  }
  
  
	for(i=0;i<iters;i++)
	{
    R_CheckUserInterrupt();
		
		//Check progress
		if(progress){
      if(!((i+1)%progress)){
			  pSampCounter[0]=i+1;
			  SETCADR(R_fcall, sampCounter);
			  eval(R_fcall, rho); //Update the progress bar
		  }
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
	  
    samples[i] = jeffmlikeNWayAov(XtCX, priorX, XtCy, ytCy, N, P, g2, incCont, logDetPrX,&cholInfo);
    if(cholInfo){
      *badSamples++;
		}else{
      if(isFirst)
		  {
			  avg = samples[i];
        isFirst=0;
		  }else{
			  avg = LogOnePlusExpX(samples[i]-avg)+avg;
			  //avg = LogOnePlusX(exp(avg-samples[i]))+samples[i];
		  }
		}

 } 
	
	UNPROTECT(2);
	return(avg-log(iters - *badSamples));
	
}

SEXP RjeffSamplerNwayAov(SEXP Riters, SEXP RXtCX, SEXP RpriorX, SEXP RXtCy, SEXP RytCy, SEXP RN, SEXP RP, SEXP RnGs, SEXP RgMap, SEXP Ra, SEXP Rb, SEXP RincCont, SEXP progressR, SEXP pBar, SEXP rho)
{
  double *a,*b,*XtCy,*XtCX,ytCy,*samples,*avg, *priorX;
	int iters,nGs,*gMap,N,P,progress,incCont, badSamples=0;

	SEXP Rsamples,Rgsamples,Ravg,returnList,RbadSamples;
	
	iters = INTEGER_VALUE(Riters);
	XtCX = REAL(RXtCX);
  priorX = REAL(RpriorX);
  XtCy = REAL(RXtCy);
	ytCy = REAL(RytCy)[0];
	nGs = INTEGER_VALUE(RnGs);
	gMap = INTEGER_POINTER(RgMap);
	a = REAL(Ra);
	b = REAL(Rb);
	N = INTEGER_VALUE(RN);
	P = INTEGER_VALUE(RP);
	progress = INTEGER_VALUE(progressR);
  incCont = INTEGER_VALUE(RincCont);

	PROTECT(Rsamples = allocVector(REALSXP,iters));
	PROTECT(Ravg = allocVector(REALSXP,1));
  PROTECT(RbadSamples = allocVector(REALSXP,1));
	PROTECT(returnList = allocVector(VECSXP,3));
	PROTECT(Rgsamples = allocMatrix(REALSXP,nGs,iters));	

	samples = REAL(Rsamples);
	avg = REAL(Ravg);
	
	GetRNGstate();
	avg[0] = jeffSamplerNwayAov(samples, REAL(Rgsamples), iters, XtCX, priorX,XtCy, ytCy, N, P, nGs, gMap, a, b, incCont, &badSamples, progress, pBar, rho);
	PutRNGstate();
	
  REAL(RbadSamples)[0] = badSamples * 1.0;
  
	SET_VECTOR_ELT(returnList, 0, Ravg);
  SET_VECTOR_ELT(returnList, 1, Rsamples);
  SET_VECTOR_ELT(returnList, 2, Rgsamples);
	SET_VECTOR_ELT(returnList, 3, RbadSamples);
  
	UNPROTECT(5);
	
	return(returnList);
}


double importanceSamplerNwayAov(double *samples, double *qsamples, int iters, double *XtCX, double *priorX, double *XtCy, double ytCy, int N, int P,int nGs, int *gMap, double *a, double *b, double *mu, double *sig, int incCont, int *badSamples, int progress, SEXP pBar, SEXP rho)
{
  int i=0,j=0, info=0, cholInfo=0, isFirst=1;
	double avg = 0, *q1, g2[P], sumq=0, sumdinvgamma=0, sumdnorm=0;
	double logDetPrX=0;
  
  *badSamples = 0;
  
	// progress stuff
	SEXP sampCounter, R_fcall;
	int *pSampCounter;
    PROTECT(R_fcall = lang2(pBar, R_NilValue));
	PROTECT(sampCounter = NEW_INTEGER(1));
	pSampCounter = INTEGER_POINTER(sampCounter);

  if(incCont){
    logDetPrX = matrixDet(priorX, incCont, incCont, 1, &info);
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

		q1 = qsamples + i*nGs;
    
	  sumq = 0;
    sumdinvgamma = 0;
    sumdnorm = 0;
		for(j=0;j<nGs;j++)
		{
			q1[j]  = rnorm(mu[j],sig[j]);
		  sumq += q1[j];
      sumdinvgamma += dgamma(exp(-q1[j]), a[j], 1/b[j], 1) - 2*q1[j];
      sumdnorm += dnorm(q1[j], mu[j], sig[j], 1);
    }
    
		for(j=0;j<P;j++)
		{
			g2[j] = exp(q1[gMap[j]]);
		}
	
    samples[i] = jeffmlikeNWayAov(XtCX, priorX, XtCy, ytCy, N, P, g2, incCont, logDetPrX, &cholInfo) + sumq + sumdinvgamma - sumdnorm;;
		if(cholInfo){
      *badSamples++;
		}else{
    	if(isFirst)
  		{
  			avg = samples[i];
        isFirst=0;
  		}else{
  			avg = LogOnePlusExpX(samples[i]-avg)+avg;
  			//avg = LogOnePlusX(exp(avg-samples[i]))+samples[i];
  		}    
		}
	}
	
	UNPROTECT(2);
	
	return(avg-log(iters-*badSamples));
	
}

SEXP RimportanceSamplerNwayAov(SEXP Riters, SEXP RXtCX, SEXP RpriorX, SEXP RXtCy, SEXP RytCy, SEXP RN, SEXP RP, SEXP RnGs, SEXP RgMap, SEXP Ra, SEXP Rb, SEXP Rmu, SEXP Rsig, SEXP RincCont, SEXP progressR, SEXP pBar, SEXP rho)
{
	double *a,*b,*XtCy,*XtCX,ytCy,*samples,*avg, *mu, *sig, *priorX;
	int iters,nGs,*gMap,N,P,progress,incCont,badSamples=0;

	SEXP Rsamples,Rgsamples,Ravg,returnList,RbadSamples;
	
	iters = INTEGER_VALUE(Riters);
	XtCX = REAL(RXtCX);
  priorX = REAL(RpriorX);
  XtCy = REAL(RXtCy);
	ytCy = REAL(RytCy)[0];
	nGs = INTEGER_VALUE(RnGs);
	gMap = INTEGER_POINTER(RgMap);
	a = REAL(Ra);
	b = REAL(Rb);
	N = INTEGER_VALUE(RN);
	P = INTEGER_VALUE(RP);
  mu = REAL(Rmu);
	sig = REAL(Rsig);
  incCont = INTEGER_VALUE(RincCont);

  progress = INTEGER_VALUE(progressR);
	
	PROTECT(Rsamples = allocVector(REALSXP,iters));
	PROTECT(Ravg = allocVector(REALSXP,1));
  PROTECT(RbadSamples = allocVector(REALSXP,1));
  PROTECT(returnList = allocVector(VECSXP,3));
	PROTECT(Rgsamples = allocMatrix(REALSXP,nGs,iters));	

	samples = REAL(Rsamples);
	avg = REAL(Ravg);
	
	GetRNGstate();
	avg[0] = importanceSamplerNwayAov(samples, REAL(Rgsamples), iters, XtCX, priorX, XtCy, ytCy, N, P, nGs, gMap, a, b, mu, sig, incCont, &badSamples, progress, pBar, rho);
	PutRNGstate();
	
  REAL(RbadSamples)[0] = badSamples * 1.0;

	SET_VECTOR_ELT(returnList, 0, Ravg);
  SET_VECTOR_ELT(returnList, 1, Rsamples);
  SET_VECTOR_ELT(returnList, 2, Rgsamples);
  SET_VECTOR_ELT(returnList, 3, RbadSamples);
 
  
	UNPROTECT(5);
	
	return(returnList);
}

SEXP RGibbsNwayAov(SEXP Riters, SEXP Ry, SEXP RX, SEXP RXtX, SEXP RpriorX, SEXP RXty, SEXP RN, SEXP RP, SEXP RnGs, SEXP RgMap, SEXP Rr, SEXP RincCont, SEXP RignoreCols, SEXP Rthin, SEXP progressR, SEXP pBar, SEXP rho)
{
	double *a,*b,*Xty,*XtX,*samples,*X, *y, *r,*priorX;
	int iters,nGs,*gMap,N,P,progress,incCont,*ignoreCols, nOutputPars=0,thin,i;
  int effectiveIterations;
  
	SEXP Rsamples;
	
	iters = INTEGER_VALUE(Riters);
	X = REAL(RX);
	y = REAL(Ry);
	XtX = REAL(RXtX);
  priorX = REAL(RpriorX);
  Xty = REAL(RXty);
	nGs = INTEGER_VALUE(RnGs);
	gMap = INTEGER_POINTER(RgMap);
	r = REAL(Rr);
	N = INTEGER_VALUE(RN);
	P = INTEGER_VALUE(RP);
	progress = INTEGER_VALUE(progressR);
  incCont = INTEGER_VALUE(RincCont); 
  ignoreCols = INTEGER_POINTER(RignoreCols);  
  thin = INTEGER_VALUE(Rthin); 
  
  effectiveIterations = iters / thin;
  
  for(i=0;i<(2 + P + nGs);i++){
    if(!ignoreCols[i]) nOutputPars++;
  }
  
	PROTECT(Rsamples = allocMatrix(REALSXP,effectiveIterations, nOutputPars));
  
	samples = REAL(Rsamples);
	  
	GetRNGstate();
	GibbsNwayAov(samples, iters, y, X, XtX, priorX, Xty, N, P, nGs, gMap, r, incCont, ignoreCols, nOutputPars, thin, progress, pBar, rho);
	PutRNGstate();
	
	UNPROTECT(1);
	
	return(Rsamples);
}

void GibbsNwayAov(double *chains, int iters, double *y, double *X, double *XtX, double *priorX, double *Xty, int N, int P, int nGs, int *gMap, double *r, int incCont, int *ignoreCols, int nOutputPars, int thin, int progress, SEXP pBar, SEXP rho)
{
	int i=0,j=0,k=0, nPars=2+P+nGs, P1Sq=(P+1)*(P+1), P1=P+1, iOne=1,nParG[nGs];
	double *g1, Sigma[P1Sq], SSq, oneOverSig2,dZero=0,dOne=1,dnegOne=-1;
	double oneChainCol[nPars];
  
	double g[nGs], beta[P+1], sig2, yTemp[N], SSqG[nGs],bTemp;
	
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
    // Remember, continuous g is always last
		// Sample beta
		Memcpy(Sigma,XtX,P1Sq);
    for(j=0;j<(P+1);j++)
		{
			// Diagonal
      if(j>incCont){
        Sigma[ j + (P+1)*j] = (Sigma[ j + (P+1)*j] + 1/g[gMap[j-1]])/sig2;
      }else if(j>0){ // only way this can happen is if incCont>0
	      Sigma[ j + (P+1)*j] = (Sigma[ j + (P+1)*j] + priorX[(j-1) + incCont*(j-1)]/g[nGs-1])/sig2;
			}else{
				Sigma[ j + (P+1)*j] = Sigma[ j + (P+1)*j]/sig2;
			}
			for(k=j+1;k<(P+1);k++)
			{	
        if(k<=incCont && j>0){ // only way this can happen is if incCont>0
          Sigma[j + (P+1)*k] = (Sigma[j + (P+1)*k] + priorX[(j-1) + incCont*(k-1)]/g[nGs-1])/sig2;
        }else{
          Sigma[j + (P+1)*k] = Sigma[j + (P+1)*k]/sig2;
        }
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
    
    // continuous g is last
    if(incCont){
      SSq += quadform2(beta+1,priorX,incCont,1,incCont)/g[nGs - 1];
    }    
    for(j=incCont;j<P;j++)
		{
			SSq += beta[j+1] * beta[j+1] / g[gMap[j]];
		}
    sig2 = 1/rgamma( (N+P)/2.0, 2.0/SSq );
    
		// Sample g
		AZERO(SSqG,nGs);
		
    // continuous g is last
    if(incCont){
      SSqG[nGs - 1] += quadform2(beta+1,priorX,incCont,1,incCont);
    }
    
    for(j=incCont;j<P;j++)
		{
      SSqG[gMap[j]] += beta[j+1]*beta[j+1];
		}
		
		for(j=0;j<nGs;j++)
		{
      bTemp = SSqG[j]/sig2 + r[j]*r[j];
      g[j]  = 1/rgamma( (nParG[j]+1)/2.0, 2.0/bTemp);
    }
	
  	Memcpy(oneChainCol, beta,P+1);	
		oneChainCol[P+1] = sig2;	
		Memcpy(oneChainCol + P + 2,g,nGs);	
    
    // Copy only requested columns to the final answer
    if(!(i % thin)){
      k=0;
      for(j=0;j<(2+P+nGs);j++){
        if(!ignoreCols[j]) chains[k++ + nOutputPars* (i/thin)] = oneChainCol[j];
      }
    } 	
	 
   } // end MCMC iterations
	
	UNPROTECT(2);
		
}


