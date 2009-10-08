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


double sampG(double x, double *data, int prior, int method, int logg, double *control);
double dG(double x, double *data, int prior, int logg);
SEXP RsampG(SEXP x, SEXP Niters, SEXP data, SEXP prior, SEXP method, SEXP logg, SEXP control, SEXP progress, SEXP pBar, SEXP rho);
SEXP RdG(SEXP x, SEXP data, SEXP prior, SEXP logg);
double gprior1(double x);


// Prior 1: standard exponential
double gprior1(double x)
{
	double lambda = 1;
	return( log(lambda) - lambda * x) ;
}

SEXP RsampG(SEXP x, SEXP Niters, SEXP data, SEXP prior, SEXP method, SEXP logg, SEXP control, SEXP progress, SEXP pBar, SEXP rho)
{
	double cx = REAL(x)[0], *cdata = REAL(data), *ccontrol = REAL(control);
	int cprior = INTEGER_VALUE(prior), clogg = INTEGER_VALUE(logg);
	int cmethod = INTEGER_VALUE(method), cNiters = INTEGER_VALUE(Niters);
	int i=0;
	 
	
	// Progress bar stuff
	SEXP sampCounter, R_fcall;
	int *pSampCounter, iProgress = INTEGER_VALUE(progress);
    PROTECT(R_fcall = lang2(pBar, R_NilValue));
	PROTECT(sampCounter = NEW_INTEGER(1));
	pSampCounter = INTEGER_POINTER(sampCounter);
	
	
	if(length(data)!=3)
		error("Length of DATA should be 3.");
	
	SEXP ans;
	PROTECT(ans = allocVector(REALSXP, cNiters));
	double *pAns = REAL(ans);
	
	pAns[0]=cx;
	
	for(i=1;i<cNiters;i++)
	{
		//Check progress
		if(iProgress && !((i+1)%iProgress)){
			pSampCounter[0]=i+1;
			SETCADR(R_fcall, sampCounter);
			eval(R_fcall, rho); //Update the progress bar
		}

	pAns[i] = sampG(pAns[i-1], cdata, cprior, cmethod, clogg, ccontrol);
	}
	
	UNPROTECT(3);
	return(ans);

}

double sampG(double x, double *data, int prior, int method, int logg, double *control)
{
	double ans=x, p, cand, v;
	GetRNGstate();
	
	
	
	switch(method)
	{
		// random walk Metropolis Hastings. control[0] is standard dev.
		case 1:
			p = rnorm(0,control[0]);
			cand = x + p;
			v = log(runif(0,1));
			if(v < dG(cand, data, prior, logg) - dG(x, data, prior, logg))
				ans = cand;
			break;
		
		default:
			error("Sampling method not implemented.");
			break;
	}

	PutRNGstate();
	return(ans);
}

// data[0] = SSE, data[1]=ybar^2, data[2] = N
double dG(double x, double *data, int prior, int logg)
{
	double SSE = data[0], ybar2 = data[1], N=data[2],ans;
	
	if(logg){
		x = exp(x);
	}else if(x<0){
		return(log(0));
	}
	
	ans = log( SSE + ybar2 * (N/(x*N + 1)) );

	switch(prior)
	{
		case 1:
			ans = -N/2*ans  + gprior1(x) - .5*log(N*x + 1);
			break;
		
		default:
			error("Prior not implemented.");
			break;
	}
	
	if(logg) ans = ans + log(x);
	return(ans);
}

SEXP RdG(SEXP x, SEXP data, SEXP prior, SEXP logg)
{
	double *cx = REAL(x), *cdata = REAL(data);
	int cprior = INTEGER_VALUE(prior), clogg = INTEGER_VALUE(logg);
	int i=0, points=length(x);
	
	if(length(data)!=3)
		error("Length of DATA should be 3.");
	
	SEXP ans;
	PROTECT(ans = allocVector(REALSXP, points));
	double *pAns = REAL(ans);
	
	
	for(i=0;i<points;i++)
	{
		pAns[i] = dG(cx[i], cdata, cprior, clogg);
	}
	
	UNPROTECT(1);
	return(ans);
}



