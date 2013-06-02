#include "BFPCL.h"



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
