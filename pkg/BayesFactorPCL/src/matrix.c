#include "BFPCL.h"



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
