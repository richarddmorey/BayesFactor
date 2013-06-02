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
