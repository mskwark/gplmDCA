#include "mex.h"

/* MEX FUNCTION */
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    /* variables */
    int N,M,*y,m,n;

    /* input */
    y = (int*)mxGetPr(prhs[0]);
    N=mxGetDimensions(prhs[0])[1];
    M=mxGetDimensions(prhs[0])[0];

    /* output */
    double* leftGapMat;
    plhs[0] = mxCreateDoubleMatrix(M, N, mxREAL);
    leftGapMat=mxGetPr(plhs[0]);
    double* rightGapMat;
    plhs[1] = mxCreateDoubleMatrix(M,N,mxREAL);
    rightGapMat=mxGetPr(plhs[1]);
    int k1; 
    int k2;
    int n2;
    for(m=0; m<M; m++){
    	k1=0;
	k2=0;
    	for(n=N-1; n>=0; n--){
			if(y[m+M*n]==0)
				k1++;
			else
			   	k1=0;
			rightGapMat[m+M*n]=k1;
			n2=N-n-1;
			if(y[m+M*n2]==0)
				k2++;
			else
				k2=0;
			leftGapMat[m+M*n2]=k2;

		}
    }
    return;
    }


