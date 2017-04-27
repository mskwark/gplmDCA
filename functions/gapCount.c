#include "mex.h"

/*
This is a reworked version of the pseudolikelihood objective contained in Mark Schmidt's thesis code (http://www.di.ens.fr/~mschmidt/Software/thesis.html). The copyright conditions for the original code is included below.
---------------------------------

Copyright 2005-2012 Mark Schmidt. All rights reserved.

Redistribution and use in source and binary forms, with or without modification, are
permitted provided that the following conditions are met:

   1. Redistributions of source code must retain the above copyright notice, this list of
      conditions and the following disclaimer.

   2. Redistributions in binary form must reproduce the above copyright notice, this list
      of conditions and the following disclaimer in the documentation and/or other materials
      provided with the distribution.

THIS SOFTWARE IS PROVIDED BY THE AUTHOR ''AS IS'' AND ANY EXPRESS OR IMPLIED
WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND
FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE AUTHOR OR
CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON
ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/

/* This gives the position of Gamma_{2}^sigma 
 * in the G array IN C INDEXING */

/* MEX FUNCTION */
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    /* variables */
    int N,M,*y,i,j,m;

    /* input */
    y = (int*)mxGetPr(prhs[0]);
    N=mxGetDimensions(prhs[0])[1];
    M=mxGetDimensions(prhs[0])[0];

    /* output */
    double* gapStretches;
    plhs[0] = mxCreateDoubleMatrix(N, N, mxREAL);
    gapStretches=mxGetPr(plhs[0]);
    double*  gapStretchesHist;
    plhs[1] = mxCreateDoubleMatrix(N,1,mxREAL);
    gapStretchesHist=mxGetPr(plhs[1]);
    for(i=0; i<N; i++)
    	for(j=0; j<N; j++)
		gapStretches[i+N*j]=0;
    for(m=0; m<M; m++){
		int k=0;
    		for(i=0; i<N; i++){
			for(j=i; j<N; j++){	
					if(y[m+j*M]==0){
						gapStretches[i+(j-i)*N]++;
					}
					else{
						break;
					}
			}
			if(y[m+i*M]==0 && i != N-1)
				k++;
			else{
				if(k>0)
				{
					gapStretchesHist[k-1]++;
					k=0;
				}
			}
				
		}
    }
    return;
}


