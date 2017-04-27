#include <math.h>
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

/* This gives the position of Gamma_{1}^sigma 
* in the G array IN C INDEXING */

int GindStart(int sigma,int N)
{
	return (sigma-1)*N-(sigma-1)*((sigma-1)+1)/2 +sigma-1;
}



/* MEX FUNCTION */
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
/* variables */
char param;
int j, i, s, t, n, nInstances, nNodes, nStates,
    *y, y1, y2, *rint,r,M, *lH, *rH;
double *weights, *grad1, *grad2, *grad3, *logPot, *z, *fval, *h_r, *J_r, *nodeBel, *lambdas, *G, *pM;

/* input */
y = (int*)mxGetPr(prhs[0]);
weights = mxGetPr(prhs[1]);
h_r = mxGetPr(prhs[2]);
J_r = mxGetPr(prhs[3]);
lambdas = mxGetPr(prhs[4]);
rint = (int*)mxGetPr(prhs[5]);
G=mxGetPr(prhs[6]);
pM=mxGetPr(prhs[7]);
M=(int) *pM;
lH = (int*)mxGetPr(prhs[8]);
rH = (int*)mxGetPr(prhs[9]);

/*int Gsize=mxGetDimensions(prhs[6])[0];
mexPrintf("Dim: %d %d %d %d %d %d", mxGetDimensions(prhs[8])[0], mxGetDimensions(prhs[8])[1],mxGetDimensions(prhs[9])[0],mxGetDimensions(prhs[9])[1], Gsize); */

/* compute sizes */
nNodes = mxGetDimensions(prhs[0])[1];
nStates = mxGetDimensions(prhs[2])[1];
nInstances = mxGetDimensions(prhs[0])[0];
float fM;
float fnNodes;
fM=(float) M;
fnNodes = (float) nNodes;
int nrGapParam;
nrGapParam=(int) fM*(fnNodes-(fM+1)/2+1);

/* allocate memory */
logPot = mxCalloc(nStates, sizeof(double));
z = mxCalloc(1, sizeof(double));
nodeBel = mxCalloc(nStates, sizeof(double));

int length;
int index;

/* output */
plhs[0] = mxCreateDoubleMatrix(1, 1, mxREAL);
fval = mxGetPr(plhs[0]);
*fval = 0;
plhs[1] = mxCreateDoubleMatrix(nStates, 1, mxREAL);
grad1 = mxGetPr(plhs[1]);
plhs[2] = mxCreateDoubleMatrix(nStates*nStates*(nNodes-1), 1, mxREAL);
grad2 = mxGetPr(plhs[2]);
plhs[3] = mxCreateDoubleMatrix(nrGapParam,1,mxREAL);
grad3= mxGetPr(plhs[3]);
for(i=0; i<nrGapParam; i++)
{
	grad3[i]=0.0;
}
	
r=rint[0]-1;

for(i=0;i < nInstances;i++) {
/*Some notes on variable names:
logPot contains, for the current sequence i, the exponentials in the pseudolikelihood: logPot(s)=h_r(s)+sum_{j!=r}J_{rj}(s,sigma^(i)_j).
nodeBel is the conditional probability P(sigma_r=s|sigma_{\r}=sigma^(i)_{\r}), i.e., nodeBel(s) = e^[ logPot(s) ] / sum_l e^[ logPot(l) ].
z is the denominator of nodeBel.*/

for(s=0;s < nStates;s++) {
    logPot[s] = h_r[s];

}

for(n = 0;n < nNodes;n++) {
    if(n!=r) {
	y2 = y[i + nInstances*n];       
	for(s=0; s<nStates; s++) {
	    logPot[s] += J_r[s+nStates*(y2+nStates*(n-(n>r)))];               
	}
    }
    
}

	/* Add GAP parameters */
	/* Restitute r by a gap */
	length=1;
	index=(r+1);
	if(r<nNodes-1)
		length=length+rH[i+nInstances*(r+1)];
	
	if(r>0){
		length=length+lH[i+nInstances*(r-1)];
		index=index-lH[i+nInstances*(r-1)];
	}
	logPot[0] += G[GindStart(length,nNodes)+index-1];

	/* Restitute r by non-gap */
	/* Look if there is now a gap to the right */
	if(r<nNodes-1)
		if(rH[i+nInstances*(r+1)] != 0)
			{
				length=rH[i+nInstances*(r+1)];
				index = (r+1)+1;
				for(s=1;s<nStates;s++)
					logPot[s] += G[GindStart(length,nNodes)+index-1];
			}
	/* Look if there is a gap to the left */
	if(r>0)
		if(lH[i+nInstances*(r-1)] != 0)
			{	
				length=lH[i+nInstances*(r-1)];
				index=(r+1)-length;
				for(s=1;s<nStates;s++)
					logPot[s] += G[GindStart(length,nNodes)+index-1];
			}


        z[0] = 0;
        for(s = 0; s < nStates; s++) {
            z[0] += exp(logPot[s]);
        }
        *fval -= weights[i]*logPot[y[i+nInstances*r]];
        *fval += weights[i]*log(z[0]);
 
        
        
	/*Gradient:*/
        
 
        for(s = 0; s < nStates; s++) {
            nodeBel[s] = exp(logPot[s] - log(z[0]));
        }
                       
        y1 = y[i + nInstances*r]; 
        grad1[y1] -= weights[i]*1;
            
        for(s=0; s < nStates; s++) {
            grad1[s] += weights[i]*nodeBel[s];
        }
   
        for(n=0;n<nNodes;n++) {
            if(n!=r) {
                y2 = y[i + nInstances*n];

                grad2[y1+nStates*(y2+nStates*(n-(n>r)))] -= weights[i];                

                for(s=0;s<nStates;s++) {
                   grad2[s+nStates*(y2+nStates*(n-(n>r)))] += weights[i]*nodeBel[s];
                }	
            }
        }


	/* Gap gradient */
	/* Restitute r by a gap */
	length=1;
	index=(r+1);
	if(r<nNodes-1)
		length=length+rH[i+nInstances*(r+1)];
	
	if(r>0){
		length=length+lH[i+nInstances*(r-1)];
		index=index-lH[i+nInstances*(r-1)];
	}
	grad3[GindStart(length,nNodes)+index-1] += weights[i]*nodeBel[0];
	if(y[i+nInstances*r]==0)
		grad3[GindStart(length,nNodes)+index-1] -= weights[i];

	/* Restitute r by non-gap */
	/* Look if there is now a gap to the right */
	if(r<nNodes-1)
		if(rH[i+nInstances*(r+1)] != 0)
			{
				length=rH[i+nInstances*(r+1)];
				index = (r+1)+1;
				for(s=1;s<nStates;s++)
					grad3[GindStart(length,nNodes)+index-1] +=weights[i]*nodeBel[s];
				if(y[i+nInstances*r] != 0)
					grad3[GindStart(length,nNodes)+index-1] -=weights[i];

			}
	/* Look if there is a gap to the left */
	if(r>0)
		if(lH[i+nInstances*(r-1)] != 0)
			{	
				length=lH[i+nInstances*(r-1)];
				index=(r+1)-length;
				for(s=1;s<nStates;s++)
					grad3[GindStart(length,nNodes)+index-1] +=weights[i]*nodeBel[s];
				if(y[i+nInstances*r] != 0)
					grad3[GindStart(length,nNodes)+index-1] -=weights[i];

			}




	}

  
    
    
    
    /*Add contributions from R_l2*/

    for(s = 0; s < nStates; s++) {
        *fval += lambdas[0]*h_r[s]*h_r[s];
        grad1[s] += lambdas[0]*2*h_r[s]; 
    }
    for(j=0; j<nrGapParam; j++)
    	{
    		*fval += lambdas[2]*G[j]*G[j];
    		grad3[j] += lambdas[2]*2*G[j];
    	}

    for(n = 0;n < nNodes;n++) {
        if(n!=r) {
            for(s = 0; s < nStates; s++) {
                for(t = 0; t < nStates; t++) {                   
                        *fval += lambdas[1]*J_r[s+nStates*(t+nStates*(n-(n>r)))]*J_r[s+nStates*(t+nStates*(n-(n>r)))];
                        grad2[s+nStates*(t+nStates*(n-(n>r)))] += lambdas[1]*2*J_r[s+nStates*(t+nStates*(n-(n>r)))];                    
                }
            }
        }
    }
             
    mxFree(logPot);
    mxFree(z);
    mxFree(nodeBel);
    return;
}


