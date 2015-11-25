#include "math.h"
#include "mex.h"
#include "matrix.h"

/* mxArray* X = reconstructSpike(mxArray* A,mxArray* t,mxArray* Templates,mxArray* trange,int onewaveformlength,int numberchanns)
 *
 * Takes list of amplitude scale factors, times, and templates as inputs 
 * to produce a reconstructed spike for subtraction / display.
 *
 * Original MATLAB code:
 *
 * function X = reconstructSpike(A,t,Templates,trange)
 *   global onewaveformlength; global numberchanns;
 *   X = zeros(length(trange),numberchanns);    
 *   for i=1:length(A)
 *       shiftspike = interp1((1:onewaveformlength)-floor(onewaveformlength/2),Templates(:,:,i),trange-t(i),'nearest',0);
 *       X = X + A(i) * shiftspike;
 *   end
 * end
 */


void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{
    // Read input variables
    mxArray* AData; double *A;
    mxArray* tData; double *t;
    double *Templates;
    mxArray* trangeData; int *trange;
    int onewaveformlength;
    int numberchanns;
    
    AData = prhs[0];
    tData = prhs[1];
    trangeData = prhs[3];
    
    A = mxGetPr(AData);
    t = mxGetPr(tData);
    Templates = mxGetPr(prhs[2]);
    trange = (int*) mxGetPr(trangeData);
        
    onewaveformlength = (int) mxGetScalar(prhs[4]);
    numberchanns = (int) mxGetScalar(prhs[5]);

    int midwaveform = floor(onewaveformlength/2);
    
    // Main loop:
    int tlen, Alen;    
    tlen = mxGetNumberOfElements(trangeData);
    Alen = mxGetNumberOfElements(AData);
    
    mxArray* X = mxCreateDoubleMatrix(tlen, numberchanns, mxREAL);
    double* XValue = mxGetPr(X);
    
    int i,j,k;
    int ti, timeL, timeR;
    for (i=0; i<Alen; i++)
    {        
        ti = (int) t[i];
        if (tlen-ti > onewaveformlength-midwaveform) 
            timeR = onewaveformlength - midwaveform;
        else
            timeR = tlen-ti;
        if (ti > midwaveform)
            timeL = midwaveform;
        else
            timeL = ti;
       
        for (j=-timeL; j<timeR; j++)
        {
            for (k=0; k<numberchanns; k++)
            {
                XValue[tlen*k + (ti + j)] = XValue[tlen*k + (ti + j)] + A[i]*Templates[onewaveformlength*numberchanns*i + onewaveformlength*k + midwaveform+j];
            }
        }
        
    }
            
    // Assign output
    plhs[0] = X;
    
    return;
}

