/***********************************************************************
 * mexrepval.c : C mex file
 *
 *   Input: nzA, omega, index
 *
 *   Output: z(index(i)+1:index(i+1)) = nzA(index(i)+1:index(i+1))*omega(i)
 ***********************************************************************/

#include <mex.h>
#include <math.h>
#include <matrix.h>


/**********************************************************
 * MAIN MEX FUNCTION
 ***********************************************************/
void mexFunction(int nlhs, mxArray  *plhs[],
        int nrhs, const mxArray  *prhs[] )
        
{    double   *nzA,  *omega,  *index;
     int      i,j,m,sizeA;
     double   *z;;   /*out put*/
     
     /* CHECK FOR PROPER NUMBER OF ARGUMENTS */
     
     if (nrhs != 3){
         mexErrMsgTxt("repval: requires 3 input argument."); }
     if (nlhs > 1){
         mexErrMsgTxt("repval: requires 1 output argument."); }   

     /***** assign pointers *****/
     nzA   = mxGetPr(prhs[0]);
     omega = mxGetPr(prhs[1]);
     index = mxGetPr(prhs[2]);
     
     sizeA = mxGetM(prhs[0]);
     m     = mxGetM(prhs[1]);
     
     /***** create return argument *****/
     plhs[0] = mxCreateDoubleMatrix(sizeA,1,mxREAL);
     z  = mxGetPr(plhs[0]);
     
     for ( i = 0; i < m; i++ ) 
     { 
         for (j = index[i]; j < index[i+1]; j ++) 
         {
            z[j] = nzA[j] * omega[i]; 
         }
     }     
   
     return;
}
/**********************************************************/
