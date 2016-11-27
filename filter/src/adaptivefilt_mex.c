#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <mex.h>

#define		GW_ABS(a)       ((a) > 0 ? (a) : -(a))			//!<	Returns the absolute value a

#define A_(i,j) A[(i)+n1*(j)]
#define B_(i,j) B[(i)+n1*(j)]
#define I_(i,j) I[(i)+n1*(j)]
#define H_(i,j,k) H[(i)+p1*(j)+p1*p2*(k)]


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray*prhs[] ) 
{ 
    int n1, n2, p1, p2, m;
    int a1, a2, q1, q2;
    double *A, *H, *I, *B;
    int i, j, h, s1, s2;
    double v;
    
	// first argument : input array
	n1 = mxGetM(prhs[0]); 
	n2 = mxGetN(prhs[0]);
    A = mxGetPr(prhs[0]);

    // secong argument : input filters
    p1 = mxGetDimensions(prhs[1])[0];
    p2 = mxGetDimensions(prhs[1])[1];
    m = mxGetDimensions(prhs[1])[2];
    H = mxGetPr(prhs[1]);
    if( (p1%2)!=1 || (p2%2)!=1 )
        mexErrMsgTxt("Filters should be of odd size."); 

    // third argument : input index
    a1 = mxGetM(prhs[2]); 
    a2 = mxGetN(prhs[2]);
    I = mxGetPr(prhs[2]);
    if( a1!=n1 || a2!=n2 )
        mexErrMsgTxt("Array A and I should be of the same size."); 

    // output results
	plhs[0] = mxCreateDoubleMatrix(n1, n2, mxREAL);
	B = mxGetPr(plhs[0]);

    q1 = (p1-1)/2;
    q2 = (p2-1)/2;

    // perform filtering
    for( j=0; j<n2; ++j ){
        for( i=0; i<n1; ++i ) {
            // number of the filter
            h = (int) (I_(i,j)-1);
            B_(i,j) = 0;
            v = 0;   // normalization value
            for( s1=-q1; s1<=q1; ++s1 )
            for( s2=-q2; s2<=q2; ++s2 ) 
                if( (i+s1>=0) && (i+s1<n1) && (j+s2>=0) && (j+s2<n2)  ) {
                    B_(i,j) += A_( i+s1,j+s2 ) * H_( s1+q1,s2+q2,h );
                    v += GW_ABS(H_( s1+q1,s2+q2,h ));
                }
            B_(i,j) /= v;
        }
    }
}
