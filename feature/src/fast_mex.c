#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <mex.h>

#include <fast.h>


#define VERBOSE 0

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
   const byte *input;   // The input image
   xy* cornerlist;
   int rows, cols;               // The dimensions of the image.
    int  numc, ncorner, i;
    int threshold;
    xy *lcorner;
    unsigned char *mcorner, *imcorner;
    char flag_nonmax;
    
    //xy (*ffast)(const byte*, int, int, int, int, int*);
    xy* (*ffast)();

    // inputs
    /*  create a pointer to the input matrix y */
    input = mxGetData(prhs[0]);;
    /*  get the dimensions of the matrix input y */
    rows = mxGetM(prhs[0]);
    cols = mxGetN(prhs[0]);
    
    if( !mxIsDouble(prhs[1]) || mxGetNumberOfElements(prhs[1])!=1) {
        mexErrMsgTxt("Input numc must be a scalar in {9,10,11,12}.");
    } else
        numc = (int)mxGetScalar(prhs[1]);
    
    if( !mxIsLogicalScalar(prhs[2])) {
        mexErrMsgTxt("Input flag_nonmax must be a logical scalar.");
    } else
        flag_nonmax = (int)(mxGetScalar(prhs[2]));

    if( !mxIsDouble(prhs[3]) || mxGetNumberOfElements(prhs[3])!=1 ) {
        mexErrMsgTxt("Input threshold must be scalar.");
    } else
        threshold = (int)mxGetScalar(prhs[3]);


    // select the FAST function
    switch(numc) {
        case 9:
        if(flag_nonmax) ffast = &fast9_detect_nonmax; 
        else            ffast = &fast9_detect;
        break;
        
        case 10:
        if(flag_nonmax) ffast = &fast10_detect_nonmax; 
        else            ffast = &fast10_detect;
        break;
        
        case 11: 
        if(flag_nonmax) ffast = &fast11_detect_nonmax; 
        else            ffast = &fast11_detect;
        break;
        
        case 12:
        if(flag_nonmax) ffast = &fast12_detect_nonmax; 
        else            ffast = &fast12_detect;
        break;
        
        default:
            mexErrMsgTxt("Input numc must be a scalar in {9,10,11,12}.");
    }
    
    // run
    lcorner = (*ffast)(input, rows, cols, rows, threshold, &ncorner);
    
    // output
    /*  set the output pointer to the output matrix */
    plhs[0] = mxCreateNumericMatrix(ncorner, 2, mxUINT8_CLASS, mxREAL);
    mcorner = mxGetData(plhs[0]);
    
    // store
    for (i=0; i<ncorner; i++) {
        mcorner[2*i] = lcorner[i].x;
        mcorner[2*i+1] = lcorner[i].y;
    }
    
    if(nlhs = 2) {
        plhs[1] = mxCreateLogicalMatrix(rows, cols);
        imcorner = mxGetData(plhs[1]);
        
        for (i=0; i<ncorner; i++)
            imcorner[lcorner[i].y * rows + lcorner[i].x] = 1;
    }
    
    free(lcorner);
}


