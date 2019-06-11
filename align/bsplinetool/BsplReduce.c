/*******************************************************************************

	This file is a mex file to wrap C functions from the
	spline pyramid package at http://bigwww.epfl.ch/ and its
	trivial 3D extension for REDUCE operator to be
	used easily in MATLAB environment.

	Coded by Se Young Chun, May 29, 2007, the University of Michigan

*******************************************************************************/


#include "mex.h"
#include "matrix.h"

#include "BIG/pyramids/configs.h"
#include "BIG/pyramids/pyramidtools.h"

#include "BIG/pyramidtools3D.h"

/* --- Defines --- */
#define MAXF 200L	       /* Maximum size of the filter */

#define SPLINE	  "Spline"	  /* Spline filter (l2-norm) */
#define SPLINE_L2       "Spline L2"       /* Spline filter (L2-norm) */
#define SPLINE_CENT     "Centered Spline" 
#define SPLINE_CENT_L2  "Centered Spline L2"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	int m, n, l;
	mwSize	dims[3];
 	float   *data1, *data2;
	double  spline;
	double  h[MAXF];     
	long    nh;      
	double  g[MAXF];     /* Coefficients of the reduce filter */
	long    ng;	  /* Number of coefficients of the reduce filter */
	short   IsCentered;  /* Equal TRUE if the filter is a centered spline, 
				      FALSE otherwise */

	if (nlhs == 0 && nrhs == 1 && mxIsChar(prhs[0])) // check
		return;

	if (nlhs != 1)
	{
		mexPrintf("need nlhs=1\n");
       		return;
	}

	if (nrhs == 2)
	{
		spline = mxGetScalar(prhs[1]);
	}
	else if (nrhs == 1)
	{
		spline = 3;
	}
	else
	{
		mexPrintf("REDUCED = BsplReduce(Val, deg);\n");
		mexPrintf("[in]\n");
		mexPrintf("\tVal     : Values\n");
		mexPrintf("\tdeg     : (opt) basis degree {default: 3}\n");
		mexPrintf("[out]\n");
		mexPrintf("\tREDUCED : Reduced values of half resolution\n");
		return;
	}

    	/* Retrieve the input data */
    	data1 = mxGetData(prhs[0]);

	if (mxGetNumberOfDimensions(prhs[0]) == 2)
	{
		m = (int)mxGetDimensions(prhs[0])[0];
		n = (int)mxGetDimensions(prhs[0])[1];
		l = 1;
		dims[0] = m/2;
		if (dims[0] < 1L) dims[0] = 1L;

		dims[1] = n/2;
		if (dims[1] < 1L) dims[1] = 1L;
	}
	else if (mxGetNumberOfDimensions(prhs[0]) == 3)
	{
		m = (int)mxGetDimensions(prhs[0])[0];
		n = (int)mxGetDimensions(prhs[0])[1];
		l = (int)mxGetDimensions(prhs[0])[2];
		dims[0] = m/2;
		if (dims[0] < 1L) dims[0] = 1L;

		dims[1] = n/2;
		if (dims[1] < 1L) dims[1] = 1L;

		dims[2] = l/2;
		if (dims[2] < 1L) dims[2] = 1L;
	}
	else 
	{
	       	mexErrMsgTxt("Input data must be 2D or 3D\n");
		return;
	}

    	if (mxIsSingle(prhs[0]) != 1)
       	{
	       	mexErrMsgTxt("Input argument must be single type\n");
	       	return;
       	}

    	/* Create an mxArray for the output data */
	plhs[0] = mxCreateNumericArray(mxGetNumberOfDimensions(prhs[0])
			, dims, mxSINGLE_CLASS, mxREAL);

	/* Retrieve the output data */
    	data2 = mxGetData(plhs[0]);

	/* Get the filter coefficients for the Spline (order = 3) filter*/
	if (GetPyramidFilter(SPLINE_L2,spline,g,&ng,h,&nh,&IsCentered)==ERROR){
		printf("Unable to load the filter coeffiients");
		return;
	}

	/* Reducing */
	if (mxGetNumberOfDimensions(prhs[0]) == 2)
	{
		Reduce_2D(data1, m, n, data2, g, ng, IsCentered);
	}
	else
	{
		Reduce_3D(data1, m, n, l, data2, g, ng, IsCentered);
	}
}
