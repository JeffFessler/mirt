/*******************************************************************************

        This file is a mex file to wrap C functions of B-spline
        gradient interpolation using mirror condition.

        Coded by Se Young Chun, Oct 30, 2007, the University of Michigan

*******************************************************************************/

#include "mex.h"
#include "matrix.h"

#include "batchinterpol3D.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	int dim, m, n, l, ntot;
	mwSize dims[3];
 	float *data1, *data2, *wx, *wy, *wz;
	double spline, nx, ny, nz, offx, offy, offz, mx, my, mz;

	if (nlhs == 0 && nrhs == 1 && mxIsChar(prhs[0])) // check
		return;

	if (nlhs != 1)
	{
		mexPrintf("need nlhs=1\n");
		return;
	}

	if (nrhs == 6)
	{
		spline = mxGetScalar(prhs[5]);
	}
	else if (nrhs == 5)
	{
		spline = 3;
	}
	else
	{
                mexPrintf("VAL = BsplCo2GdZMirr(COEFF, [nx ny nz], [offx offy"
			" offz], [mx my mz], {wx wy wz}, deg);\n");
                mexPrintf("[in]\n");
                mexPrintf("\tCOEFF              : bspline coefficients -"
			" single type\n");
                mexPrintf("\t[nx ny nz]         : output dimension\n");
                mexPrintf("\t[offx offy offz]   : offset of the origin\n");
                mexPrintf("\t[mx my mz]         : magnification factor for"
			" each direction\n");
                mexPrintf("\t{wx wy wz}         : deformed value vectors\n");
                mexPrintf("\tdeg                : (opt) spline basis degree"
			" {default: 3}\n");
                mexPrintf("[out]\n");
                mexPrintf("\tVAL                : interpolated grad z"
			" values\n");
		return;
	}

    	if (mxIsSingle(prhs[0]) != 1)
        {
                mexErrMsgTxt("First argument must be of single type\n");
                return;
        }

    	/* Retrieve the input data */
    	data1 = mxGetData(prhs[0]);

	if (mxGetNumberOfDimensions(prhs[0]) == 3)
	{
		m = mxGetDimensions(prhs[0])[0];
		n = mxGetDimensions(prhs[0])[1];
		l = mxGetDimensions(prhs[0])[2];
	}
	else 
	{
                mexErrMsgTxt("First argument must be 3D\n");
		return;
	}

	ntot = mxGetNumberOfElements(prhs[1]);
	dim = ntot;
	if (ntot == 3)
	{
		nx = *((double*) mxGetData(prhs[1]));
		ny = *((double*) mxGetData(prhs[1]) + 1);
		nz = *((double*) mxGetData(prhs[1]) + 2);
		dims[0] = (int)nx;
		dims[1] = (int)ny;
		dims[2] = (int)nz;
	}
	else
	{
                mexErrMsgTxt("Second argument should be [nx ny nz]\n");
		return;
	}

	ntot = mxGetNumberOfElements(prhs[2]);
	if (ntot == 3)
	{
		offx = *((double*) mxGetData(prhs[2]));
		offy = *((double*) mxGetData(prhs[2]) + 1);
		offz = *((double*) mxGetData(prhs[2]) + 2);
	}
	else
	{
                mexErrMsgTxt("Third argument should be [offx offy offz]\n");
		return;
	}

	ntot = mxGetNumberOfElements(prhs[3]);
	if (ntot == 3)
	{
		mx = *((double*) mxGetData(prhs[3]));
		my = *((double*) mxGetData(prhs[3]) + 1);
		mz = *((double*) mxGetData(prhs[3]) + 2);
	}
	else
	{
                mexErrMsgTxt("Fourth argument should be [mx my mz]\n");
		return;
	}

	ntot = mxGetNumberOfElements(prhs[4]);
	if (ntot == 0)
	{
	}
	else if (ntot == 3)
	{
		wx = (float*) mxGetData(mxGetCell(prhs[4], 0));
		wy = (float*) mxGetData(mxGetCell(prhs[4], 1));
		wz = (float*) mxGetData(mxGetCell(prhs[4], 2));
	}
	else
	{
                mexErrMsgTxt("Fifth argument should be either {} or"
			" {wx wy wz}\n");
		return;
	}

    	/* Create an mxArray for the output data */
	plhs[0] = mxCreateNumericArray(mxGetNumberOfElements(prhs[1]), dims, 
			mxSINGLE_CLASS, mxREAL);

	/* Retrieve the output data */
    	data2 = mxGetData(plhs[0]);

	if (ntot == 0)
	{
		BatchInterpolatedGradZ3DMirr(data2, data1, m, n, l, nx, offx, 
				mx, ny, offy, my, nz, offz, mz, (long)spline); 
	}
	else
	{
		BatchInterpolatedGradZ3DMirrWarp(data2, data1, m, n, l, nx, 
				offx, mx, wx, ny, offy, my, wy, nz, offz, mz, 
				wz, (long)spline); 
	}
}

