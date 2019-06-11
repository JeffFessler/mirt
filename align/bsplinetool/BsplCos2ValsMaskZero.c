/*******************************************************************************

        This file is a mex file to wrap C functions of B-spline 
	interpolation using zero end condition.

        Coded by Se Young Chun, Mar 1, 2008, University of Michigan

*******************************************************************************/

#include "mex.h"
#include "matrix.h"
#include "batchinterpol3D.h"

void mexFunction(int nlhs, mxArray *plhs[], 
    int nrhs, const mxArray *prhs[])
{
	int dim, m, n, l, ntot;
	mwSize dims[3];
 	float *data1, *data2, *data20, *data21, *data22, *data23,*wx, *wy, *wz;
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
		mexPrintf("[Gxx Gyx Gzx Gxy Gyy Gzy Gxz Gyz Gzz] = "
		    "BsplCos2ValsMaskZero(COEFFx, COEFFy, COEFFz, "
		    "[nx ny nz], [offx offy offz], [mx my mz], MASK, deg);\n");
		mexPrintf("[in]\n");
		mexPrintf("\tCOEFFx COEFFy COEFFz: bspline coefficients - "
			"single type\n");
		mexPrintf("\t[nx ny nz]          : output dimension\n");
		mexPrintf("\t[offx offy offz]    : offset of the origin\n");
		mexPrintf("\t[mx my mz]          : magnification factor for "
			"each direction\n");
		mexPrintf("\tMASK                : MASK, non zeros");
		mexPrintf("\tdeg           	: (opt) spline basis degree "
			"{default: 3}\n");
		mexPrintf("[out]\n");
		mexPrintf("\tGx Gy Gz, x y z		: gradient values\n");
		return;
	}

    	if (mxIsSingle(prhs[0]) != 1)
        {
                mexErrMsgTxt("First argument must be of single type\n");
                return;
        }

    	/* Retrieve the input data */
    	data1 = mxGetData(prhs[0]);

	if (mxGetNumberOfDimensions(prhs[0]) == 2)
	{
		m = mxGetDimensions(prhs[0])[0];
		n = mxGetDimensions(prhs[0])[1];
		l = 1;
	}
	else if (mxGetNumberOfDimensions(prhs[0]) == 3)
	{
		m = mxGetDimensions(prhs[0])[0];
		n = mxGetDimensions(prhs[0])[1];
		l = mxGetDimensions(prhs[0])[2];
	}
	else 
	{
                mexErrMsgTxt("First argument must be 2D or 3D\n");
		return;
	}

	ntot = mxGetNumberOfElements(prhs[1]);
	dim = ntot;
	if (ntot == 2)
	{
		nx = *((double*) mxGetData(prhs[1]));
		ny = *((double*) mxGetData(prhs[1]) + 1);
		dims[0] = (int)nx;
		dims[1] = (int)ny;
	}
	else if (ntot == 3)
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
                mexErrMsgTxt("Second argument should be either [nx ny nz] or "
			"[nx ny]\n");
		return;
	}

	ntot = mxGetNumberOfElements(prhs[2]);
	if (ntot == 2)
	{
		offx = *((double*) mxGetData(prhs[2]));
		offy = *((double*) mxGetData(prhs[2]) + 1);
	}
	else if (ntot == 3)
	{
		offx = *((double*) mxGetData(prhs[2]));
		offy = *((double*) mxGetData(prhs[2]) + 1);
		offz = *((double*) mxGetData(prhs[2]) + 2);
	}
	else
	{
		mexErrMsgTxt("Third argument should be either [offx offy offz]"
			" or [offx offy]\n");
		return;
	}

	ntot = mxGetNumberOfElements(prhs[3]);
	if (ntot == 2)
	{
		mx = *((double*) mxGetData(prhs[3]));
		my = *((double*) mxGetData(prhs[3]) + 1);
	}
	else if (ntot == 3)
	{
		mx = *((double*) mxGetData(prhs[3]));
		my = *((double*) mxGetData(prhs[3]) + 1);
		mz = *((double*) mxGetData(prhs[3]) + 2);
	}
	else
	{
       		mexErrMsgTxt("Fourth argument should be either [mx my mz] or"
			" [mx my]\n");
       		return;
	}

        ntot = mxGetNumberOfElements(prhs[4]);
	if (ntot == 0)
	{
	}
        else if (ntot == 2)
        {
                wx = (float*) mxGetData(mxGetCell(prhs[4], 0));
                wy = (float*) mxGetData(mxGetCell(prhs[4], 1));
        }
        else if (ntot == 3)
        {
                wx = (float*) mxGetData(mxGetCell(prhs[4], 0));
                wy = (float*) mxGetData(mxGetCell(prhs[4], 1));
                wz = (float*) mxGetData(mxGetCell(prhs[4], 2));
        }
        else
        {
           	mexErrMsgTxt("Fifth argument should be {}, {wx wy wz} or "
			"{wx wy}\n");
           	return;
        }


    	/* Create an mxArray for the output data */
	if (dim == 3)
	{
		if (nlhs == 4)
		{
			plhs[0] = mxCreateNumericArray(
					mxGetNumberOfElements(prhs[1]), dims, 
					mxSINGLE_CLASS, mxREAL);
			plhs[1] = mxCreateNumericArray(
					mxGetNumberOfElements(prhs[1]), dims, 
					mxSINGLE_CLASS, mxREAL);
			plhs[2] = mxCreateNumericArray(
					mxGetNumberOfElements(prhs[1]), dims, 
					mxSINGLE_CLASS, mxREAL);
			plhs[3] = mxCreateNumericArray(
					mxGetNumberOfElements(prhs[1]), dims, 
					mxSINGLE_CLASS, mxREAL);

			data20 = mxGetData(plhs[0]);
			data21 = mxGetData(plhs[1]);
			data22 = mxGetData(plhs[2]);
			data23 = mxGetData(plhs[3]);
			if (ntot == 0) /* no warp */
			{
				BatchInterpolatedAll3DZero (
					data20, data21, data22, data23, data1,
					m, n, l, nx, offx, mx, ny, offy, my, 
					nz, offz, mz, (long)spline); 
			}
			else
			{
				BatchInterpolatedAll3DZeroWarp (
					data20, data21, data22, data23, data1,
					m, n, l, nx, offx, mx, wx, ny, offy, my,
 					wy, nz, offz, mz, wz, (long)spline); 
			}
		}
		else
		{
			plhs[0] = mxCreateNumericArray(
					mxGetNumberOfElements(prhs[1]), dims, 
					mxSINGLE_CLASS, mxREAL);
    	
			data2 = mxGetData(plhs[0]);
			if (ntot == 0) /* no warp */
			{
				BatchInterpolatedValue3DZero(data2, data1, 
					m, n, l, nx, offx, mx, ny, offy, my, 
					nz, offz, mz, (long)spline); 
			}
			else
			{
				BatchInterpolatedValue3DZeroWarp(data2, data1, 
					m, n, l, nx, offx, mx, wx, ny, offy, my,
 					wy, nz, offz, mz, wz, (long)spline); 
			}
		}
	}
	else
	{
		if (nlhs == 3)
		{
			plhs[0] = mxCreateNumericArray(
					mxGetNumberOfElements(prhs[1]), dims, 
					mxSINGLE_CLASS, mxREAL);
			plhs[1] = mxCreateNumericArray(
					mxGetNumberOfElements(prhs[1]), dims, 
					mxSINGLE_CLASS, mxREAL);
			plhs[2] = mxCreateNumericArray(
					mxGetNumberOfElements(prhs[1]), dims, 
					mxSINGLE_CLASS, mxREAL);

			data20 = mxGetData(plhs[0]);
			data21 = mxGetData(plhs[1]);
			data22 = mxGetData(plhs[2]);
			if (ntot == 0) /* no warp */
			{
				BatchInterpolatedValue2DZero(data20,data1, m, n,
					 nx, offx, mx, ny, offy, my, 
					(long)spline); 
				BatchInterpolatedGradX2DZero(data21,data1, m, n,
					 nx, offx, mx, ny, offy, my, 
					(long)spline); 
				BatchInterpolatedGradY2DZero(data22,data1, m, n,
					 nx, offx, mx, ny, offy, my, 
					(long)spline); 
			}
			else
			{
				BatchInterpolatedValue2DZeroWarp(data20, data1, 
					m, n, nx, offx, mx, wx, ny, offy, my,wy,
					(long)spline); 
				BatchInterpolatedGradX2DZeroWarp(data21, data1, 
					m, n, nx, offx, mx, wx, ny, offy, my,wy,
					(long)spline); 
				BatchInterpolatedGradY2DZeroWarp(data22, data1, 
					m, n, nx, offx, mx, wx, ny, offy, my,wy,
					(long)spline); 
			}
		}
		else
		{
			plhs[0] = mxCreateNumericArray(
					mxGetNumberOfElements(prhs[1]), dims, 
					mxSINGLE_CLASS, mxREAL);

			data2 = mxGetData(plhs[0]);
			if (ntot == 0) /* no warp */
			{
				BatchInterpolatedValue2DZero(data2, data1, m, n,
					 nx, offx, mx, ny, offy, my, 
					(long)spline); 
			}
			else
			{
				BatchInterpolatedValue2DZeroWarp(data2, data1, 
					m, n, nx, offx, mx, wx, ny, offy, my,wy,
					(long)spline); 
			}
		}
	}
}




