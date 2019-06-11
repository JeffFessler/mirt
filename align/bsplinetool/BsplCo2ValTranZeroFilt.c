/*******************************************************************************

        This file is a mex file to wrap C functions of B-spline
        transpose interpolation using filtering technique with
	zero end condition.

        Coded by Se Young Chun, Oct 30, 2007, the University of Michigan

*******************************************************************************/

#include "mex.h"
#include "matrix.h"

#include "filtinterpol2Dtran.h"
#include "filtinterpol3Dtran.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	int m, n, l, ntot;
	mwSize dims[3];
 	double *data1, *data2, *kernelx, *kernely, *kernelz;
	int nx, ny, nz, offx, offy, offz, mx, my, mz, nkx, nky, nkz;

	if (nlhs == 0 && nrhs == 1 && mxIsChar(prhs[0])) // check
		return;

	if (nlhs != 1)
	{
		mexPrintf("need nlhs=1\n");
		return;
	}

	if (nrhs == 5)
	{
		ntot = mxGetNumberOfElements(prhs[4]);
        	if (ntot == 2)
        	{
                	kernelx = (double*) mxGetData(mxGetCell(prhs[4], 0));
                	kernely = (double*) mxGetData(mxGetCell(prhs[4], 1));
			kernelz = NULL;
			nkx = (int)(mxGetDimensions(mxGetCell(prhs[4], 0))[0] 
				* mxGetDimensions(mxGetCell(prhs[4], 0))[1]);
			nky = (int)(mxGetDimensions(mxGetCell(prhs[4], 1))[0] 
				* mxGetDimensions(mxGetCell(prhs[4], 1))[1]);
			nkz = 0;
        	}
        	else if (ntot == 3)
        	{
                	kernelx = (double*) mxGetData(mxGetCell(prhs[4], 0));
                	kernely = (double*) mxGetData(mxGetCell(prhs[4], 1));
                	kernelz = (double*) mxGetData(mxGetCell(prhs[4], 2));
			nkx = (int)(mxGetDimensions(mxGetCell(prhs[4], 0))[0] 
				* mxGetDimensions(mxGetCell(prhs[4], 0))[1]);
			nky = (int)(mxGetDimensions(mxGetCell(prhs[4], 1))[0] 
				* mxGetDimensions(mxGetCell(prhs[4], 1))[1]);
			nkz = (int)(mxGetDimensions(mxGetCell(prhs[4], 2))[0] 
				* mxGetDimensions(mxGetCell(prhs[4], 2))[1]);
        	}
        	else
        	{
                	mexErrMsgTxt("Fifth argument should be {KERNELx KERNELy"
				" (KERNELz)}\n");
                	return;
        	}
	}
	else
	{
                mexPrintf("IMG = BsplCo2ValTranZeroFilt(COEFF, [nx ny (nz)], "
			"[offx offy (offz)], [mx my (mz)], {KERNELx KERNELy "
			"(KERNELz)});\n");
                mexPrintf("[in]\n");
                mexPrintf("\tCOEFF              : bspline coefficients - single"
			" type\n");
                mexPrintf("\t[nx ny (nz)]       : output dimension\n");
                mexPrintf("\t[offx offy (offz)] : offset of the origin\n");
                mexPrintf("\t[mx my (mz)]       : magnification factor for each"
			" direction\n");
                mexPrintf("\tKERNEL             : FIR kernel "
			"filter coefficients, bm^n\n");
                mexPrintf("[out]\n");
                mexPrintf("\tIMG                : interpolated values\n");
		return;
	}

    	if (mxIsDouble(prhs[0]) != 1)
        {
                mexErrMsgTxt("First argument must be of double type\n");
                return;
        }

    	/* Retrieve the input data */
    	data1 = mxGetData(prhs[0]);

	if (mxGetNumberOfDimensions(prhs[0]) == 2)
	{
		m = (int)mxGetDimensions(prhs[0])[0];
		n = (int)mxGetDimensions(prhs[0])[1];
		l = 1;
	}
	else if (mxGetNumberOfDimensions(prhs[0]) == 3)
	{
		m = (int)mxGetDimensions(prhs[0])[0];
		n = (int)mxGetDimensions(prhs[0])[1];
		l = (int)mxGetDimensions(prhs[0])[2];
	}
	else 
	{
                mexErrMsgTxt("First argument must be 2D or 3D\n");
		return;
	}

	ntot = mxGetNumberOfElements(prhs[1]);
	if (ntot == 2)
	{
		nx = (int)*((double*) mxGetData(prhs[1]));
		ny = (int)*((double*) mxGetData(prhs[1]) + 1);
		nz = 1;
		dims[0] = nx;
		dims[1] = ny;
	}
	else if (ntot == 3)
	{
		nx = (int)*((double*) mxGetData(prhs[1]));
		ny = (int)*((double*) mxGetData(prhs[1]) + 1);
		nz = (int)*((double*) mxGetData(prhs[1]) + 2);
		dims[0] = nx;
		dims[1] = ny;
		dims[2] = nz;
	}
	else
	{
                mexErrMsgTxt("Second argument should be either [nx ny nz] or"
			" [nx ny]\n");
		return;
	}

	ntot = mxGetNumberOfElements(prhs[2]);
	if (ntot == 2)
	{
		offx = (int)*((double*) mxGetData(prhs[2]));
		offy = (int)*((double*) mxGetData(prhs[2]) + 1);
		offz = 0;
	}
	else if (ntot == 3)
	{
		offx = (int)*((double*) mxGetData(prhs[2]));
		offy = (int)*((double*) mxGetData(prhs[2]) + 1);
		offz = (int)*((double*) mxGetData(prhs[2]) + 2);
	}
	else
	{
                mexErrMsgTxt("Third argument should be either [offx offy"
			" offz] or [offx offy]\n");
		return;
	}

	ntot = mxGetNumberOfElements(prhs[3]);
	if (ntot == 2)
	{
		mx = (int)*((double*) mxGetData(prhs[3]));
		my = (int)*((double*) mxGetData(prhs[3]) + 1);
		mz = 0;
	}
	else if (ntot == 3)
	{
		mx = (int)*((double*) mxGetData(prhs[3]));
		my = (int)*((double*) mxGetData(prhs[3]) + 1);
		mz = (int)*((double*) mxGetData(prhs[3]) + 2);
	}
	else
	{
                mexErrMsgTxt("Fourth argument should be either [mx my mz]"
			" or [mx my]\n");
		return;
	}

    	/* Create an mxArray for the output data */
	plhs[0] = mxCreateNumericArray(mxGetNumberOfElements(prhs[1]), dims, 
			mxDOUBLE_CLASS, mxREAL);

	/* Retrieve the output data */
    	data2 = mxGetData(plhs[0]);

	if (ntot == 3)
	{
		FiltInterpolatedValue3DTranZero(data2, data1, m, n, l, nx, offx,
 			mx, ny, offy, my, nz, offz, mz, 
			kernelx, nkx, kernely, nky, kernelz, nkz); 
	}
	else
	{
		FiltInterpolatedValue2DTranZero(data2, data1, m, n, nx, offx, 
			mx, ny, offy, my, kernelx, nkx, kernely, nky); 
	}
}
