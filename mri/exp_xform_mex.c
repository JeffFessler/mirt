// exp_xform_mex.c
// Copyright 2004-9-23, Jeff Fessler, University of Michigan

#include "mex.h"
#include <math.h>
#include <stdio.h>

#define Usage \
"Usage: y = function(x, u, v);\n\
in:\n\
	x	[N L]	complex vector(s)\n\
	u	[D N]	complex vectors\n\
	v	[D M]	complex vectors\n\
out:\n\
	y	[M L]	complex vector\n\
	y(m,l) = sum_n x(n,l) exp(-sum(u(:,n) .* v(:,m)))\n\
	This is the 'slow' 'exact' transform model for MRI.\n\
\n\
	All 3 inputs must be the same type (single or double).\n\
	Output type will be same as input.\n\
\n"

#ifndef Fail
#	define Fail(msg)	{ \
	(void)fprintf(stderr, "FAIL %s %d: %s\n", __FILE__, __LINE__, msg); \
	(void)fflush(stderr); \
	return 0; }
#endif


#define dtype float
static void exp_xform_float
#include "exp_xform_template.c"
#undef dtype

#define dtype double
static void exp_xform_double
#include "exp_xform_template.c"
#undef dtype


// exp_xform_mex()
// secondary gateway
static bool exp_xform_mex(
mxArray *plhs[],
const mxArray *mxx,
const mxArray *mxu,
const mxArray *mxv)
{
	const int is_single = mxIsSingle(mxx);

	if (is_single ^ mxIsSingle(mxu))	Fail("x & u type differ")
	if (is_single ^ mxIsSingle(mxv))	Fail("x & v type differ")

	if (!mxIsComplex(mxx))	Fail("x must be complex")
	if (!mxIsComplex(mxu))	Fail("u must be complex")
	if (!mxIsComplex(mxv))	Fail("v must be complex")

	// input sizes
	const int NN = mxGetM(mxx); // N
	const int LL = mxGetN(mxx); // L

	// check input size
	if (NN != (int) mxGetN(mxu))
	{
		printf("NN=%d size(u,2)=%d\n", NN, (int) mxGetN(mxu));
		Fail("ncol of 2nd input != 1st dim of 1st input")
	}
	const int DD = mxGetM(mxu);

	const int MM = mxGetN(mxv);
	if (DD != (int) mxGetM(mxv))
	{
		printf("DD=%d size(v,1)=%d\n", DD, (int) mxGetM(mxv));
		Fail("2nd & 3rd inputs must have same # of rows!")
	}

#if 0
	printf("MM=%d LL=%d NN=%d DD=%d\n", MM, LL, NN, DD);
#endif

	// create output matrix [M L]
	if (!(plhs[0] = mxCreateNumericMatrix (MM, LL,
		is_single ? mxSINGLE_CLASS : mxDOUBLE_CLASS,
		mxCOMPLEX)))
		Fail("create output failed")

	if (is_single)
	{
#define dtype float
		const dtype *xr = (const dtype *) mxGetData(mxx);
		const dtype *xi = (const dtype *) mxGetImagData(mxx);
		const dtype *ur = (const dtype *) mxGetData(mxu);
		const dtype *ui = (const dtype *) mxGetImagData(mxu);
		const dtype *vr = (const dtype *) mxGetData(mxv);
		const dtype *vi = (const dtype *) mxGetImagData(mxv);
		dtype *yr = (dtype *) mxGetData(plhs[0]);
		dtype *yi = (dtype *) mxGetImagData(plhs[0]);
#undef dtype

		for (int ll=0; ll < LL; ++ll)
			exp_xform_float(yr+ll*MM, yi+ll*MM, xr+ll*NN, xi+ll*NN,
					ur, ui, vr, vi, DD, MM, NN);
	}

	else
	{
#define dtype double
		const dtype *xr = (const dtype *) mxGetData(mxx);
		const dtype *xi = (const dtype *) mxGetImagData(mxx);
		const dtype *ur = (const dtype *) mxGetData(mxu);
		const dtype *ui = (const dtype *) mxGetImagData(mxu);
		const dtype *vr = (const dtype *) mxGetData(mxv);
		const dtype *vi = (const dtype *) mxGetImagData(mxv);
		dtype *yr = (dtype *) mxGetData(plhs[0]);
		dtype *yi = (dtype *) mxGetImagData(plhs[0]);
#undef dtype

		for (int ll=0; ll < LL; ++ll)
			exp_xform_double(yr+ll*MM, yi+ll*MM, xr+ll*NN, xi+ll*NN,
					ur, ui, vr, vi, DD, MM, NN);
	}

	return 1;
}


// Gateway routine - Interface with Matlab
void mexFunction(
int nlhs, mxArray *plhs[],
int nrhs, const mxArray *prhs[])
{
	if (!nlhs && !nrhs)
	{
		printf(Usage);
		return;
	}

	if (nlhs == 0 && nrhs == 1 && mxIsChar(prhs[0])) // check
		return;
	if (nlhs > 1 || nrhs != 3)
		fprintf(stderr, Usage);
	if (nlhs > 1)
		mexErrMsgTxt("One output argument required.");
	if (nrhs != 3)
		mexErrMsgTxt("Three input arguments required");
	if (!exp_xform_mex(plhs, prhs[0], prhs[1], prhs[2]))
		mexErrMsgTxt("exp_xform_mex() failed");
}

#undef Usage
