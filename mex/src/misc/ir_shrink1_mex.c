// ir_shrink1_mex.c
// Matlab mex gateway for iterative 1D shrinkage algorithm
//
// min_x 1/2 (y - x)^2 + reg * pot(x)
// user must provide tabulated dpot to be used with linear interpolation.
//
// Copyright 2012-06-25, Jeff Fessler, University of Michigan

#include "jf,mex,def.h"
#include "jf,thread1.h"
// #include <unistd.h>

// #if !defined(Need_ir_shrink1_mex_gateway)
#include "ir_shrink1.h"
// #endif

#define Usage "usage error. see above"

static void ir_shrink1_mex_help(void)
{
	printf("\n\
\n\
	min_x 1/2 (y - x)^2 + reg * pot(x)\n\
\n\
	x = function(y, reg, table, niter, thresh, nthread, chat)\n\
	y: (single) [(N)]\n\
	reg: (single) [(N)]\n\
	table: (single) [K] samples (1:K)*dt\n\
	dt: (single) table spacing\n\
	niter: (int32) # of iterations (max)\n\
	thresh: (single) change in x threshold for stopping rule\n\
	nthread: (int32) # of threads\n\
	chat: (int32) be verbose with threading\n\
\n");
}


// ir_shrink1_p_mex()
// x = function(y, reg, table, nthread, chat)
static sof ir_shrink1_p_mex(
mxArray *plhs[], // [(N)]
Cmx mx_y, // [(N)]
Cmx mx_reg, // [(N)]
Cmx mx_table, // [K]
Cmx mx_dt, // [1]
Cmx mx_niter, // [1]
Cmx mx_thresh, // [1]
Cmx mx_nthread, // [1]
Cmx mx_chat)
{
	Call(mxIsRealSingle, (mx_y))
	cint N = mxGetM(mx_y) * mxGetN(mx_y);

	Call(mxIsRealSingle, (mx_reg))
	if (N != (int) (mxGetM(mx_reg) * mxGetN(mx_reg)))
		Fail("y and reg must have same size")

	Call(mxIsRealSingle, (mx_table))
	cint K = mxGetM(mx_table) * mxGetN(mx_table);

	Call(mxIsScalarSingle, (mx_dt))
	cfloat dt = mxGetSingle(mx_dt);

	Call(mxIsScalarInt32, (mx_niter))
	cint niter = mxGetInt(mx_niter);

	Call(mxIsScalarSingle, (mx_thresh))
	cfloat thresh = mxGetSingle(mx_thresh);

	Call(mxIsScalarInt32, (mx_nthread))
	cint nthread = mxGetInt(mx_nthread);

	Call(mxIsScalarInt32, (mx_chat))
	cint chat = mxGetInt(mx_chat);

	// output
	Call(plhs[0] = mxCreateNumericMatrix,
		(mxGetM(mx_y), mxGetN(mx_y), mxSINGLE_CLASS, mxREAL))

	Call(ir_shrink1_p, ( (float *) mxGetData(plhs[0]),
			mxGetPr_cfloat(mx_y), N,
			mxGetPr_cfloat(mx_reg),
			mxGetPr_cfloat(mx_table), dt, K,
			thresh, niter, nthread, chat))
	Ok
}


// intermediate gateway routine
#if defined(Need_ir_shrink1_mex_gateway)
static
#endif
sof ir_shrink1_mex(
int nlhs, mxArray *plhs[],
int nrhs, Const mxArray *prhs[])
{
	// for "check" option needed by test_all_mex
	if (nrhs == 1 && mxIsChar(prhs[0]))
		Ok

	if (nlhs != 1 || nrhs < 8)
	{
		ir_shrink1_mex_help();
		Call(mxu_arg, (nrhs, prhs))
		Fail(Usage)
	}

	Call(ir_shrink1_p_mex, (plhs,
		prhs[0], prhs[1], prhs[2], prhs[3],
		prhs[4], prhs[5], prhs[6], prhs[7]))

	Ok
}


#if defined(Need_ir_shrink1_mex_gateway)
// gateway routine
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, Const mxArray *prhs[])
{
	if (!nlhs && !nrhs)
	{
		ir_shrink1_mex_help();
		return;
	}
	if (!ir_shrink1_mex(nlhs, plhs, nrhs, prhs))
		mexErrMsgTxt("ir_shrink1_mex");
}
#endif
