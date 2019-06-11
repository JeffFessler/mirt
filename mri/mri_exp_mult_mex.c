// mri_exp_mult_mex.c
// Sangwoo Lee, University of Michigan, Jun, 2004
// 2004-6-21 modified by JF

#include "mex.h"
#include <math.h>
#include <stdio.h>

#define Usage \
"Usage: D = function(A, u, v);\n\
in:\n\
	A	[N L]	complex matrix\n\
	u	[N 1]	vector\n\
	v	[M 1]	vector\n\
			one (and only one) of u and v must be complex!\n\
out:\n\
	D	[L M]	complex vector, D = A' * exp(-u * v.')\n\
			D_lm = sum_n A_nl^* exp(-u_n v_m)\n\
\n\
 This function is a memory efficient and fast implementation\n\
 of AA matrix computation in int_tim_seg.m function in NUFFT package.\n\
\n"

#ifndef Fail
#	define Fail(msg) { \
	(void)fprintf(stderr, "FAIL %s %d: %s\n", __FILE__, __LINE__, msg); \
	(void)fflush(stderr); \
	return 0; }
#endif

#define dtype double
static bool mri_exp_mult_double
#include "mri_exp_mult_template.c"
#undef dtype

#define dtype float
static bool mri_exp_mult_float
#include "mri_exp_mult_template.c"
#undef dtype


// mri_exp_mult_mex()
// secondary gateway
static bool mri_exp_mult_mex(
mxArray *plhs[],
const mxArray *mx_A,
const mxArray *mx_u,
const mxArray *mx_v)
{
	const int is_single = mxIsSingle(mx_A);
	double *work_r, *work_i;

	if (is_single ^ mxIsSingle(mx_u)) Fail("A & u type differ")
	if (is_single ^ mxIsSingle(mx_v)) Fail("A & v type differ")

	if (!mxIsComplex(mx_A))  Fail("A must be complex")
	if (!(mxIsComplex(mx_u) ^ mxIsComplex(mx_v)))
		Fail("need exactly one of u and v to be complex")

	// input sizes
	const int NN = mxGetM(mx_A); // N
	const int LL = mxGetN(mx_A); // L

	// check input size
	if (NN != (int) mxGetM(mx_u))
		Fail("size of 2nd input != 1st dim of 1st input")
	if (1 != mxGetN(mx_u))
		Fail("2nd input argument must be a column vector!")
	const int MM = mxGetM(mx_v);
	if (1 != mxGetN(mx_v))
		Fail("3rd input argument must be a column vector!")

	// create output matrix [L M]
	if (!(plhs[0] = mxCreateNumericMatrix (LL, MM,
		is_single ? mxSINGLE_CLASS : mxDOUBLE_CLASS,
		mxCOMPLEX)))
		Fail("allocate output")

	if (!(work_r = (double *) calloc(NN, sizeof(*work_r)))) Fail("alloc")
	if (!(work_i = (double *) calloc(NN, sizeof(*work_i)))) Fail("alloc")

	if (is_single)
	{
#define dtype float
		const dtype *Ar = (const dtype *) mxGetData(mx_A);
		const dtype *Ai = (const dtype *) mxGetImagData(mx_A);
		const dtype *ur = (const dtype *) mxGetData(mx_u);
		const dtype *ui = (const dtype *) mxGetImagData(mx_u);
		const dtype *vr = (const dtype *) mxGetData(mx_v);
		const dtype *vi = (const dtype *) mxGetImagData(mx_v);
		dtype *Dr = (dtype *) mxGetData(plhs[0]);
		dtype *Di = (dtype *) mxGetImagData(plhs[0]);
#undef dtype

		if (!mri_exp_mult_float(Dr, Di, work_r, work_i,
			Ar, Ai, ur, ui, vr, vi, LL, MM, NN))
			Fail("mri_exp_mult")
	}

	else
	{
#define dtype double
		const dtype *Ar = (const dtype *) mxGetData(mx_A);
		const dtype *Ai = (const dtype *) mxGetImagData(mx_A);
		const dtype *ur = (const dtype *) mxGetData(mx_u);
		const dtype *ui = (const dtype *) mxGetImagData(mx_u);
		const dtype *vr = (const dtype *) mxGetData(mx_v);
		const dtype *vi = (const dtype *) mxGetImagData(mx_v);
		dtype *Dr = (dtype *) mxGetData(plhs[0]);
		dtype *Di = (dtype *) mxGetImagData(plhs[0]);
#undef dtype

		if (!mri_exp_mult_double(Dr, Di, work_r, work_i,
			Ar, Ai, ur, ui, vr, vi, LL, MM, NN))
			Fail("mri_exp_mult")
	}

	free(work_r);
	free(work_i);

	return 1;
}


// gateway routine - interface with Matlab
void mexFunction(
int nlhs, mxArray *plhs[],
int nrhs, const mxArray *prhs[])
{
	if (!nlhs && !nrhs)
	{
#if 0 // defined(Date)
		printf("compiled %s\n", Date);
#endif
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
	if (!mri_exp_mult_mex(plhs, prhs[0], prhs[1], prhs[2]))
		mexErrMsgTxt("mri_exp_mult_mex() failed");
}

#undef Usage
