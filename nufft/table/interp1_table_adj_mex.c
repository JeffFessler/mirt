// interp1_table_adj_mex.c
// Mex file for *adjoint* of 1D periodic interpolation using table lookup.
//
// forward direction: (for m = 1,...,M)
// f(t_m) = \sum_{k=0}^{K-1} c_k h( (t_m - k) mod K )
//
// adjoint direction: (for k=0,...,K-1) (note complex conjugate!)
// c_k = \sum_{m=1}^M f(t_m) h^*( (t_m - k) mod K )
//
// The interpolator h is nonzero (and tabulated) for -J/2 <= t <= J/2.
//
// Copyright 2004-4-1 Jeff Fessler and Yingying Zhang, University of Michigan

#include "def,table.h"
#include "def,table1.h"


// interp1_table_adj_mex()
// usage: ck = function(fm, h_table, J, L, tm, K, order, flips)
static int interp1_table_adj_mex(
mxArray *plhs[],
const mxArray *mx_fm,
const mxArray *mx_h1,
const mxArray *mx_J,
const mxArray *mx_L,
const mxArray *mx_tm,
const mxArray *mx_K,
const mxArray *mx_order, // optional: may be NULL
const mxArray *mx_flips) // optional: may be NULL
{
	const int M = mxGetM(mx_fm); // # of time samples
	const int N = mxGetN(mx_fm); // # of realizations, prod(size(2:end))

	const int J1 = *((const int *) mxGetData(mx_J));
	const int K1 = *((const int *) mxGetData(mx_K));
	const int L1 = *((const int *) mxGetData(mx_L));

	const double *p_tm = mxGetPr(mx_tm);
	const double *r_fm = mxGetPr(mx_fm);
	const double *i_fm = mxGetPi(mx_fm);
	const double *r_h1 = mxGetPr(mx_h1);

	const int order = mx_order ? *((const int *) mxGetData(mx_order)) : 0;
	const int flip1 = mx_flips ? *((const int *) mxGetData(mx_flips)) : 0;

//	if (N != 1)
//		fprintf(stderr, "Caution: multiple columns?");

	Call(mxIsComplexDouble, (mx_fm))
	Call(mxIsRealDouble, (mx_tm))

	// J L K must be scalar
	if (!mxIsScalarInt32(mx_J))
		Fail("J must be scalar int32")
	if (!mxIsScalarInt32(mx_K))
		Fail("K must be scalar int32")
	if (!mxIsScalarInt32(mx_L))
		Fail("L must be scalar int32")

	if (mx_order && !mxIsScalarInt32(mx_order))
		Fail("order must be scalar int32 (0 | 1)")
	if (mx_flips && !mxIsScalarInt32(mx_flips))
		Fail("flips must be scalar int32 (0 | 1)")

	// check h1 table size
	if ((int) mxGetM(mx_h1) != J1*L1+1)
	{
		fprintf(stderr, "J=%d L=%d tablelength=%d\n",
			J1, L1, (int) mxGetM(mx_h1));
		Fail("h size problem")
	}
	if (mxGetN(mx_h1) != 1)
		Fail("h must be col vector")

	if (M != (int) mxGetM(mx_tm) || 1 != mxGetN(mx_tm))
		Fail("t_m must be Mx1 col vector")

	// create a new array and set the output pointer to it
	Call(plhs[0] = mxCreateDoubleMatrix, (K1, N, mxCOMPLEX))
	double *r_ck = mxGetPr(plhs[0]);
	double *i_ck = mxGetPi(plhs[0]);

	// call the C subroutine N times; once for each realization
	if (mxIsComplexDouble(mx_h1))
	{
		const double *i_h1 = mxGetPi(mx_h1);
		for (int nn=0; nn < N; ++nn)
		{
			interp1_table0_complex_per_adj(r_ck, i_ck, K1,
				r_h1, i_h1, J1, L1, p_tm, M, r_fm, i_fm);
			r_ck += K1; i_ck += K1;
			r_fm += M; i_fm += M;
		}
	}

	else if (mxIsRealDouble(mx_h1))
	{
		interp1_table_real_per_adj_t *fun;
		if (order == 0)		fun = interp1_table0_real_per_adj;
		else if (order == 1)	fun = interp1_table1_real_per_adj;
		else Fail("bad order")

#ifndef Provide_flip
		if (flip1)
			Fail("flip not compiled")
#endif

		for (int nn=0; nn < N; ++nn)
		{
			fun(r_ck, i_ck, K1, r_h1,
#ifdef Provide_flip
				flip1,
#endif
				J1, L1, p_tm, M, r_fm, i_fm);
			r_ck += K1; i_ck += K1;
			r_fm += M; i_fm += M;
		}
	}

	else
		Fail("h must be real or complex double (preferably real)")

	return 1;
}


// gateway
// Usage: ck = function(fm, h_table, J, L, tm, K)
void mexFunction(
const int nlhs, mxArray *plhs[],
const int nrhs, const mxArray *prhs[])
{
	if (nlhs == 0 && nrhs == 1 && mxIsChar(prhs[0])) // check
		return;

	if (nrhs < 6 || nrhs > 8)
		mexFail("8 inputs needed: (f, h, J, L, t, K, [order, flips])")
	if (nlhs > 1)
		mexFail("Less than one output arguments.")

	if (!interp1_table_adj_mex(plhs,
		prhs[0], prhs[1], prhs[2], prhs[3], prhs[4], prhs[5],
		nrhs > 6 ? prhs[6] : NULL,
		nrhs > 7 ? prhs[7] : NULL))
		mexFail("interp1_table_adj_mex() failed")
}
