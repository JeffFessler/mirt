// interp4_table_adj_mex.c
// see interp4_table1_adj.c
// Mex file for *adjoint* of 2D periodic interpolation using table lookup.
// ck = interp4_table_adj_mex(fm, h1_table, h2_table J, L, tm, K)
//
// Copyright 2004-4-2 Jeff Fessler and Yingying Zhang, University of Michigan
// Extended 2013-08-02 to 4D by David Johnson, The Ohio State University Wexner Medical Center


#include "def,table.h"
#include "def,table4.h"


// interp4_table_adj_mex()
// Usage: ck = function(fm, h1, h2, h3, h4, J, L, tm, K)
static int interp4_table_adj_mex(
mxArray *plhs[],
const mxArray *mx_fm,
const mxArray *mx_h1,
const mxArray *mx_h2,
const mxArray *mx_h3,
const mxArray *mx_h4,
const mxArray *mx_J,
const mxArray *mx_L,
const mxArray *mx_tm,
const mxArray *mx_K,
const mxArray *mx_order, // optional: may be NULL
const mxArray *mx_flips) // optional: may be NULL
{
	const int M = mxGetM(mx_fm); // # of time samples
	const int N = mxGetN(mx_fm); // # of realizations, prod(size(2:end))

	const int *Jd = (const int *) mxGetData(mx_J);
	const int *Ld = (const int *) mxGetData(mx_L);
	const int *Kd = (const int *) mxGetData(mx_K);

	const double *p_tm = mxGetPr(mx_tm);
	const double *r_fm = mxGetPr(mx_fm);
	const double *i_fm = mxGetPi(mx_fm);
	const double *r_h1 = mxGetPr(mx_h1);
	const double *r_h2 = mxGetPr(mx_h2);
	const double *r_h3 = mxGetPr(mx_h3);
	const double *r_h4 = mxGetPr(mx_h4);

	const int order = mx_order ? *((int *) mxGetData(mx_order)) : 0;
	const int *flips = mx_flips ? ((int *) mxGetData(mx_flips)) : NULL;

	Call(mxIsComplexDouble, (mx_fm))
	Call(mxIsRealDouble, (mx_tm))

	if (!mxIsInt32n(mx_J, 4))
		Fail("J must be int32 [1,4]")
	if (!mxIsInt32n(mx_L, 4))
		Fail("L must be int32 [1,4]")
	if (!mxIsInt32n(mx_K, 4))
		Fail("K must be int32 [1,4]")
	if (mx_order && !mxIsScalarInt32(mx_order))
		Fail("order must be scalar int32 (0 | 1)")
	if (mx_flips && !mxIsInt32n(mx_flips, 4))
		Fail("flips must be [1 4] int32 (0 | 1)")

	/* check h1,h2,h3,h4 tables' sizes */
	if ((int) mxGetM(mx_h1) != Jd[0]*Ld[0]+1 || (mxGetN(mx_h1) != 1))
	{
		fprintf(stderr, "J1=%d L1=%d tablelength=%d\n",
			Jd[0], Ld[0], (int) mxGetM(mx_h1));
		Fail("h1 size problem")
	}

	if ((int) mxGetM(mx_h2) != Jd[1]*Ld[1]+1 || (mxGetN(mx_h2) != 1))
	{
		fprintf(stderr, "J2=%d L2=%d tablelength=%d\n",
			Jd[1], Ld[1], (int) mxGetM(mx_h2));
		Fail("h2 size problem")
	}

	if ((int) mxGetM(mx_h3) != Jd[2]*Ld[2]+1 || (mxGetN(mx_h3) != 1))
	{
		fprintf(stderr, "J3=%d L3=%d tablelength=%d\n",
			Jd[2], Ld[2], (int) mxGetM(mx_h3));
		Fail("h3 size problem")
	}

	if ((int) mxGetM(mx_h4) != Jd[3]*Ld[3]+1 || (mxGetN(mx_h4) != 1))
	{
		fprintf(stderr, "J4=%d L4=%d tablelength=%d\n",
			Jd[3], Ld[3], (int) mxGetM(mx_h4));
		Fail("h4 size problem")
	}

	if (M != (int) mxGetM(mx_tm) || 4 != mxGetN(mx_tm))
		Fail("t_m must be Mx4 matrix")

	// create a new array and set the output pointer to it
	Call(plhs[0] = mxCreateDoubleMatrix, (Kd[0]*Kd[1]*Kd[2]*Kd[3], N, mxCOMPLEX))
	double *r_ck = mxGetPr(plhs[0]);
	double *i_ck = mxGetPi(plhs[0]);

	// call the C subroutine N times; once for each realization
	if (	mxIsComplexDouble(mx_h1) &&
		mxIsComplexDouble(mx_h2) &&
		mxIsComplexDouble(mx_h3) &&
		mxIsComplexDouble(mx_h4))
	{
		const double *i_h1 = mxGetPi(mx_h1);
		const double *i_h2 = mxGetPi(mx_h2);
		const double *i_h3 = mxGetPi(mx_h3);
		const double *i_h4 = mxGetPi(mx_h4);
		if (order)
			Fail("only 0th order implemented for complex")

		for (int nn=0; nn < N; ++nn)
		{
			interp4_table0_complex_per_adj(r_ck, i_ck,
				Kd[0], Kd[1], Kd[2], Kd[3],
				r_h1, i_h1, r_h2, i_h2, r_h3, i_h3, r_h4, i_h4,
				Jd[0], Jd[1], Jd[2], Jd[3], Ld[0], Ld[1], Ld[2], Ld[3],
				p_tm, M, r_fm, i_fm);
			r_ck += Kd[0]*Kd[1]*Kd[2]*Kd[3]; i_ck += Kd[0]*Kd[1]*Kd[2]*Kd[3];
			r_fm += M; i_fm += M;
		}
	}

	else if (mxIsRealDouble(mx_h1) &&
		mxIsRealDouble(mx_h2) &&
		mxIsRealDouble(mx_h3) &&
		mxIsRealDouble(mx_h4))
	{
		interp4_table_real_per_adj_t *fun;
		if (order == 0)		fun = interp4_table0_real_per_adj;
		else if (order == 1)	fun = interp4_table1_real_per_adj;
		else Fail("bad order")

#ifndef Provide_flip
		if (flips && (flips[0] || flips[1] || flips[2] || flips[3]))
			Fail("flip not compiled")
#endif

		for (int nn=0; nn < N; ++nn)
		{
			fun(r_ck, i_ck, Kd[0], Kd[1], Kd[2], Kd[3], r_h1, r_h2, r_h3, r_h4,
#ifdef Provide_flip
				flips ? flips[0] : 0,
				flips ? flips[1] : 0,
				flips ? flips[2] : 0,
#endif
				Jd[0], Jd[1], Jd[2], Jd[3], Ld[0], Ld[1], Ld[2], Ld[3],
				p_tm, M, r_fm, i_fm);
			r_ck += Kd[0]*Kd[1]*Kd[2]*Kd[3]; i_ck += Kd[0]*Kd[1]*Kd[2]*Kd[3];
			r_fm += M; i_fm += M;
		}
	}

	else
		Fail("h must be real or complex double (preferably real)")

	return 1;
}


// gateway
// ck = function(fm, h1_table, h2_table, h3_table, h4_table, J, L, tm, K)
void mexFunction(
const int nlhs, mxArray *plhs[],
const int nrhs, const mxArray *prhs[])
{
	// check for the proper number of arguments

	if (nlhs == 0 && nrhs == 1 && mxIsChar(prhs[0])) // check
		return;

	if (nrhs < 9 || nrhs > 11)
		mexFail("11 inputs needed: (f, h1, h2, h3, h4, J, L, t, K, [order, flips])")
	if (nlhs > 1)
		mexFail("Less than one output arguments.")

	if (!interp4_table_adj_mex(plhs,
		prhs[0], prhs[1], prhs[2], prhs[3], prhs[4],
		prhs[5], prhs[6], prhs[7], prhs[8],
		nrhs > 9 ? prhs[9] : NULL,
		nrhs > 10 ? prhs[10] : NULL))
		mexFail("interp4_table_adj_mex() failed")
}
