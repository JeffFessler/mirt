// interp2_table_mex.c
// Mex file for 2D periodic interpolation using tabulated interpolator.
// see interp2_table1_for.c
// Copyright 2004-03-30 Yingying Zhang and Jeff Fessler, University of Michigan

#include "def,table.h"
#include "def,table2.h"


// interp2_table_mex()
static int interp2_table_mex(
mxArray *plhs[],
const mxArray *mx_ck, // [K1 K2] DFT coefficients
const mxArray *mx_h1,
const mxArray *mx_h2,
const mxArray *mx_J,
const mxArray *mx_L,
const mxArray *mx_tm, // [M]
const mxArray *mx_order, // optional: may be NULL
const mxArray *mx_flips) // optional: may be NULL
{
	const int ndim = mxGetNumberOfDimensions(mx_ck);
	const int K1 = (ndim > 0) ? (mxGetDimensions(mx_ck))[0] : 1;
	const int K2 = (ndim > 1) ? (mxGetDimensions(mx_ck))[1] : 1;
	const int N = (ndim > 2) ? (mxGetDimensions(mx_ck))[2] : 1;
	const int M = mxGetM(mx_tm); // # of time samples

	const int *Jd = (const int *) mxGetData(mx_J);
	const int *Ld = (const int *) mxGetData(mx_L);

	const double *p_tm = mxGetPr(mx_tm);
	const double *r_h1 = mxGetPr(mx_h1);
	const double *r_h2 = mxGetPr(mx_h2);
	const double *r_ck = mxGetPr(mx_ck);
	const double *i_ck = mxGetPi(mx_ck);

	const int order = mx_order ? *((const int *) mxGetData(mx_order)) : 0;
	const int *flips = mx_flips ? ((const int *) mxGetData(mx_flips)) : NULL;

	Call(mxIsComplexDouble, (mx_ck))
	Call(mxIsRealDouble, (mx_tm))

	// J L must be [1,2]
	if (!mxIsInt32n(mx_J, 2))
		Fail("J must be [1,2]")
	if (!mxIsInt32n(mx_L, 2))
		Fail("L must be [1,2]")
	if (mx_order && !mxIsScalarInt32(mx_order))
		Fail("order must be scalar int32 (0 | 1)")
	if (mx_flips && !mxIsInt32n(mx_flips, 2))
                Fail("flips must be [1 2] int32 (0 | 1)")

	// check size, type of tables
	if ((int) mxGetM(mx_h1) != Jd[0]*Ld[0]+1 || (mxGetN(mx_h1) != 1))
	{
		fprintf(stderr, "J=%d L=%d tablelength=%d\n",
			Jd[0], Ld[0], (int) mxGetM(mx_h1));
		Fail("h1 size problem")
	}

	if ((int) mxGetM(mx_h2) != Jd[1]*Ld[1]+1 || mxGetN(mx_h2) != 1)
	{
		fprintf(stderr, "J=%d L=%d tablelength=%d\n",
			Jd[1], Ld[1], (int) mxGetM(mx_h2));
		Fail("h2 size problem")
	}

	if (mxGetN(mx_tm) != 2)
		Fail("tm must have 2 columns")

	// create a new array and set the output pointer to it
	Call(plhs[0] = mxCreateDoubleMatrix, (M, N, mxCOMPLEX))
	double *r_fm = mxGetPr(plhs[0]);
	double *i_fm = mxGetPi(plhs[0]);

	// call the C subroutine N times; once for each realization
	if (mxIsComplexDouble(mx_h1) && mxIsComplexDouble(mx_h2))
	{
		const double *i_h1 = mxGetPi(mx_h1);
		const double *i_h2 = mxGetPi(mx_h2);
		if (order)
			Fail("only 0th order implemented for complex")

		for (int nn=0; nn < N; ++nn)
		{
			interp2_table0_complex_per(r_ck, i_ck, K1, K2,
				r_h1, i_h1, r_h2, i_h2,
				Jd[0], Jd[1], Ld[0], Ld[1],
				p_tm, M, r_fm, i_fm);
			r_ck += K1*K2; i_ck += K1*K2;
			r_fm += M; i_fm += M;
		}
	}

	else if (mxIsRealDouble(mx_h1) && mxIsRealDouble(mx_h2))
	{
		interp2_table_real_per_t *fun;
		if (order == 0)		fun = interp2_table0_real_per;
		else if (order == 1)	fun = interp2_table1_real_per;
		else Fail("bad order")

#ifndef Provide_flip
                if (flips && (flips[0] || flips[1]))
			Fail("flip not compiled")
#endif

		for (int nn=0; nn < N; ++nn)
		{
			fun(r_ck, i_ck, K1, K2, r_h1, r_h2,
#ifdef Provide_flip
				flips ? flips[0] : 0,
				flips ? flips[1] : 0,
#endif
				Jd[0], Jd[1], Ld[0], Ld[1],
				p_tm, M, r_fm, i_fm);
			r_ck += K1*K2; i_ck += K1*K2;
			r_fm += M; i_fm += M;
		}
	}

	else
		Fail("h must be real or complex double (preferably real)")

	return 1;
}


// gateway
void mexFunction(
int nlhs, mxArray *plhs[],
int nrhs, const mxArray *prhs[])
{
	if (nlhs == 0 && nrhs == 1 && mxIsChar(prhs[0])) // check
		return;

	if (nrhs < 6 || nrhs > 8)
		mexFail("8 inputs required: (ck, h1, h2, J, L, tm, [order, flips])")
	if (nlhs > 1)
		mexFail("Less than one output arguments.")

	if (!interp2_table_mex(plhs,
		prhs[0], prhs[1], prhs[2], prhs[3], prhs[4], prhs[5],
		nrhs > 6 ? prhs[6] : NULL,
		nrhs > 7 ? prhs[7] : NULL))
		mexFail("interp2_table_mex failed")
}
