// interp1_table1_for.c
// 1D periodic interpolation using table lookup.
//
// forward direction: (for m = 1,...,M)
// f(t_m) = \sum_{k=0}^{K-1} c_k h( (t_m - k) mod K )
//
// The interpolator h is nonzero (and tabulated) for -J/2 <= t <= J/2.
//
// Copyright 2004-3 Yingying Zhang and Jeff Fessler, University of Michigan

#include <math.h>
#include <stdio.h>
#include "def,table1.h"

/*
* interp1_table0_complex_per()
* 1D, 0th order, complex, periodic
*/
void interp1_table0_complex_per(
const double *r_ck,	/* [K1,1] in */
const double *i_ck,
const int K1,
const double *r_h1,	/* [J1*L1+1,1] in */
const double *i_h1,	/* imaginary part of complex interpolator */
const int J1,
const int L1,
const double *p_tm,	/* [M,1] in */
const int M,
double *r_fm,		/* [M,1] out */
double *i_fm)
{
	int mm;

	/* trick: shift table pointer to center */
	{
	const int ncenter1 = floor(J1 * L1/2.);
	r_h1 += ncenter1;
	i_h1 += ncenter1;
	}

	/* interp */
    for (mm=0; mm < M; mm++) {
	int jj1;
	const double t1 = *p_tm++;
	register double sum1r = 0;
	register double sum1i = 0;
	int k1 = 1 + floor(t1 - J1 / 2.);

	for (jj1=0; jj1 < J1; jj1++, k1++) {
		const double p1 = (t1 - k1) * L1;
		const int n1 = /* ncenter1 + */ iround(p1);
		register const double coef1r = r_h1[n1];
		register const double coef1i = i_h1[n1];
		const int k1mod = mymod(k1, K1);

		/* sum1 += coef1 * ck */
		sum1r += coef1r * r_ck[k1mod] - coef1i * i_ck[k1mod];
		sum1i += coef1r * i_ck[k1mod] + coef1i * r_ck[k1mod];
	}
	*r_fm++ = sum1r;
	*i_fm++ = sum1i;
    }
}


/*
* interp1_table0_real_per()
* 1D, 0th-order, real, periodic
*/
void interp1_table0_real_per(
const double *r_ck,	/* [K,1] in */
const double *i_ck,
const int K1,
const double *r_h1,	/* [J1*L1+1,1] in (real) */
#ifdef Provide_flip
const int flip1,	/* sign flips every K? */
#endif
const int J1,
const int L1,
const double *p_tm,	/* [M,1] in */
const int M,
double *r_fm,		/* [M,1] out */
double *i_fm)
{
	int mm;

	/* trick: shift table pointer to center */
	{
	const int ncenter1 = floor(J1 * L1/2.);
	r_h1 += ncenter1;
	}

	/* interp */
    for (mm=0; mm < M; mm++) {
	int jj1;
	const double t1 = *p_tm++;
	register double sum1r = 0;
	register double sum1i = 0;
	int k1 = 1 + floor(t1 - J1 / 2.);

	for (jj1=0; jj1 < J1; jj1++, k1++) {
		const double p1 = (t1 - k1) * L1;
		const int n1 = /* ncenter1 + */ iround(p1);
		register double coef1r = r_h1[n1];
		const int wrap1 = floor(k1 / (double) K1);
		const int k1mod = k1 - K1 * wrap1;
#ifdef Provide_flip
		if (flip1 && (wrap1 % 2))
			coef1r = -coef1r; /* trick: sign flip */
#endif

		/* sum1 += coef1 * ck */
		sum1r += coef1r * r_ck[k1mod];
		sum1i += coef1r * i_ck[k1mod];
	}
	*r_fm++ = sum1r;
	*i_fm++ = sum1i;
    }
}


/*
* interp1_table1_real_per()
* 1D, 1st-order, real, periodic
*/
void interp1_table1_real_per(
const double *r_ck,	/* [K,1] in */
const double *i_ck,
const int K1,
const double *r_h1,	/* [J1*L1+1,1] in (real) */
#ifdef Provide_flip
const int flip1,	/* sign flips every K? */
#endif
const int J1,
const int L1,
const double *p_tm,	/* [M,1] in */
const int M,
double *r_fm,		/* [M,1] out */
double *i_fm)
{
	int mm;

	/* trick: shift table pointer to center */
	{
	const int ncenter1 = floor(J1 * L1/2.);
	r_h1 += ncenter1;
	}

	/* interp */
    for (mm=0; mm < M; mm++) {
	int jj1;
	const double t1 = *p_tm++;
	register double sum1r = 0;
	register double sum1i = 0;
	int k1 = 1 + floor(t1 - J1 / 2.);

	for (jj1=0; jj1 < J1; jj1++, k1++) {
		const double p1 = (t1 - k1) * L1;
		const int n1 = floor(p1);
		const double alf1 = p1 - n1;
		register const double *ph1 = r_h1 + n1;
		register double coef1r = (1 - alf1) * *ph1 + alf1 * *(ph1+1);
		const int wrap1 = floor(k1 / (double) K1);
		const int k1mod = k1 - K1 * wrap1;
#ifdef Provide_flip
		if (flip1 && (wrap1 % 2))
			coef1r = -coef1r; /* trick: sign flip */
#endif

		/* sum1 += coef1 * ck */
		sum1r += coef1r * r_ck[k1mod];
		sum1i += coef1r * i_ck[k1mod];
	}
	*r_fm++ = sum1r;
	*i_fm++ = sum1i;
    }
}
