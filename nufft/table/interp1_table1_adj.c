// interp1_table1_adj.c
// adjoint of 1D periodic interpolation using table lookup.
//
// forward direction: (for m = 1,...,M)
// f(t_m) = \sum_{k=0}^{K-1} c_k h( (t_m - k) mod K )
//
// adjoint direction: (for k=0,...,K-1) (note complex conjugate!)
// c_k = \sum_{m=1}^M f(t_m) h^*( (t_m - k) mod K )
//
// The interpolator h is nonzero (and tabulated) for -J/2 <= t <= J/2.
//
// Copyright 2004-3 Yingying Zhang and Jeff Fessler, University of Michigan

#include <math.h>
#include "def,table1.h"


// interp1_table0_complex_per_adj()
// 1D, 0th order, complex, periodic
void interp1_table0_complex_per_adj(
double *r_ck, // [K1 1] out
double *i_ck,
const int K1,
const double *r_h1,	// [J1*L1+1,1] in
const double *i_h1,	// imaginary part of complex interpolator
const int J1,
const int L1,
const double *p_tm,	// [M 1] in
const int M,
const double *r_fm,	// [M 1] in
const double *i_fm)
{
	// initialize output to zero
	(void) memset((void *) r_ck, 0, K1*sizeof(*r_ck));
	(void) memset((void *) i_ck, 0, K1*sizeof(*i_ck));

	// trick: shift table pointer to center
	{
	const int ncenter1 = floor(J1 * L1/2.);
	r_h1 += ncenter1;
	i_h1 += ncenter1;
	}

	// interp
    for (int mm=0; mm < M; mm++) {
	const double t1 = *p_tm++;
	const double fmr = *r_fm++;
	const double fmi = *i_fm++;
	int k1 = 1 + floor(t1 - J1 / 2.);

	for (int jj1=0; jj1 < J1; jj1++, k1++)
	{
		const double p1 = (t1 - k1) * L1;
		const int n1 = /* ncenter1 + */ iround(p1);
		register const double coef1r = r_h1[n1];
		register const double coef1i = i_h1[n1];
		const int k1mod = mymod(k1, K1);

		// instead of f = h c, we have c += h^* f
		r_ck[k1mod] += coef1r * fmr + coef1i * fmi;
		i_ck[k1mod] += coef1r * fmi - coef1i * fmr;
	}
    }
}


// interp1_table0_real_per()
// 1D, 0th-order, real, periodic
void interp1_table0_real_per_adj(
double *r_ck, // [K1,1] out
double *i_ck,
const int K1,
const double *r_h1,	// [J1*L1+1,1] in (real)
#ifdef Provide_flip
const int flip1,	// sign flips every K?
#endif
const int J1,
const int L1,
const double *p_tm,	// [M,1] in
const int M,
const double *r_fm,	// [M,1] in
const double *i_fm)
{
	// initialize output to zero
	(void) memset((void *) r_ck, 0, K1*sizeof(*r_ck));
	(void) memset((void *) i_ck, 0, K1*sizeof(*i_ck));

	// trick: shift table pointer to center
	{
	const int ncenter1 = floor(J1 * L1/2.);
	r_h1 += ncenter1;
	}

	// interp
    for (int mm=0; mm < M; mm++) {
	const double t1 = *p_tm++;
	const double fmr = *r_fm++;
	const double fmi = *i_fm++;
	int k1 = 1 + floor(t1 - J1 / 2.);

	for (int jj1=0; jj1 < J1; jj1++, k1++)
	{
		const double p1 = (t1 - k1) * L1;
		const int n1 = /* ncenter1 + */ iround(p1);
		register double coef1r = r_h1[n1];
		const int wrap1 = floor(k1 / (double) K1);
		const int k1mod = k1 - K1 * wrap1;
#ifdef Provide_flip
		if (flip1 && (wrap1 % 2))
			coef1r = -coef1r; // trick: sign flip
#endif

		// instead of f = h c, we have c += h^* f
		r_ck[k1mod] += coef1r * fmr;
		i_ck[k1mod] += coef1r * fmi;
	}
    }
}


// interp1_table1_real_per()
// 1D, 1st-order, real, periodic
void interp1_table1_real_per_adj(
double *r_ck,		// [K1 1] out
double *i_ck,
const int K1,
const double *r_h1,	// [J1*L1+1,1] in (real)
#ifdef Provide_flip
const int flip1,	// sign flips every K?
#endif
const int J1,
const int L1,
const double *p_tm,	// [M 1] in
const int M,
const double *r_fm,	// [M 1] in
const double *i_fm)
{
	// initialize output to zero
	(void) memset((void *) r_ck, 0, K1*sizeof(*r_ck));
	(void) memset((void *) i_ck, 0, K1*sizeof(*i_ck));

	// trick: shift table pointer to center
	{
	const int ncenter1 = floor(J1 * L1/2.);
	r_h1 += ncenter1;
	}

	// interp
    for (int mm=0; mm < M; mm++) {
	const double t1 = *p_tm++;
	const double fmr = *r_fm++;
	const double fmi = *i_fm++;
	int k1 = 1 + floor(t1 - J1 / 2.);

	for (int jj1=0; jj1 < J1; jj1++, k1++)
	{
		const double p1 = (t1 - k1) * L1;
		const int n1 = floor(p1);
		const double alf1 = p1 - n1;
		register const double *ph1 = r_h1 + n1;
		register double coef1r = (1 - alf1) * *ph1 + alf1 * *(ph1+1);
		const int wrap1 = floor(k1 / (double) K1);
		const int k1mod = k1 - K1 * wrap1;
#ifdef Provide_flip
		if (flip1 && (wrap1 % 2))
			coef1r = -coef1r; // trick: sign flip
#endif

		// instead of f = h c, we have c += h^* f
		r_ck[k1mod] += coef1r * fmr;
		i_ck[k1mod] += coef1r * fmi;
	}
    }
}
