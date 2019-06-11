// interp3_table1_adj.c
// 3D periodic interpolation using table lookup: adjoint.
//
// forward direction: (for m = 1,...,M)
// f(t_m) = \double_sum_{k1,2,3=0}^{K1,2,3-1} c[k1,k2,k3]
//	h_1( (t_m1 - k1) mod K1 ) h_2( (t_m2 - k2) mod K2 ) h_3( (t_m3 - k3) mod K3 )
//
// adjoint direction: (for k1,2,3=0,...,K1,2,3-1) (note complex conjugate!)
// c[k1,k2,k3] = \sum_{m=1}^M f(t_m) h_1^*( (t_m1 - k1) mod K1 )
//		h_2^*( (t_m2 - k2) mod K2 ) h_3^*( (t_m3 - k3) mod K3 )
//
// Interpolators h1,2,3 are nonzero (and tabulated) for -J/2 <= t <= J/2.
//
// Copyright 2004-4-2 Jeff Fessler and Yingying Zhang, University of Michigan

#include <math.h>
#include <string.h>
#include "def,table3.h"


/*
* interp3_table0_complex_per_adj()
*/
void interp3_table0_complex_per_adj(
double *r_ck,		/* [K1,K2,K3] out */
double *i_ck,
const int K1,
const int K2,
const int K3,
const double *r_h1,	/* [J1*L1+1,1] in */
const double *i_h1,
const double *r_h2,	/* [J2*L2+1,1] in */
const double *i_h2,
const double *r_h3,	/* [J3*L3+1,1] in */
const double *i_h3,
const int J1,
const int J2,
const int J3,
const int L1,
const int L2,
const int L3,
const double *p_tm,	/* [M,3] in */
const int M,
const double *r_fm,	/* [M,1] in */
const double *i_fm)
{
	int mm;

	/* trick: shift table pointer to center */
	{
	const int ncenter1 = floor(J1 * L1/2.);
	r_h1 += ncenter1;
	i_h1 += ncenter1;
	}
	{
	const int ncenter2 = floor(J2 * L2/2.);
	r_h2 += ncenter2;
	i_h2 += ncenter2;
	}
	{
	const int ncenter3 = floor(J3 * L3/2.);
	r_h3 += ncenter3;
	i_h3 += ncenter3;
	}

	/* initialize output to zero */
	(void) memset((void *) r_ck, 0, K1*K2*K3*sizeof(*r_ck));
	(void) memset((void *) i_ck, 0, K1*K2*K3*sizeof(*i_ck));

	/* interp */
    for (mm=0; mm < M; mm++) {
	int jj1, jj2, jj3;
	const double t3 = p_tm[2*M];
	const double t2 = p_tm[M];
	const double t1 = *p_tm++;
	const double fmr = *r_fm++;
	const double fmi = *i_fm++;
	const int koff1 = 1 + floor(t1 - J1 / 2.);
	const int koff2 = 1 + floor(t2 - J2 / 2.);
	int k3 = 1 + floor(t3 - J3 / 2.);

	for (jj3=0; jj3 < J3; jj3++, k3++) {
		const double p3 = (t3 - k3) * L3;
		const int n3 = /* ncenter3 + */ iround(p3);
		const double coef3r = r_h3[n3];
		const double coef3i = i_h3[n3];
		const int k3mod = mymod(k3, K3);

		const double v3r = coef3r * fmr + coef3i * fmi;
		const double v3i = coef3r * fmi - coef3i * fmr;
		int k2 = koff2;

	for (jj2=0; jj2 < J2; jj2++, k2++) {
		const double p2 = (t2 - k2) * L2;
		const int n2 = /* ncenter2 + */ iround(p2);
		const double coef2r = r_h2[n2];
		const double coef2i = i_h2[n2];
		const int k2mod = mymod(k2, K2);
		const int k23mod = (k3mod * K2 + k2mod) * K1;

		const double v2r = coef2r * v3r + coef2i * v3i;
		const double v2i = coef2r * v3i - coef2i * v3r;
		int k1 = koff1;

	for (jj1=0; jj1 < J1; jj1++, k1++) {
		const double p1 = (t1 - k1) * L1;
		const int n1 = /* ncenter1 + */ iround(p1);
		register const double coef1r = r_h1[n1];
		register const double coef1i = i_h1[n1];
		const int k1mod = mymod(k1, K1);
		const int kk = k23mod + k1mod; /* 3D array index */

		r_ck[kk] += coef1r * v2r + coef1i * v2i;
		i_ck[kk] += coef1r * v2i - coef1i * v2r;

	} /* j1 */
	} /* j2 */
	} /* j3 */
    }
}


/*
* interp3_table0_real_per_adj()
* 3D, 0th-order, real, periodic
*/
void interp3_table0_real_per_adj(
double *r_ck,		/* [K1,K2,K3] out */
double *i_ck,
const int K1,
const int K2,
const int K3,
const double *r_h1,	/* [J1*L1+1,1] in */
const double *r_h2,	/* [J2*L2+1,1] in */
const double *r_h3,	/* [J3*L3+1,1] in */
#ifdef Provide_flip
const int flip1,	/* sign flips every K? */
const int flip2,
const int flip3,
#endif
const int J1,
const int J2,
const int J3,
const int L1,
const int L2,
const int L3,
const double *p_tm,	/* [M,3] in */
const int M,
const double *r_fm,	/* [M,1] in */
const double *i_fm)
{
	int mm;

	/* trick: shift table pointer to center */
	{
	const int ncenter1 = floor(J1 * L1/2.);
	r_h1 += ncenter1;
	}
	{
	const int ncenter2 = floor(J2 * L2/2.);
	r_h2 += ncenter2;
	}
	{
	const int ncenter3 = floor(J3 * L3/2.);
	r_h3 += ncenter3;
	}

	/* initialize output to zero */
	(void) memset((void *) r_ck, 0, K1*K2*K3*sizeof(*r_ck));
	(void) memset((void *) i_ck, 0, K1*K2*K3*sizeof(*i_ck));

	/* interp */
    for (mm=0; mm < M; mm++) {
	int jj1, jj2, jj3;
	const double t3 = p_tm[2*M];
	const double t2 = p_tm[M];
	const double t1 = *p_tm++;
	const double fmr = *r_fm++;
	const double fmi = *i_fm++;
	const int koff1 = 1 + floor(t1 - J1 / 2.);
	const int koff2 = 1 + floor(t2 - J2 / 2.);
	int k3 = 1 + floor(t3 - J3 / 2.);

	for (jj3=0; jj3 < J3; jj3++, k3++) {
		const double p3 = (t3 - k3) * L3;
		const int n3 = /* ncenter3 + */ iround(p3);
		double coef3r = r_h3[n3];
		const int wrap3 = floor(k3 / (double) K3);
		const int k3mod = k3 - K3 * wrap3;

#ifdef Provide_flip
		if (flip3 && (wrap3 % 2))
			coef3r = -coef3r; /* trick: sign flip */
#endif

		const double v3r = coef3r * fmr;
		const double v3i = coef3r * fmi;
		int k2 = koff2;

	for (jj2=0; jj2 < J2; jj2++, k2++) {
		const double p2 = (t2 - k2) * L2;
		const int n2 = /* ncenter2 + */ iround(p2);
		double coef2r = r_h2[n2];
		const int wrap2 = floor(k2 / (double) K2);
		const int k2mod = k2 - K2 * wrap2;
		const int k23mod = (k3mod * K2 + k2mod) * K1;

#ifdef Provide_flip
		if (flip2 && (wrap2 % 2))
			coef2r = -coef2r; /* trick: sign flip */
#endif

		const double v2r = coef2r * v3r;
		const double v2i = coef2r * v3i;
		int k1 = koff1;

	for (jj1=0; jj1 < J1; jj1++, k1++) {
		const double p1 = (t1 - k1) * L1;
		const int n1 = /* ncenter1 + */ iround(p1);
		register double coef1r = r_h1[n1];
		const int wrap1 = floor(k1 / (double) K1);
		const int k1mod = k1 - K1 * wrap1;
		const int kk = k23mod + k1mod;

#ifdef Provide_flip
		if (flip1 && (wrap1 % 2))
			coef1r = -coef1r; /* trick: sign flip */
#endif

		r_ck[kk] += coef1r * v2r;
		i_ck[kk] += coef1r * v2i;
	} /* j1 */
	} /* j2 */
	} /* j3 */
    }
}


/*
* interp3_table1_real_per_adj()
* 3D, 1st-order, real, periodic
*/
void interp3_table1_real_per_adj(
double *r_ck,		/* [K1,K2,K3] out */
double *i_ck,
const int K1,
const int K2,
const int K3,
const double *r_h1,	/* [J1*L1+1,1] in */
const double *r_h2,	/* [J2*L2+1,1] in */
const double *r_h3,	/* [J3*L3+1,1] in */
#ifdef Provide_flip
const int flip1,	/* sign flips every K? */
const int flip2,
const int flip3,
#endif
const int J1,
const int J2,
const int J3,
const int L1,
const int L2,
const int L3,
const double *p_tm,	/* [M,3] in */
const int M,
const double *r_fm,	/* [M,1] in */
const double *i_fm)
{
	int mm;

	/* trick: shift table pointer to center */
	{
	const int ncenter1 = floor(J1 * L1/2.);
	r_h1 += ncenter1;
	}
	{
	const int ncenter2 = floor(J2 * L2/2.);
	r_h2 += ncenter2;
	}
	{
	const int ncenter3 = floor(J3 * L3/2.);
	r_h3 += ncenter3;
	}

	/* initialize output to zero */
	(void) memset((void *) r_ck, 0, K1*K2*K3*sizeof(*r_ck));
	(void) memset((void *) i_ck, 0, K1*K2*K3*sizeof(*i_ck));

	/* interp */
    for (mm=0; mm < M; mm++) {
	int jj1, jj2, jj3;
	const double t3 = p_tm[2*M];
	const double t2 = p_tm[M];
	const double t1 = *p_tm++;
	const double fmr = *r_fm++;
	const double fmi = *i_fm++;
	const int koff1 = 1 + floor(t1 - J1 / 2.);
	const int koff2 = 1 + floor(t2 - J2 / 2.);
	int k3 = 1 + floor(t3 - J3 / 2.);

	for (jj3=0; jj3 < J3; jj3++, k3++) {
		const double p3 = (t3 - k3) * L3;
		const int n3 = floor(p3);
		const double alf3 = p3 - n3;
		register const double *ph3 = r_h3 + n3;
		double coef3r = (1 - alf3) * *ph3 + alf3 * *(ph3+1);
		const int wrap3 = floor(k3 / (double) K3);
		const int k3mod = k3 - K3 * wrap3;

#ifdef Provide_flip
		if (flip3 && (wrap3 % 2))
			coef3r = -coef3r; /* trick: sign flip */
#endif

		const double v3r = coef3r * fmr;
		const double v3i = coef3r * fmi;
		int k2 = koff2;

	for (jj2=0; jj2 < J2; jj2++, k2++) {
		const double p2 = (t2 - k2) * L2;
		const int n2 = floor(p2);
		const double alf2 = p2 - n2;
		register const double *ph2 = r_h2 + n2;
		double coef2r = (1 - alf2) * *ph2 + alf2 * *(ph2+1);
		const int wrap2 = floor(k2 / (double) K2);
		const int k2mod = k2 - K2 * wrap2;
		const int k23mod = (k3mod * K2 + k2mod) * K1;

#ifdef Provide_flip
		if (flip2 && (wrap2 % 2))
			coef2r = -coef2r; /* trick: sign flip */
#endif

		const double v2r = coef2r * v3r;
		const double v2i = coef2r * v3i;
		int k1 = koff1;

	for (jj1=0; jj1 < J1; jj1++, k1++) {
		const double p1 = (t1 - k1) * L1;
		const int n1 = floor(p1);
		const double alf1 = p1 - n1;
		register const double *ph1 = r_h1 + n1;
		register double coef1r = (1 - alf1) * *ph1 + alf1 * *(ph1+1);
		const int wrap1 = floor(k1 / (double) K1);
		const int k1mod = k1 - K1 * wrap1;
		const int kk = k23mod + k1mod;

#ifdef Provide_flip
		if (flip1 && (wrap1 % 2))
			coef1r = -coef1r; /* trick: sign flip */
#endif

		r_ck[kk] += coef1r * v2r;
		i_ck[kk] += coef1r * v2i;
	} /* j1 */
	} /* j2 */
	} /* j3 */
    }
}
