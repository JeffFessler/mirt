// interp4_table1_adj.c
// 4D periodic interpolation using table lookup: adjoint.
//
// forward direction: (for m = 1,...,M)
// f(t_m) = \double_sum_{k1,2,3,4=0}^{K1,2,3,4-1} c[k1,k2,k3,k4]
//	h_1( (t_m1 - k1) mod K1 ) h_2( (t_m2 - k2) mod K2 ) h_3( (t_m3 - k3) mod K3 ) h_4( (t_m4 - k4) mod K4 )
//
// adjoint direction: (for k1,2,3,4=0,...,K1,2,3,4-1) (note complex conjugate!)
// c[k1,k2,k3,k4] = \sum_{m=1}^M f(t_m) h_1^*( (t_m1 - k1) mod K1 )
//		h_2^*( (t_m2 - k2) mod K2 ) h_3^*( (t_m3 - k3) mod K3 ) h_4^*( (t_m4 - k4) mod K4 )
//
// Interpolators h1,2,3,4 are nonzero (and tabulated) for -J/2 <= t <= J/2.
//
// Copyright 2004-4-2 Jeff Fessler and Yingying Zhang, University of Michigan

// Extended 2013-08-02 to 4D by David Johnson, The Ohio State University Wexner Medical Center

#include <math.h>
#include <string.h>
#include "def,table4.h"


#ifdef Is_pc
# include "def,mwSizePatch.h" // for mwIndex on PC
#else
# define mwIndex long
#endif


// interp4_table0_complex_per_adj()
void interp4_table0_complex_per_adj(
double *r_ck,		/* [K1,K2,K3,K4] out */
double *i_ck,
const int K1,
const int K2,
const int K3,
const int K4,
const double *r_h1,	/* [J1*L1+1,1] in */
const double *i_h1,
const double *r_h2,	/* [J2*L2+1,1] in */
const double *i_h2,
const double *r_h3,	/* [J3*L3+1,1] in */
const double *i_h3,
const double *r_h4,	/* [J4*L4+1,1] in */
const double *i_h4,
const int J1,
const int J2,
const int J3,
const int J4,
const int L1,
const int L2,
const int L3,
const int L4,
const double *p_tm,	/* [M,4] in */
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
	{
	const int ncenter4 = floor(J4 * L4/2.);
	r_h4 += ncenter4;
	i_h4 += ncenter4;
	}

	/* initialize output to zero */
	(void) memset((void *) r_ck, 0, K1*K2*K3*K4*sizeof(*r_ck));
	(void) memset((void *) i_ck, 0, K1*K2*K3*K4*sizeof(*i_ck));

	/* interp */
    for (mm=0; mm < M; mm++) {
	int jj1, jj2, jj3, jj4;
	const double t4 = p_tm[3*M];
	const double t3 = p_tm[2*M];
	const double t2 = p_tm[M];
	const double t1 = *p_tm++;
	const double fmr = *r_fm++;
	const double fmi = *i_fm++;
	const int koff1 = 1 + floor(t1 - J1 / 2.);
	const int koff2 = 1 + floor(t2 - J2 / 2.);
	const int koff3 = 1 + floor(t3 - J3 / 2.);
	int k4 = 1 + floor(t4 - J4 / 2.);

	for (jj4=0; jj4 < J4; jj4++, k4++) {
		const double p4 = (t4 - k4) * L4;
		const int n4 = /* ncenter4 + */ iround(p4);
		const double coef4r = r_h4[n4];
		const double coef4i = i_h4[n4];
		const int k4mod = mymod(k4, K4);
		const mwIndex k4Idx = k4mod*K3*K2*K1;

		const double v4r = coef4r * fmr + coef4i * fmi;
		const double v4i = coef4r * fmi - coef4i * fmr;
		int k3 = koff3;

	for (jj3=0; jj3 < J3; jj3++, k3++) {
		const double p3 = (t3 - k3) * L3;
		const int n3 = /* ncenter3 + */ iround(p3);
		const double coef3r = r_h3[n3];
		const double coef3i = i_h3[n3];
		const int k3mod = mymod(k3, K3);
		const mwIndex k3Idx = k3mod*K2*K1;

		const double v3r = coef3r * v4r + coef3i * v4i;
		const double v3i = coef3r * v4i - coef3i * v4r;
		int k2 = koff2;

	for (jj2=0; jj2 < J2; jj2++, k2++) {
		const double p2 = (t2 - k2) * L2;
		const int n2 = /* ncenter2 + */ iround(p2);
		const double coef2r = r_h2[n2];
		const double coef2i = i_h2[n2];
		const int k2mod = mymod(k2, K2);
		const mwIndex k2Idx = k2mod*K1;

		const double v2r = coef2r * v3r + coef2i * v3i;
		const double v2i = coef2r * v3i - coef2i * v3r;
		int k1 = koff1;

	for (jj1=0; jj1 < J1; jj1++, k1++) {
		const double p1 = (t1 - k1) * L1;
		const int n1 = /* ncenter1 + */ iround(p1);
		register const double coef1r = r_h1[n1];
		register const double coef1i = i_h1[n1];
		const int k1mod = mymod(k1, K1);
		const mwIndex kk = k4Idx + k3Idx + k2Idx + k1mod; /* 4D array index */

		r_ck[kk] += coef1r * v2r + coef1i * v2i;
		i_ck[kk] += coef1r * v2i - coef1i * v2r;

	} /* j1 */
	} /* j2 */
	} /* j3 */
	} /* j4 */
    }
}


/*
* interp4_table0_real_per_adj()
* 4D, 0th-order, real, periodic
*/
void interp4_table0_real_per_adj(
double *r_ck,		/* [K1,K2,K3,K4] out */
double *i_ck,
const int K1,
const int K2,
const int K3,
const int K4,
const double *r_h1,	/* [J1*L1+1,1] in */
const double *r_h2,	/* [J2*L2+1,1] in */
const double *r_h3,	/* [J3*L3+1,1] in */
const double *r_h4,	/* [J4*L4+1,1] in */
#ifdef Provide_flip
const int flip1,	/* sign flips every K? */
const int flip2,
const int flip3,
const int flip4,
#endif
const int J1,
const int J2,
const int J3,
const int J4,
const int L1,
const int L2,
const int L3,
const int L4,
const double *p_tm,	/* [M,4] in */
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
	{
	const int ncenter4 = floor(J4 * L4/2.);
	r_h4 += ncenter4;
	}

	/* initialize output to zero */
	(void) memset((void *) r_ck, 0, K1*K2*K3*K4*sizeof(*r_ck));
	(void) memset((void *) i_ck, 0, K1*K2*K3*K4*sizeof(*i_ck));

	/* interp */
    for (mm=0; mm < M; mm++) {
	int jj1, jj2, jj3, jj4;
	const double t4 = p_tm[3*M];
	const double t3 = p_tm[2*M];
	const double t2 = p_tm[M];
	const double t1 = *p_tm++;
	const double fmr = *r_fm++;
	const double fmi = *i_fm++;
	const int koff1 = 1 + floor(t1 - J1 / 2.);
	const int koff2 = 1 + floor(t2 - J2 / 2.);
	const int koff3 = 1 + floor(t3 - J3 / 2.);
	int k4 = 1 + floor(t4 - J4 / 2.);

	for (jj4=0; jj4 < J4; jj4++, k4++) {
		const double p4 = (t4 - k4) * L4;
		const int n4 = /* ncenter4 + */ iround(p4);
		double coef4r = r_h4[n4];
		const int wrap4 = floor(k4 / (double) K4);
		const int k4mod = k4 - K4 * wrap4;
		const mwIndex k4Idx = k4mod*K3*K2*K1;

#ifdef Provide_flip
		if (flip4 && (wrap4 % 2))
			coef4r = -coef4r; /* trick: sign flip */
#endif

		const double v4r = coef4r * fmr;
		const double v4i = coef4r * fmi;
		int k3 = koff3;

	for (jj3=0; jj3 < J3; jj3++, k3++) {
		const double p3 = (t3 - k3) * L3;
		const int n3 = /* ncenter3 + */ iround(p3);
		double coef3r = r_h3[n3];
		const int wrap3 = floor(k3 / (double) K3);
		const int k3mod = k3 - K3 * wrap3;
		const mwIndex k3Idx = k3mod*K2*K1;

#ifdef Provide_flip
		if (flip3 && (wrap3 % 2))
			coef3r = -coef3r; /* trick: sign flip */
#endif

		const double v3r = coef3r * v4r;
		const double v3i = coef3r * v4i;
		int k2 = koff2;

	for (jj2=0; jj2 < J2; jj2++, k2++) {
		const double p2 = (t2 - k2) * L2;
		const int n2 = /* ncenter2 + */ iround(p2);
		double coef2r = r_h2[n2];
		const int wrap2 = floor(k2 / (double) K2);
		const int k2mod = k2 - K2 * wrap2;
		const mwIndex k2Idx = k2mod*K1;

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
		const mwIndex kk = k4Idx + k3Idx + k2Idx + k1mod; /* 4D array index */

#ifdef Provide_flip
		if (flip1 && (wrap1 % 2))
			coef1r = -coef1r; /* trick: sign flip */
#endif

		r_ck[kk] += coef1r * v2r;
		i_ck[kk] += coef1r * v2i;
	} /* j1 */
	} /* j2 */
	} /* j3 */
	} /* j4 */
    }
}


/*
* interp4_table1_real_per_adj()
* 4D, 1st-order, real, periodic
*/
void interp4_table1_real_per_adj(
double *r_ck,		/* [K1,K2,K3,K4] out */
double *i_ck,
const int K1,
const int K2,
const int K3,
const int K4,
const double *r_h1,	/* [J1*L1+1,1] in */
const double *r_h2,	/* [J2*L2+1,1] in */
const double *r_h3,	/* [J3*L3+1,1] in */
const double *r_h4,	/* [J4*L4+1,1] in */
#ifdef Provide_flip
const int flip1,	/* sign flips every K? */
const int flip2,
const int flip3,
const int flip4,
#endif
const int J1,
const int J2,
const int J3,
const int J4,
const int L1,
const int L2,
const int L3,
const int L4,
const double *p_tm,	/* [M,4] in */
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
	{
	const int ncenter4 = floor(J4 * L4/2.);
	r_h4 += ncenter4;
	}

	/* initialize output to zero */
	(void) memset((void *) r_ck, 0, K1*K2*K3*K4*sizeof(*r_ck));
	(void) memset((void *) i_ck, 0, K1*K2*K3*K4*sizeof(*i_ck));

	/* interp */
    for (mm=0; mm < M; mm++) {
	int jj1, jj2, jj3, jj4;
	const double t4 = p_tm[3*M];
	const double t3 = p_tm[2*M];
	const double t2 = p_tm[M];
	const double t1 = *p_tm++;
	const double fmr = *r_fm++;
	const double fmi = *i_fm++;
	const int koff1 = 1 + floor(t1 - J1 / 2.);
	const int koff2 = 1 + floor(t2 - J2 / 2.);
	const int koff3 = 1 + floor(t3 - J3 / 2.);
	int k4 = 1 + floor(t4 - J4 / 2.);

	for (jj4=0; jj4 < J4; jj4++, k4++) {
		const double p4 = (t4 - k4) * L4;
		const int n4 = floor(p4);
		const double alf4 = p4 - n4;
		register const double *ph4 = r_h4 + n4;
		double coef4r = (1 - alf4) * *ph4 + alf4 * *(ph4+1);
		const int wrap4 = floor(k4 / (double) K4);
		const int k4mod = k4 - K4 * wrap4;
		const mwIndex k4Idx = k4mod*K3*K2*K1;

#ifdef Provide_flip
		if (flip4 && (wrap4 % 2))
			coef4r = -coef4r; /* trick: sign flip */
#endif

		const double v4r = coef4r * fmr;
		const double v4i = coef4r * fmi;
		int k3 = koff3;

	for (jj3=0; jj3 < J3; jj3++, k3++) {
		const double p3 = (t3 - k3) * L3;
		const int n3 = floor(p3);
		const double alf3 = p3 - n3;
		register const double *ph3 = r_h3 + n3;
		double coef3r = (1 - alf3) * *ph3 + alf3 * *(ph3+1);
		const int wrap3 = floor(k3 / (double) K3);
		const int k3mod = k3 - K3 * wrap3;
		const mwIndex k3Idx = k3mod*K2*K1;

#ifdef Provide_flip
		if (flip3 && (wrap3 % 2))
			coef3r = -coef3r; /* trick: sign flip */
#endif

		const double v3r = coef3r * v4r;
		const double v3i = coef3r * v4i;
		int k2 = koff2;

	for (jj2=0; jj2 < J2; jj2++, k2++) {
		const double p2 = (t2 - k2) * L2;
		const int n2 = floor(p2);
		const double alf2 = p2 - n2;
		register const double *ph2 = r_h2 + n2;
		double coef2r = (1 - alf2) * *ph2 + alf2 * *(ph2+1);
		const int wrap2 = floor(k2 / (double) K2);
		const int k2mod = k2 - K2 * wrap2;
		const mwIndex k2Idx = k2mod*K1;

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
		const mwIndex kk = k4Idx + k3Idx + k2Idx + k1mod; /* 4D array index */

#ifdef Provide_flip
		if (flip1 && (wrap1 % 2))
			coef1r = -coef1r; /* trick: sign flip */
#endif

		r_ck[kk] += coef1r * v2r;
		i_ck[kk] += coef1r * v2i;
	} /* j1 */
	} /* j2 */
	} /* j3 */
	} /* j4 */
    }
}
