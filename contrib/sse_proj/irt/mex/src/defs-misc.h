/*
* defs-misc.h
*
* Copyright 1993, Jeff Fessler, The University of Michigan
*/
#ifndef DefMisc
#define DefMisc

#include "defs-env.h"
#include <limits.h>

#define PoissonBreak   40.
#define Poisson(mean)   ((mean) ? \
	(((mean) < PoissonBreak) ? poisson1(mean) : poisson2(mean)) : 0)



/*
* FFT often uses smallest power of 2 that is >= n
*/
#define Power2(n)	((n) <= 1 ? 1 : Rint(pow(2., Ceil(log(n-.1) / log(2)))))


/*
* Complex numbers for FFT's
*/
#ifndef ComplexPrec
#define ComplexPrec	8
#endif

#if ComplexPrec == 4
#	define ComplexType	float
#elif ComplexPrec == 8
#	define ComplexType	double
#else
	This is a bug
#endif

typedef struct {
	ComplexType	r;	/* real	*/
	ComplexType	i;	/* imaginary */
} Complex;

#ifndef Fourn
#if ComplexPrec == 4
#define Fourn(cx, nn, nd, pm)   fourn_float(((float *) (cx))-1, nn-1, nd, pm)
#elif ComplexPrec == 8
#define Fourn(cx, nn, nd, pm)   fourn_double(((double *) (cx))-1, nn-1, nd, pm)
#else
	This is a bug
#endif
#endif

#define f2power(x,p)	\
	( (p == -1.) ? (1. / (x)) : ((p == -1./2.) ? sqrt(1./(x)) : pow(x,p)) )

/*
* Pair of 2D filters
*/
typedef	struct {
	int	nx;
	int	ny;
	Complex	*p;	/* [nx,ny] */
} Filters2;

/*
* Pair of 3D filters
*/
typedef	struct {
	int	nx;
	int	ny;
	int	nz;
	Complex	*p;	/* [nx,ny,nz] */
} Filters3;

/*
* How to handle negative (or zero) coefficients for A'A
*/
typedef enum {
	FFTnegZero,	/* zero resulting filter coefficient */
	FFTnegTruncate,	/* truncate A'A part to zero, use 1 / (0 + R) */
	FFTnegMinPos	/* use minimum positive value */
} FFTnegType;

/* complex.c */
extern void cx_put_real(Complex *pc, cfloat *pr, int nn);
extern void cx_put_imag(Complex *pc, cfloat *pi, int nn);
extern void cx_get_real(float *pr, Const Complex *pc, int nn);
extern void cx_get_imag(float *pi, Const Complex *pc, int nn);

extern void cx_put_real_shift2_pad(Complex *, float *, cint, cint, cfloat *, cint, cint);
extern void cx_put_imag_shift2_pad(Complex *, float *, cint, cint, cfloat *, cint, cint);

extern void cx_put3(Complex *, cint, cint, cfloat *, cfloat *, cint, cint, cint);
extern void cx_get3(float *, float *, cint, cint, cint, Const Complex *, cint, cint);
extern void cx_put_real_shift3(Complex *pc, cfloat *pr, cint nx, cint ny, cint);
extern void cx_put_imag_shift3(Complex *pc, cfloat *pr, cint nx, cint ny, cint);
extern void cx_put_real_shift3_pad(Complex *, float *, cint, cint, cint, cfloat *, cint, cint, cint);
extern void cx_put_imag_shift3_pad(Complex *, float *, cint, cint, cint, cfloat *, cint, cint, cint);
#if 0
extern void cx_get_real_shift3(float *pr, Complex *pc, cint nx, cint ny, cint);
#endif

/* f2,filt.c */
extern Filters2 *f2_alloc(cint nx, cint ny, cint chat);
extern bool f2_free(Filters2 *f2);
extern bool f2_set_filters2(Filters2 *, cfloat *, cfloat *, cint, cint, cbool);
extern bool f2_filter_image2(float *, float *, cint, cint, Filters2 *, cint);

/* f2,op.c */
extern bool f2_negtype(FFTnegType *, cchar *);
extern void f2_mul2(Complex *, Complex *, Complex *, cint, cint);
extern int  f2_inv2(Complex *, Complex *, cdouble, cint, cint);
extern bool	f2_inv2_betas(Complex *, Complex *, cdouble, cint, cint, cdouble, cdouble, FFTnegType, cint);
extern bool f2_div2(Complex *, Complex *, Complex *, cint, cint);

/* f2,qpuls.c */
extern bool f2_quad_penalty(float *reg, cint nx, cint ny, cint hood, cdouble beta);
extern bool f2_qpuls_filters2(Filters2 *, FFTnegType, cfloat *, cint, cint, cdouble, cdouble, cint, cdouble, cint);
extern bool f2_qpuls_decon(float *, cfloat *, cint, cint, cfloat *, cint, cint,
	cdouble, cint, FFTnegType, cint);

/* f3,filt.c */
extern void f3_psf_sym(float *, cint, cint, cint);
extern bool f3_psf_sym_check(cfloat *, cint, cint, cint);
extern Filters3 *f3_alloc(cint, cint, cint, cint chat);
extern bool f3_free(Filters3 *f3);
extern bool f3_set_filters2(Filters3 *, cfloat *, cfloat *, cint, cint, cint, cbool);
extern bool f3_filter_image2(float *, float *, cint, cint, cint, Filters3 *, cint);

/* f3,op.c */
extern void f3_mul2(Complex *, Complex *, Complex *, cint, cint, cint);
extern int  f3_inv2(Complex *, Complex *, cdouble, cint, cint, cint);
extern bool	f3_inv2_betas(Complex *, Complex *, cdouble, cint, cint, cint, cdouble, cdouble, FFTnegType, cint);

/* f3,qpuls.c */
extern bool f3_qpuls_filters2(Filters3 *, FFTnegType, cfloat *, cint, cint, cint,
	cint, cdouble, cdouble, cdouble, cdouble, cdouble, cint);
extern bool f3_qpuls_decon(float *, cfloat *, cint, cint, cint,
	cfloat *, cint, cint, cint, cint, cdouble, cdouble, FFTnegType, cint);


/* four1.c fourn.c */
extern void four1(float *data, int nn, int isign);
extern void fourn_float(float *, cint *, cint ndim, cint isign);
extern void fourn_double(double *, cint *, cint ndim, cint isign);

/* filter-sym.c */
extern void filter_sym_x(float *, cfloat *, cfloat *, int, int, int);
extern void filter_sym_x_per(float *, cfloat *, cfloat *, int, int, int);
extern void filter_sym_y(float *, cfloat *, cfloat *, int, int, int);
extern void filter_sym_y_per(float *, cfloat *, cfloat *, int, int, int);
extern void filter_sym_edge_z(float *, cfloat *, cfloat *, cint, cint, cint, cint);
extern void filter_sym_z_zero(float *, cfloat *, cfloat *, cint, cint, cint, cint);
extern void filter_sym_z_renorm(float *, cfloat *, cfloat *, cint, cint, cint, cint);

#if 0
extern void filter_pad_input(float *, float *, float *, int, int);
extern void filter_pad_output(float *, float *, float *, int, int);
extern void filter_pad(float *, float *, float *, int, int, bool);
#endif

/* fwhm.c */
extern bool vect_imax_float(int *p_imax, cfloat *data, cint nn);
extern bool fwhm_print_1d(FILE *, cfloat *, cfloat *, cdouble, cint, cchar *, cint);
extern bool fwhm_print_2d(FILE *, cfloat *, cfloat *, cdouble, cint, cint, cchar *, cint);
extern bool fwhm_print_3d(FILE *, cfloat *, cfloat *, cdouble,
	cint, cint, cint, cchar *, cint);


/* gauss.c */
extern double gauss(void);

/* gauss-psf.c */
extern bool gauss_half_psf(float *, cint, cdouble, cint);

/* nearest.c */
extern bool nearest_class(byte *, float *, cfloat *, cfloat *, int, cint);

/* penalty,diff.c */
#if 0
extern bool penalty_diff_wk_tight_offset(float *, cfloat *, cint, cint, cint, cint, cdouble);
extern bool penalty_diff_wk_leak_offset(float *, cfloat *, cint, cint, cint, cint, cdouble);
extern bool penalty_diff1_wk_tight(float *, cfloat *, cint *, cint, cint *, cint, cdouble);
extern bool penalty_diff1_wk_leak(float *, cfloat *, cint *, cint, cint *, cint, cdouble);
extern bool penalty_diff2_wk_tight(float *, cfloat *, cint *, cint, cint *, cint, cdouble);
extern bool penalty_diff2_wk_leak(float *, cfloat *, cint *, cint, cint *, cint, cdouble);
#endif
extern bool penalty_diff_wk_tight(float *, cfloat *, cint *, cint, cint *, cint, cdouble, cint);
extern bool penalty_diff_wk_leak(float *, cfloat *, cint *, cint, cint *, cint, cdouble, cint);
#if 0
extern void penalty_diff1_forw_offset(float *, cfloat *, cint, cint);
extern void penalty_diff1_forw2_offset(float *, cfloat *, cint, cint);
extern void penalty_diff1_back_offset(float *, cfloat *, cint, cint);
extern void penalty_diff1_back2_offset(float *, cfloat *, cint, cint);
extern void penalty_diff1_forw1(float *, cfloat *, cint, cint *, cint, cint);
extern void penalty_diff2_forw1(float *, cfloat *, cint, cint *, cint, cint);
extern void penalty_diff1_forw2(float *, cfloat *, cint, cint *, cint, cint);
extern void penalty_diff2_forw2(float *, cfloat *, cint, cint *, cint, cint);
extern void penalty_diff1_back1(float *, cfloat *, cint, cint *, cint, cint);
extern void penalty_diff2_back1(float *, cfloat *, cint, cint *, cint, cint);
extern void penalty_diff1_back2(float *, cfloat *, cint, cint *, cint, cint);
extern void penalty_diff2_back2(float *, cfloat *, cint, cint *, cint, cint);
#endif
extern bool penalty_diff_forw(float *, cfloat *, cint, cint *, cint, cint, cint, cint);
extern bool penalty_diff_back(float *, cfloat *, cint, cint *, cint, cint, cint, cint);

/* pick-rect.c */
extern void pick_rect_float(float *,	cfloat *, cint, cint, cint, cint, cint, cint, cint, cdouble);
extern void pick_rect_double(double *,	cfloat *, cint, cint, cint, cint, cint, cint, cint, cdouble);

/* poisson.c */
extern void poisson(float *, cfloat *, int, cdouble, cdouble);
extern unsigned long poisson1(cdouble);
extern unsigned long poisson2(cdouble);

/* rotate0.c */
extern void rotate0(float *, cfloat *, cint, cint, cdouble);

/* shuffle.c */
extern void transpose_0_1(float *po, cfloat *pi, cint nx, cint ny, cint nz);
extern void transpose_0_2(float *po, cfloat *pi, cint nx, cint ny, cint nz);
extern void transpose_1_2(float *po, cfloat *pi, cint nx, cint ny, cint nz);
extern void transpose_1_2_inplace(float *, float *, cint, cint, cint);
extern bool transpose_2_3(float *po, cfloat *pi, cint, cint, cint, cint);
extern bool transpose3(float *, cfloat *, cint, cint, cint, cchar *);
extern bool flip_data(byte *data, cint np, cint, cint, cint, cchar *dir);

/* subsample.c */
extern void subsample_byte(byte *, cbyte *,
	int, int, int, int, int, int, int, int);
extern void subsample_short(short *, cshort *,
	int, int, int, int, int, int, int, int);
extern void subsample_int(int *, cint *,
	int, int, int, int, int, int, int, int);
extern void subsample_float(float *, cfloat *,
	int, int, int, int, int, int, int, int);
extern void subsample_double(double *, cdouble *,
	int, int, int, int, int, int, int, int);

/* vol2mat.c */
extern void vol2mat(byte *, cbyte *, cint, cint, int, cint, cint, cint);
extern void vote_byte(byte *, cbyte *, int, int, int, int, int, int, int);

#endif	/* DefMisc */
