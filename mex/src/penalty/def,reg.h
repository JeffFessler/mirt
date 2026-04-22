// def,reg.h
// Copyright 2006, Jeff Fessler, University of Michigan

#ifndef DefReg
#define DefReg

#include "defs-env.h"


// reg,dxyz

extern sof penalty_diff_dxyz(
int *p_dx, int *p_dy, int *p_dz,
cint nx, cint ny, cint nz,
cint offset);


// reg,misc

#if 0
typedef struct
{
	int K; // # of samples
	float dt; // t spacing
	cfloat *dk; // [1:K] samples of dpot (constant because user provided)
	float *ck; // [K+1] curvatures 0:K
	float *sk; // [K] cumulative sums for pot(t)
} reg_pot_table;
#endif

typedef struct
{
	float delta; // for optional use
	float param[3]; // [n_param <= 3] parameters
//	int n_param; // #
//	struct reg_pot_table *rpt; // used only for tabulated potentials
} reg_pot_data;

typedef Const reg_pot_data c_reg_pot_data;

// typedef float pot_t(cfloat x, cfloat *param);
typedef float pot_t(cfloat x, c_reg_pot_data *);

typedef struct
{
	pot_t *fpot; // potential function
	pot_t *dpot; // derivative
	pot_t *wpot; // weighting
	reg_pot_data rpd;
//	float param[3]; // [n_param <= 3] parameters
//	int n_param; // #
//	struct reg_pot_table *rpt; // used only for tabulated potentials
} reg_pot;

extern sof pot_set(
reg_pot *pot,
cfloat *pot_param,
cint n_pot_param,
cchar *pot_name);


// reg,cost.c
extern sof penalty_cost_offset(
double *p_value, // [1]
cint nn, // nx*ny*nz
cfloat *x, // [nn]
cfloat *kappa, // [nn]
cint *offset, // [n_offset]
cint n_offset,
cfloat *beta, // [n_offset] or [n_offset nn]
ctruf is_beta_array,
cchar *pot_name,
cfloat *pot_param, // [n_pot_param]
cint n_pot_param,
cint order, // 1st or 2nd order differences
cint control,
cint nthread,
cint chat);


// reg,diff.c
#if 0
extern sof penalty_diff_wk_tight_offset(float *, cfloat *, cint, cint, cint, cint, cdouble);
extern sof penalty_diff_wk_leak_offset(float *, cfloat *, cint, cint, cint, cint, cdouble);
extern sof penalty_diff1_wk_tight(float *, cfloat *, cint *, cint, cint *, cint, cdouble);
extern sof penalty_diff1_wk_leak(float *, cfloat *, cint *, cint, cint *, cint, cdouble);
extern sof penalty_diff2_wk_tight(float *, cfloat *, cint *, cint, cint *, cint, cdouble);
extern sof penalty_diff2_wk_leak(float *, cfloat *, cint *, cint, cint *, cint, cdouble);
#endif
extern sof penalty_diff_wk_tight
		(float *, cfloat *, cint *, cint, cint *, cint, cdouble, cint);
extern sof penalty_diff_wk_leak
		(float *, cfloat *, cint *, cint, cint *, cint, cdouble, cint);

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

extern sof penalty_diff_forw
		(float *, cfloat *, cint, cint *, cint, cint, cint, cint);
extern sof penalty_diff_back
		(float *, cfloat *, cint, cint *, cint, cint, cint, cint);



// reg,grad.c

extern sof penalty_cgrad_offset(
float *cgrad, // [nn] same size as image
cint nn, // nx*ny*nz
cfloat *x, // [nn]
cfloat *kappa, // [nn]
cint *offset, // [n_offset]
cint n_offset,
cfloat *beta, // [n_offset] or [n_offset nn]
ctruf is_beta_array,
cchar *pot_name,
cfloat *pot_param, // [n_pot_param]
cint n_pot_param,
cint order, // 1st or 2nd order differences
cint control,
cint nthread,
cint chat);


// reg,denom

extern sof penalty_denom_offset(
float *denom, // [nn] same size as image
cint nn, // nx*ny*nz
cfloat *x, // [nn]
cfloat *kappa, // [nn]
cint *offset, // [n_offset]
cint n_offset,
cfloat *beta, // [n_offset] or [n_offset nn]
ctruf is_beta_array,
cchar *pot_name,
cfloat *pot_param, // [n_pot_param]
cint n_pot_param,
cint order, // 1st or 2nd order differences
cint control,
cint nthread,
cint chat);

#if 0
extern sof penalty_cgrad(cint, cint,
	float *, cfloat *, cfloat *, cfloat *, cint *,
	char *, float *, cint, cint, cint);
extern sof penalty_denom(cint, cint,
	float *, cfloat *, cfloat *, cfloat *, cint *,
	char *, float *, cint, cint);
#endif
extern sof penalty_cgrad_denom(cint, cint,
	float *, cfloat *, cfloat *, cfloat *, cint *,
	cchar *, cfloat *, cint, cint);
extern sof penalty_cgrad_denom_par(cint, cint, cint,
	float *, cfloat *, cfloat *, cfloat *, cint *,
	cchar *, cfloat *, cint, cint);
extern sof penalty_cgden_par1(void *in, cint id, cint nthread);


// reg,mask2

typedef enum // z end conditions for regularizer
{
	reg_end_unknown,
	reg_end_replicate,
	reg_end_periodic
} reg_end_t;

typedef struct
{
	cbyte *mask2; // [nx ny] binary mask
	enum
	{
		reg_mask2_unknown, // unknown border condition
		reg_mask2_goes_to_edge, // mask goes to edge
		reg_mask2_border1, // verified to have 1 pixel border around
		reg_mask2_border2 // verified to have 2 pixel border around
	} type;
} reg_mask2;

extern sof reg_mask2_set(reg_mask2 *rm, cbyte *mask2, cint nx, cint ny);


// reg,pow

extern sof reg_set_distance_factor(
float *distance_factor, // [n_offset]
cint n1,
cint n2,
cint n3,
cint *offset, // [n_offset]
cint n_offset,
cdouble distance_power);


// reg,zxy,cost

extern sof reg_cost_zxy(
double *p_value, // [1] cost value of penalty function
cint nx,
cint ny,
cint nz,
cfloat *x, // [nz nx ny]
cfloat *kappa, // [nz nx ny]
cint *offset, // [n_offset]
cint n_offset,
cfloat *beta, // [n_offset]
cchar *pot_name,
cfloat *pot_param, // [n_pot_param]
cint n_pot_param,
Const reg_mask2 *rm,
cint order,
reg_end_t reg_end,
cint nthread,
cint chat);


// reg,zxy,cgrad

extern sof reg_zxy_d1(
int *d1, // [n_offset]
cint *offset, // [n_offset]
cint n_offset,
cint n1, // usually nz
cint n2, // usually nx
cint n3); // usually ny

extern sof reg_dis_max(
cint *offset, // [n_offset]
cint n_offset,
cint n1,
cint n2,
cint n3);

extern sof reg_cgrad_zxy(
float *cgrad, // [nz nx ny] gradient
cint nx,
cint ny,
cint nz,
cfloat *x, // [nz nx ny]
cfloat *kappa, // [nz nx ny]
cint *offset, // [n_offset]
cint n_offset,
cfloat *beta, // [n_offset]
cchar *pot_name,
cfloat *pot_param, // [n_pot_param]
cint n_pot_param,
Const reg_mask2 *rm,
cint order,
reg_end_t reg_end,
cint nthread,
cint chat);


// reg,zxy

extern sof reg_cgrad_denom_zxy(
float *cgrad, // [nz nx ny] gradient
float *denom, // [nz nx ny] denominator
cint nx,
cint ny,
cint nz,
cfloat *x, // [nz nx ny]
cfloat *kappa, // [nz nx ny]
cint *offset, // [n_offset]
cint n_offset,
cfloat *beta, // [n_offset]
cchar *pot_name,
cfloat *pot_param, // [n_pot_param]
cint n_pot_param,
Const reg_mask2 *rm,
cint order,
reg_end_t reg_end,
cint nthread,
cint chat);

#endif // DefReg
