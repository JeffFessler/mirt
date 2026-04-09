/*
* defs-misc.h
*
* Copyright 1993, Jeff Fessler, The University of Michigan
*/
#ifndef DefMisc
#define DefMisc

#include "defs-env.h"

/* penalty,diff.c */
#if 0
extern sof penalty_diff_wk_tight_offset(float *, cfloat *, cint, cint, cint, cint, cdouble);
extern sof penalty_diff_wk_leak_offset(float *, cfloat *, cint, cint, cint, cint, cdouble);
extern sof penalty_diff1_wk_tight(float *, cfloat *, cint *, cint, cint *, cint, cdouble);
extern sof penalty_diff1_wk_leak(float *, cfloat *, cint *, cint, cint *, cint, cdouble);
extern sof penalty_diff2_wk_tight(float *, cfloat *, cint *, cint, cint *, cint, cdouble);
extern sof penalty_diff2_wk_leak(float *, cfloat *, cint *, cint, cint *, cint, cdouble);
#endif
extern sof penalty_diff_wk_tight(float *, cfloat *, cint *, cint, cint *, cint, cdouble, cint);
extern sof penalty_diff_wk_leak(float *, cfloat *, cint *, cint, cint *, cint, cdouble, cint);
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
extern sof penalty_diff_forw(float *, cfloat *, cint, cint *, cint, cint, cint, cint);
extern sof penalty_diff_back(float *, cfloat *, cint, cint *, cint, cint, cint, cint);

#endif	/* DefMisc */
