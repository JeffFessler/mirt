/*
* def,fdk.h
* Copyright 2005-6-27, Jeff Fessler, University of Michigan
*/
#include "cbct,def.h"

extern void fdk_st_help(void);
extern void fdk_ts_help(void);

extern void fdk_st_back1(
float *image,
cint	nx,
cint	ny,
cint	nz,
cdouble	dx,
cdouble	dy,
cdouble	dz,
cdouble	offset_x,
cdouble	offset_y,
cdouble	offset_z,
cbyte	*mask,
cdouble	dso,
cdouble	dsd,
cint	ns,
cint	nt,
cdouble	ds,
cdouble	dt,
cdouble	offset_s,
cdouble	offset_t,
cfloat	*proj,
cdouble beta,
cjool is_flat);


extern void fdk_ts_back1(
float *image,
cint	nx,
cint	ny,
cint	nz,
cdouble	dx,
cdouble	dy,
cdouble	dz,
cdouble	offset_x,
cdouble	offset_y,
cdouble	offset_z,
cbyte	*mask,
cbyte	mask_id,
cdouble	dso,
cdouble	dsd,
cdouble	dfs,
cint	ns,
cint	nt,
cdouble	ds,
cdouble	dt,
cdouble	offset_s,
cdouble	offset_t,
cfloat	*proj,
cdouble beta);


// fdk,ts,t.c
extern jool fdk_ts_back_t(
float *image,
const cbct_ig *ig,
const cbct_cg *cg,
cint na,
cfloat  *proj,
cdouble *beta,
cint nthread,
cint chat);
