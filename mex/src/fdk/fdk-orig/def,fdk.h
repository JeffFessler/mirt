// def,fdk.h
// Copyright 2005-6-27, Jeff Fessler, University of Michigan

#include "cbct,def.h"

extern void fdk_st_help(void);
extern void fdk_ts_help(void);

extern void fdk_st_back1(
float *image,
cint nx,
cint ny,
cint nz,
cfloat dx,
cfloat dy,
cfloat dz,
cfloat offset_x,
cfloat offset_y,
cfloat offset_z,
cbyte *mask,
cfloat dso,
cfloat dsd,
cfloat dfs,
cint ns,
cint nt,
cfloat ds,
cfloat dt,
cfloat offset_s,
cfloat offset_t,
cfloat *proj,
cfloat beta);


extern sof fdk_ts_back1(
float *image,
cint nx,
cint ny,
cint nz,
cfloat dx,
cfloat dy,
cfloat dz,
cfloat offset_x,
cfloat offset_y,
cfloat offset_z,
cbyte *mask,
cbyte mask_id,
cfloat dso,
cfloat dsd,
cfloat dfs,
cint ns,
cint nt,
cfloat ds,
cfloat dt,
cfloat offset_s,
cfloat offset_t,
cfloat beta,
cfloat *view_work);


// fdk,ts,t.c

extern void fdk_ts_put_view(
float *view_work, // [(nt+2) (ns+2)], assume outer border already 0
cfloat *proj, // [nt ns]
cint ns,
cint nt);

extern sof fdk_ts_back_t(
float *image,
const cbct_ig *ig,
const cbct_cg *cg,
cint na,
cfloat *proj,
cdouble *beta,
cint nthread,
cint chat);
