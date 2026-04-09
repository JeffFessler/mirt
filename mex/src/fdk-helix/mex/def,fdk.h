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
cfloat source_zs, // DK
cfloat source_zs_min, // DK
cfloat source_zs_max, // DK
cint iz_min_beta, // DK
cint iz_max_beta, // DK
cint cone_par, // DK
cint w3d, // DK
cdouble source_dz_per_rad, // DK
cfloat pitch, // DK
cfloat *pow_tab, // DK
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
cdouble *source_zs, // DK
cint na1_half, // DK
cint cone_par, // DK
cint w3d, // DK
cdouble  source_dz_per_rad, // DK
cfloat pitch, // DK
cint nthread,
cint chat);



// fdk,w3d,pow,tab.c

extern sof fdk_w3d_pow_tab_init(
float *pow_tab,
cint nt,
cint nz,
cfloat dt,
cfloat dz,
cfloat dsd,
cint w3d,
cfloat pitch);

extern sof fdk_w3d_pow_tab(
float *ww,
cfloat *pow_tab,
cint nt,
cfloat t_bin,
cfloat t_bin_c,
cfloat wt,
cint iz,
cfloat wz,
cfloat pitch);
 
