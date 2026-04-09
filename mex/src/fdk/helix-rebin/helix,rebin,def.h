// helix,rebin,def.h
// Copyright 2010-07-26, Jeff Fessler, University of Michigan

#ifndef jf_helix_ssrb_def_h
#define jf_helix_ssrb_def_h

#include "defs-env.h"

extern sof helix_rebin_dim(
int *p_na1,
float *p_orbit_short,
cfloat dso,
cfloat dsd,
cfloat dfs,
cfloat ds,
cfloat offset_s,
cfloat orbit,
cint ns,
cint na);

extern sof helix_rebin_ssrb_z(
float *sino, // [ns na_rebin]
float *p_z_orbit_start,
cint na_rebin,
cfloat zmid,
cfloat dso,
cfloat dsd,
cfloat dfs,
cint ns,
cint nt,
cint na,
cfloat ds, 
cfloat dt, 
cfloat offset_s, 
cfloat offset_t, 
cfloat pitch,
cfloat source_z0, // first z position of source
cfloat orbit, 
cfloat orbit_start, 
cfloat *proj); // [ns nt na] I: projection

extern sof helix_rebin_ssrb(
float *sino, // [ns na_rebin nz]
float *z_orbit_start, // [nz]
cint na_rebin,
cint nz,
cfloat dz,
cfloat offset_z,
cfloat dso,
cfloat dsd,
cfloat dfs,
cint ns,
cint nt,
cint na,
cfloat ds,
cfloat dt,
cfloat offset_s,
cfloat offset_t,
cfloat pitch,
cfloat source_z0,
cfloat orbit,
cfloat orbit_start,
cfloat *proj);

#endif // jf_helix_rebin_def_h
