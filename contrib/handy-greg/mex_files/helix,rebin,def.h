// helix,rebin,def.h
// Copyright 2010-07-26, Jeff Fessler, University of Michigan

#ifndef jf_helix_ssrb_def_h
#define jf_helix_ssrb_def_h

#include "defs-env.h"

extern sof helix_rebin_ssrb(
float *sino, // [ns na_rebin nz]
float *z_orbit_start, // [num_betas]
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
