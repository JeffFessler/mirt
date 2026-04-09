// helix,fbp,def.h
// Copyright 2010-07-26, Jeff Fessler, University of Michigan

#ifndef jf_helix_fbp_def_h
#define jf_helix_fbp_def_h

#include "defs-env.h"

extern sof helix_fbp_ssrb(
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
cbyte *mask2, // [nx ny]
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

#endif // jf_helix_fbp_def_h
