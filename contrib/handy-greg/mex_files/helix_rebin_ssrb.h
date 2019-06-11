/*
 *  helix_rebin_ssrb.h
 *  
 *
 *  Created by Gregory Handy on 7/15/10.
 *  Copyright 2010 UMBC. All rights reserved.
 *
 */

#include "defs-env.h"

extern sof helix_rebin_ssrb(
float *sino,
float *z_orbit_start,
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
