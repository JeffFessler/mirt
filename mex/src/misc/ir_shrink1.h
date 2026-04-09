// ir_shrink1.h
// Copyright 2012-06-25, Jeff Fessler, University of Michigan

#ifndef ir_shrink1_h
#define ir_shrink1_h

#include "defs-env.h"


// ir_shrink1_p()
extern sof ir_shrink1_p(
float *x, // [N]
cfloat *y, // [N]
cint N,
cfloat *reg, // [N]
cfloat *table, // [K]
cfloat dt,
cint K,
cfloat thresh,
cint niter,
cint nthread,
cint chat);

#endif // ir_shrink1_h
