// jf,mex,def.h
// Copyright 2005-12-6, Jeff Fessler, University of Michigan

#include "def,mexarg.h"

#define jf_mex_proto \
	int nlhs, mxArray *plhs[], int nrhs, Const mxArray *prhs[]

extern sof back2_mex(jf_mex_proto);
extern sof bspline_mex(jf_mex_proto);
extern sof cbct_mex(jf_mex_proto);
extern sof dtft_mex(jf_mex_proto);
extern sof fbp_fan_mex(jf_mex_proto);
extern sof fdk_mex(jf_mex_proto);
extern sof kde_mex(jf_mex_proto);
extern sof ptab_mex(jf_mex_proto);
extern sof moj2_mex(jf_mex_proto);
extern sof ncore_mex(jf_mex_proto);
extern sof jmh_sse_mex(jf_mex_proto);

#undef jf_mex_proto
