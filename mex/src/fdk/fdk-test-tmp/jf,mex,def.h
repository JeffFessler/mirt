/*
* jf,mex,def.h
*
* Copyright 2005-12-6, Jeff Fessler, University of Michigan
*/
#include "def,mexarg.h"

#define jf_mex_proto \
	int nlhs, mxArray *plhs[], int nrhs, Const mxArray *prhs[]

extern jool back2_mex(jf_mex_proto);
extern jool bspline_mex(jf_mex_proto);
extern jool cbct_mex(jf_mex_proto);
extern jool dtft_mex(jf_mex_proto);
extern jool fbp_fan_mex(jf_mex_proto);
extern jool fdk_mex(jf_mex_proto);
extern jool ptab_mex(jf_mex_proto);
extern jool moj2_mex(jf_mex_proto);
extern jool jmh_sse_mex(jf_mex_proto);

#undef jf_mex_proto
