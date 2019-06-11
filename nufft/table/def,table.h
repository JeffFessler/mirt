// def,table.h
// header for table-based NUFFT interpolators

#include <stdio.h>
#include "mex.h"

#ifndef Fail
#define Fail(msg) { \
	(void)fprintf(stderr, "FAIL %s %d: %s\n", __FILE__, __LINE__, msg); \
	(void)fflush(stderr); \
	return 0; }
#endif

#ifndef Note1
#define Note1(msg, arg) { \
	(void) fprintf(stdout, "Note %s %d: ", __FILE__, __LINE__); \
	(void) fprintf(stdout, msg, arg); \
	(void) fprintf(stdout, "\n"); \
	(void) fflush(stdout); \
	}
#endif

#ifndef Call
#define Call(fun, arg)	{ if (!(fun arg)) Fail(#fun) }
#endif

#define mexFail(str)	mexErrMsgTxt(str);
// {mexEvalString("sprintf('fail: " #str "')"); return; }

#define mxIsInt32n(mx, n) \
	( (n == mxGetM(mx) * mxGetN(mx)) && mxIsInt32(mx) )
#define mxIsScalarInt32(mx) mxIsInt32n(mx, 1)

#define mxIsComplexSingle(mx) \
	(mxIsSingle(mx) && mxIsComplex(mx))
#define mxIsRealSingle(mx) \
	(mxIsSingle(mx) && !mxIsComplex(mx))
#define mxIsComplexDouble(mx) \
	(mxIsDouble(mx) && mxIsComplex(mx))
#define mxIsRealDouble(mx) \
	(mxIsDouble(mx) && !mxIsComplex(mx))
#define mxIsScalarSingle(mx) \
	( (1 == mxGetM(mx)) && (1 == mxGetN(mx)) && mxIsRealSingle(mx) )
#define mxIsScalarDouble(mx) \
	( (1 == mxGetM(mx)) && (1 == mxGetN(mx)) && mxIsRealDouble(mx) )
