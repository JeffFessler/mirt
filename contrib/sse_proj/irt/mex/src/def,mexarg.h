/*
* def,mexarg.h
*/
#define NO_BOOL         /* stupid Mathworks provides bool in matrix.h */
#undef CountAlloc
#ifndef bool
# define bool int
# include "defs-env.h"
# undef bool
#else
# include "defs-env.h"
#endif

#if defined(Mmex)
#include "mex.h"
typedef Const mxArray *Cmx;
extern char *mx_string(Const mxArray *mx, cchar *);
extern bool mx_string_free(char *);
extern bool mexarg(cint, Const mxArray *[]);

#define mxGetInt(mx) (*((cint *) mxGetData(mx)))

#define mxIsScalar(mx) \
	( (2 == mxGetNumberOfDimensions(mx)) \
		&& (1 == mxGetM(mx)) && (1 == mxGetN(mx)) )
#define mxIsScalarInt32(mx) \
	( mxIsScalar(mx) && mxIsInt32(mx) )
#define mxIsComplexSingle(mx) \
	(mxIsSingle(mx) && mxIsComplex(mx))
#define mxIsComplexDouble(mx) \
	(mxIsDouble(mx) && mxIsComplex(mx))
#define mxIsRealSingle(mx) \
	(mxIsSingle(mx) && !mxIsComplex(mx))
#define mxIsRealDouble(mx) \
	(mxIsDouble(mx) && !mxIsComplex(mx))
#define mxIsScalarSingle(mx) \
	( mxIsScalar(mx) && mxIsRealSingle(mx) )
#define mxIsScalarDouble(mx) \
	( mxIsScalar(mx) && mxIsRealDouble(mx) )
#endif
