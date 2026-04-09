/*
* def,mexarg.h
* all purpose header file for my mex routines
*/

#if defined(Mmex)

#include "mex.h"

#ifndef NO_BOOL
#define NO_BOOL /* stupid Mathworks provides bool in matrix.h */
#endif
#undef CountAlloc
#include "defs-env.h"

#ifndef jf_clock
#define jf_clock_declare clock_t jf_clock_store;
#define jf_clock0	jf_clock_store = clock();
#define jf_clock	( (clock() - jf_clock_store) / (double) CLOCKS_PER_SEC )
#endif


/*
* in Matlab 2006b release, they changed many things
* to use mwSize and mwIndex, to support 64bit.
* I include these for backwards compat.
*/
#if defined(Need_mwIndex_int)
	typedef	int mwSize;
	typedef	int mwIndex;
#endif

// mexarg.c
typedef Const mxArray *Cmx;
extern char *mxu_string(Cmx mx, cchar *);
extern jool mxu_string_free(char *);
extern jool mxu_arg(cint, Const mxArray *[]);
extern int mxu_numel(Cmx mx);

#define mxGetPr_cint32(mx) ((cint *) mxGetPr(mx))
#define mxGetPr_cfloat(mx) ((cfloat *) mxGetPr(mx))
#define mxGetPr_cdouble(mx) ((cdouble *) mxGetPr(mx))

#define mxGetInt(mx) (*((cint *) mxGetData(mx)))
#define mxGetDouble(mx) (*((cdouble *) mxGetData(mx)))
#define mxGetSingle(mx) (*((cfloat *) mxGetData(mx)))

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

#endif /* Mmex */
