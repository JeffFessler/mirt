// mexarg.c
// Show matlab mex arguments (for debugging and error messages)
// Copyright 2001-04-23, Jeff Fessler, University of Michigan

#include "def,mexarg.h"

#ifdef Mmex // needed for wt.c -> wtfmex.c

// mxu_streq()
// check equality
sof mxu_streq(Const mxArray *mx, cchar *arg, cchar *pattern)
{
	char *string;
	Call(string = mxu_string, (mx, arg))
	sof out = Streq(string, pattern);
	Call(mxu_string_free, (string))
	return out;
}


// mxu_string()
// caller must free using mxu_string_free()
char *mxu_string(Const mxArray *mx, cchar *arg)
{
	char *string;
	int n = mxGetM(mx) * mxGetN(mx) + 1;

	if (!mxIsChar(mx))
		Fail1("%s must be char array", arg)

#if 1
	Call(string = (char *) mxCalloc, (n, sizeof(char)))
#else
	Mem0(string, n)
#endif
	if (mxGetString(mx, string, n))
		Warn("bug with mxGetString")
	return string;
}


// mxu_string_free()
sof mxu_string_free(char *s)
{
#if 1
	mxFree(s);
#else
	Free0(s)
#endif
	Ok
}


// mxu_showdim()
static sof mxu_showdim(Const mxArray *mx)
{
	int id, ndim;
	Const mwSize *dims;

	Call(ndim = mxGetNumberOfDimensions, (mx))
	Call(dims = mxGetDimensions, (mx))

	printf("dims");
	for (id=0; id < ndim; ++id)
		printf(" %d", (int) dims[id]);
#if 0
	if (ndim == 1)	Note1("dims %d", dims[0])
	if (ndim == 2)	Note2("dims %d %d", dims[0], dims[1])
	if (ndim == 3)	Note3("dims %d %d %d", dims[0], dims[1], dims[2])
	if (ndim == 4)	Note4("dims %d %d %d %d", dims[0], dims[1], dims[2], dims[3])
	if (ndim > 4)	Note4("dims %d %d %d ... %d", dims[0], dims[1], dims[2], dims[ndim-1])
#endif
	Ok
}


// mxu_arg()
// Show arguments
sof mxu_arg(cint nmx, Const mxArray *pmx[])
{
	Note1("narg=%d", nmx)

	for (int ii=0; ii < nmx; ++ii)
	{
		Const mxArray *mx = pmx[ii];

		printf("arg %d, ", ii);
		Call(mxu_showdim, (mx))
		printf(", ");

		if (mxIsChar(mx))
		{
			char *arg;
			Call(arg = mxu_string, (mx, ""))
			printf("char, '%s'", arg);
			Call(mxu_string_free, (arg))
		}

		else if (mxIsUint8(mx))
		{
			byte val = *((cbyte *) mxGetData(mx));
			printf("uint8, val[0] = %d", (int) val);
		}

		else if (mxIsInt32(mx))
		{
			int val = mxGetInt(mx);
			printf("int32, val[0] = %d", val);
		}

		else if (mxIsSingle(mx))
		{
			float val = *((cfloat *) mxGetData(mx));
			printf("single, val[0] = %g", val);
		}

		else if (mxIsDouble(mx))
		{
			double val = *((cdouble *) mxGetData(mx));
			printf("double, val[0] = %g", val);
		}

		else if (mxIsCell(mx))
			printf("cell");

		else if (mxIsStruct(mx))
			printf("struct");

		else
			printf("UNKNOWN!?");

		printf("\n");
	}

	Ok
}


// mxu_numel()
// # of elements in a matlab array
int mxu_numel(Cmx mx)
{
	int ndim;
	Const mwSize *dims;

	Call(ndim = mxGetNumberOfDimensions, (mx))
	Call(dims = mxGetDimensions, (mx))

	int numel = 1;
	int ii; // for cuda
	for (ii=0; ii < ndim; ++ii)
		numel *= dims[ii];

	return numel;
}


#endif // Mmex
