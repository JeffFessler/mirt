/*
* mexarg.c
* Show matlab mex arguments (for debugging and error messages)
*
* Copyright 01-04-23, Jeff Fessler, The University of Michigan
*/
#include "def,mexarg.h"

#if defined(Mmex)

/*
* mx_string()
* caller must free using mx_string_free()
*/
char *mx_string(Const mxArray *mx, cchar *arg)
{
	char	*string;
	int	n = mxGetM(mx) * mxGetN(mx) + 1;

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

bool mx_string_free(char *s)
{
#if 1
	mxFree(s);
#else
	Free0(s)
#endif
	Ok
}

static bool mex_showdim(Const mxArray *mx)
{
	int ndim;
	cint *dims;
	Call(ndim = mxGetNumberOfDimensions, (mx))
	Call(dims = mxGetDimensions, (mx))
	if (ndim == 1)	Note1("dims %d", dims[0])
	if (ndim == 2)	Note2("dims %d %d", dims[0], dims[1])
	if (ndim == 3)	Note3("dims %d %d %d", dims[0], dims[1], dims[2])
	if (ndim == 4)	Note4("dims %d %d %d %d", dims[0], dims[1], dims[2], dims[3])
	if (ndim > 4)	Note4("dims %d %d %d ... %d", dims[0], dims[1], dims[2], dims[ndim-1])
	Ok
}


/*
* mexarg()
* Show arguments
*/
bool mexarg(cint nmx, Const mxArray *pmx[])
{
	int ii;

	Note1("narg=%d", nmx)

	for (ii=0; ii < nmx; ++ii) {
		Const mxArray *mx = pmx[ii];

		if (mxIsChar(mx)) {
			char *arg;
			Call(arg = mx_string, (mx, ""))
			Note2("arg %d, char, '%s'", ii, arg)
			Call(mx_string_free, (arg))
		}

		else if (mxIsInt32(mx)) {
			int val = *((int *) mxGetData(mx));
			Note2("arg %d, int32, val[0]=%d", ii, val)
		}

		else if (mxIsSingle(mx)) {
			float val = *((cfloat *) mxGetData(mx));
			Note2("arg %d, single, val[0]=%g", ii, val)
		}

		else if (mxIsDouble(mx)) {
			double val = *((cdouble *) mxGetData(mx));
			Note2("arg %d, double, val[0]=%g", ii, val)
		}

		else
			Note1("arg %d UNKNOWN!?", ii)

		Call(mex_showdim, (mx))
	}

	Ok
}

#endif
