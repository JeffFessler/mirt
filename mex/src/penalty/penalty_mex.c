// penalty_mex.c
// Matlab mex gateway routine for roughness penalty differencing routines
// Copyright 2004-5-14, Jeff Fessler, University of Michigan
#include "def,mexarg.h"
#include "defs-misc.h"

#define Usage	"Usage: penalty_mex(?)\n"
#define Max_ndim	8	// maximum # of input dimensions

#if defined(Mmex)

static void penalty_mex_help()
{
	printf("Usage for penalty_mex:\n\
	d = function('diff1,forw1', x, offsets, [lastdim]);\n\
		d = C * x\n\
		offsets is int32 array such as [1 nx nx-1 nx+1].\n\
		Usually lastdim is absent, and x is treated as single 'image',\n\
		and output d is of size [size(x) length(offsets)].\n\
		If lastdim is int32, then it specifies which dimension of x\n\
		is the last one associated with one 'image data'; other dims\n\
		that follow are further 'realizations', and output d has size\n\
		[size(x)(1:lastdim) length(offsets) size(s)(lastdim+1:end)].\n\
\n\
	d = function('diff2,forw1', x, offsets, [lastdim]);\n\
		like 'diff1,forw1' but with 2nd differences.\n\
\n\
	d = function('diff1,forw2', ...);\n\
		d = C.^2 * x\n\
		like 'diff1,forw1' but with squared c_kj values.\n\
\n\
	x = function('diff1,back1', d, offsets, [lastdim]);\n\
		offsets is int32 array such as [-1 -nx -nx-1 -nx+1]\n\
		x = C' * d\n\
		output x is of size [length(d(:))/length(offsets)]\n\
\n\
	d = function('diff2,back1', x, offsets, [lastdim]);\n\
		like 'diff1,back1' but with 2nd differences.\n\
\n\
	x = function('diff1,back2', wk, offsets, [lastdim]);\n\
		like diff1,back1 but with squared c_kj values.\n\
		useful for diagonal of Hessian of penalty.\n\
		x = C' .^2 * wk \n\
\n\
	d = function('wk,tight,1', kappa, offsets, distance_power);\n\
		kappa is floats, often just a binary support mask in floats.\n\
		output is [size(kappa) length(offsets)]\n\
		contains weights {w_k} for each set of differences\n\
		where both neighbors are within the mask.\n\
			wk = kappa_j * kappa_k / distance^distance_power\n\
		Use distance_power=1 for the classical choice, leading\n\
		to 1/dis |a-b|^2, where dis=sqrt(2) for diagonal neighbors.\n\
		Use distance_power=2 for the possible improved choice,\n\
		leading to |(a-b)/dis|^2.\n\
	d = function('wk,leak,1', kappa, offsets, distance_power);\n\
		same but the 'leak' version that includes mask neighbors.\n\
	Use wk,*,2 for 2nd-order differences.\n\
	\n");
}


//
// penalty_wk_mex()
// function('wk,tight,?', kappa, offsets, distance_power)
// function( 'wk,leak,?', kappa, offsets, distance_power)
//
static sof penalty_wk_mex(
int nlhs, mxArray *plhs[],
int nrhs, Const mxArray *prhs[])
{
	char *arg;
	Const mxArray *mx_kappa;
	Const mxArray *mx_offsets;
	Const mxArray *mx_power;
	int n_offset;
	int nod_i;
	int nod_o;
	cint *dim_i;
	int dim_o[Max_ndim+1];

	if (nlhs != 1 || nrhs != 4) {
		penalty_mex_help();
		Call(mxu_arg, (nrhs, prhs))
		Fail(Usage)
	}

	mx_kappa = prhs[1];
	mx_offsets = prhs[2];
	mx_power = prhs[3];
	if (!mxIsRealSingle(mx_kappa)) Fail("mask must be real single")
	if (!mxIsInt32(mx_offsets)) Fail("offsets must be int32")
	if (!mxIsScalarDouble(mx_power)) Fail("power must be real scalar double")

	nod_i = mxGetNumberOfDimensions(mx_kappa);
	Call(dim_i = mxGetDimensions, (mx_kappa))

	n_offset = mxGetM(mx_offsets) * mxGetN(mx_offsets);
	nod_o = nod_i + 1;
	if (nod_i > Max_ndim) Fail("ndim limit exceeded!?")
	Bcopy(dim_i, dim_o, nod_i)
	dim_o[nod_o-1] = n_offset;

	Call(plhs[0] = mxCreateNumericArray,
		(nod_o, dim_o, mxSINGLE_CLASS, mxREAL))

	Call(arg = mxu_string, (prhs[0], "1st argument"))

	if (Streqn(arg, "wk,tight", 8)) {
		int order;
		Call(1 == sscanf, (arg, "wk,tight,%d", &order))
		Call(penalty_diff_wk_tight, (
			(float *) mxGetData(plhs[0]),
			(cfloat *) mxGetData(mx_kappa),
			dim_i, nod_i,
			(cint *) mxGetData(mx_offsets), n_offset,
			*mxGetPr(mx_power), order))
	}

	else if (Streqn(arg, "wk,leak", 7)) {
		int order;
		Call(1 == sscanf, (arg, "wk,leak,%d", &order))
		Call(penalty_diff_wk_leak, (
			(float *) mxGetData(plhs[0]),
			(cfloat *) mxGetData(mx_kappa),
			dim_i, nod_i,
			(cint *) mxGetData(mx_offsets), n_offset,
			*mxGetPr(mx_power), order))
	}

	else
		Fail1("unknown arg %s", arg)

	Call(mxu_string_free, (arg))

	Ok
}


/*
* penalty_diff_forw_mex()
* function('diff1,forw1', x, offsets, [lastdim])
* also diff1,forw2 and diff2,forw1
*/
static sof penalty_diff_forw_mex(
int nlhs, mxArray *plhs[],
int nrhs, Const mxArray *prhs[],
cint order, // 2 for 2nd differences
cint power)
{
	Const mxArray *mx_x;
	Const mxArray *mx_offsets;
	int ii;
	int nn = 1; /* image size */
	int nr = 1; /* # of realizations */
	int n_offset;
	int nod_i;
	int nod_o;
	cint *dim_i;
	int dim_o[Max_ndim+1];

	if (nlhs != 1 || nrhs < 3 || nrhs > 4) {
		penalty_mex_help();
		Call(mxu_arg, (nrhs, prhs))
		Fail(Usage)
	}

	mx_x = prhs[1];
	mx_offsets = prhs[2];
	if (!mxIsRealSingle(mx_x)) Fail("x must be real single")
	if (!mxIsInt32(mx_offsets)) Fail("offsets must be int32")

	nod_i = mxGetNumberOfDimensions(mx_x);
	Call(dim_i = mxGetDimensions, (mx_x))

	if (nod_i > Max_ndim) Fail("ndim limit exceeded!?")

	n_offset = mxGetM(mx_offsets) * mxGetN(mx_offsets);
	nod_o = nod_i + 1;

	/*
	* here x is a single realization (no "lastdim" arg)
	*/
	if (nrhs == 3) {
		Bcopy(dim_i, dim_o, nod_i)
		dim_o[nod_o-1] = n_offset;

		for (ii=0; ii < nod_i; ++ii)
			nn *= dim_i[ii];
	}

	/*
	* here x is [n1,...,nl, nl+1,...,nd]
	* output d is [n1,...,nl, n_offset, nl+1,...,nd]
	*/
	else {
		Const mxArray *mx_lastdim = prhs[3];
		int lastdim;

		if (!mxIsInt32(mx_lastdim))
			Fail("lastdim must be int32")
		if (2 != mxGetNumberOfDimensions(mx_lastdim)
			|| 1 != mxGetM(mx_lastdim) * mxGetN(mx_lastdim))
			Fail("lastdim must be [1 1]")

		lastdim = *( (cint *) mxGetData(mx_lastdim) );
		if (lastdim < 1 || lastdim > nod_i)
			Fail1("bad lastdim %d", lastdim)

		Bcopy(dim_i, dim_o, lastdim)
		dim_o[lastdim] = n_offset;
		Bcopy(dim_i+lastdim, dim_o+lastdim+1, nod_i-lastdim)

		for (ii=0; ii < lastdim; ++ii)
			nn *= dim_i[ii];

		for (ii=lastdim; ii < nod_i; ++ii)
			nr *= dim_i[ii];
	}

	Call(plhs[0] = mxCreateNumericArray,
		(nod_o, dim_o, mxSINGLE_CLASS, mxREAL))

	Call(penalty_diff_forw, (
		(float *) mxGetData(plhs[0]),
		(cfloat *) mxGetData(mx_x),
		nn,
		(cint *) mxGetData(mx_offsets), n_offset, nr,
		order, power))

	Ok
}


/*
* penalty_diff_back_mex()
* function('diff1,back1', d, offsets, [lastdim])
* also diff1,back2 and diff2,back1
*/
static sof penalty_diff_back_mex(
int nlhs, mxArray *plhs[],
int nrhs, Const mxArray *prhs[],
cint order,
cint power)
{
	Const mxArray *mx_d;
	Const mxArray *mx_offsets;
	int ii;
	int nn = 1; /* image size */
	int nr = 1; /* # of realizations */
	int n_offset;
	int nod_i;
	cint *dim_i;
	int nod_o;
	int dim_o[Max_ndim];

	if (nlhs != 1 || nrhs < 3 || nrhs > 4) {
		penalty_mex_help();
		Call(mxu_arg, (nrhs, prhs))
		Fail(Usage)
	}

	mx_d = prhs[1];
	mx_offsets = prhs[2];
	if (!mxIsRealSingle(mx_d)) Fail("d must be single")
	if (!mxIsInt32(mx_offsets)) Fail("offsets must be int32")

	Call(dim_i = mxGetDimensions, (mx_d))
	nod_i = mxGetNumberOfDimensions(mx_d);

	n_offset = mxGetM(mx_offsets) * mxGetN(mx_offsets);

	/* nod_i-1 would work since nod_i >= 2 says matlab, but be safe: */
	nod_o = Max(nod_i - 1, 1);
	if (nod_o > Max_ndim) Fail("ndim limit exceeded!?")

	/* single realization */
	if (nrhs == 3) {
		if (n_offset == 1)
			nod_o = nod_i; /* trick: since matlab cuts final '1' */
		else if (dim_i[nod_i-1] != n_offset) {
			for (ii=0; ii < nod_i; ++ii)
				Warn2("dim input [%d] = %d", ii, dim_i[ii])
			Fail2("last input dim is %d != length(offsets)=%d", 
				dim_i[nod_i-1], n_offset)
		}

		Bcopy(dim_i, dim_o, nod_o)

		for (ii=0; ii < nod_o; ++ii)
			nn *= dim_o[ii];
	}

	/*
	* multiple realizations (possibly)
	* input d is [n1,...,nl, n_offset, n_{l+1},...,nd]
	* output x is [n1,...,nl, n_{l+1},...,nd]
	* so nn = n1 * ... * nl and nr = n_{l+1} * ... * nd.
	* When n_offset = 1, d could be simply [n1,...,nl], not [n1,...,nl,1].
	*/
	else {
		Cmx mx_lastdim = prhs[3];
		int lastdim;

		Call(mxIsScalarInt32, (mx_lastdim))
		lastdim = *( (cint *) mxGetData(mx_lastdim) );
		if (lastdim < 1) Fail1("bad lastdim=%d < 1", lastdim)

		if (n_offset == 1 && lastdim >= nod_i) {
			nod_o = Min(nod_i, lastdim);
			Bcopy(dim_i, dim_o, Min(nod_i, lastdim))
		}
		else {
			if (lastdim > nod_o)
				Fail3("bad lastdim=%d vs nod_o=%d, nod_i=%d",
					lastdim, nod_o, nod_i)

			if (dim_i[lastdim] != n_offset)
				Fail3("size(d,%d)=%d != length(offset)=%d",
					lastdim, dim_i[lastdim], n_offset)

			Bcopy(dim_i, dim_o, lastdim)
			Bcopy(dim_i+lastdim+1, dim_o+lastdim, nod_i-lastdim)

			for (ii=lastdim+1; ii < nod_i; ++ii)
				nr *= dim_i[ii];
		}

		for (ii=0; ii < Min(nod_i, lastdim); ++ii)
			nn *= dim_i[ii];
	}

	Call(plhs[0] = mxCreateNumericArray,
		(nod_o, dim_o, mxSINGLE_CLASS, mxREAL))

	Call(penalty_diff_back, (
		(float *) mxGetData(plhs[0]),
		(cfloat *) mxGetData(mx_d),
		nn,
		(cint *) mxGetData(mx_offsets), n_offset, nr,
		order, power))

	Ok
}


/* intermediate gateway routine */
static sof penalty_mex(
int nlhs, mxArray *plhs[],
int nrhs, Const mxArray *prhs[])
{
	char *arg;

	if (nrhs < 1 || !mxIsChar(prhs[0])) {
		penalty_mex_help();
		Call(mxu_arg, (nrhs, prhs))
		Fail(Usage)
	}

	/* 'help' */
	if (nrhs == 1 && mxIsChar(prhs[0])) {
		penalty_mex_help();
		Ok
	}

	Call(arg = mxu_string, (prhs[0], "1st argument"))

	/*
	* diff forward: C * x
	*/
	if (Streqn(arg, "diff1,forw", 10))
		Call(penalty_diff_forw_mex,
		(nlhs, plhs, nrhs, prhs, 1, 1+Streq(arg, "diff1,forw2")))

	else if (Streqn(arg, "diff2,forw", 10))
		Call(penalty_diff_forw_mex,
		(nlhs, plhs, nrhs, prhs, 2, 1+Streq(arg, "diff2,forw2")))

	/*
	* diff back: C' * d
	*/
	else if (Streqn(arg, "diff1,back", 10))
		Call(penalty_diff_back_mex,
		(nlhs, plhs, nrhs, prhs, 1, 1+Streq(arg, "diff1,back2")))

	else if (Streqn(arg, "diff2,back", 10))
		Call(penalty_diff_back_mex,
		(nlhs, plhs, nrhs, prhs, 2, 1+Streq(arg, "diff2,back2")))

	/*
	* wk factors for \sum_k w_k \pot([Cx]_k)
	*/
	else if (Streqn(arg, "wk,", 3))
		Call(penalty_wk_mex, (nlhs, plhs, nrhs, prhs))

	else
		Fail1("command '%s' unknown", arg)

	Call(mxu_string_free, (arg))
	Ok
}


/* gateway routine */
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, Const mxArray *prhs[])
{
	if (!nlhs && !nrhs) {
		penalty_mex_help();
		return;
	}
	if (!penalty_mex(nlhs, plhs, nrhs, prhs))
		mexErrMsgTxt("penalty_mex");
}

#endif
