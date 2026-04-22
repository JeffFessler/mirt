// penalty_mex.c
// Matlab mex gateway routine for roughness penalty differencing routines
//
// Copyright 2004-5-14, Jeff Fessler, University of Michigan

#include "def,mexarg.h"
#include "def,reg.h"
#include "jf,time.h"
#include "jf,thread1.h"

#define Usage	"Usage: penalty_mex(?)\n"
#define Max_ndim	8 // maximum # of input dimensions

#define FailUse { \
		penalty_mex_help(); \
		Call(mxu_arg, (nrhs, prhs)) \
		Fail(Usage) }

static void penalty_mex_help()
{
	printf("Usage for penalty_mex:\n\
		string: 'diff{1|2},{forw|back}{1|2|A}'\n\
	d = function('diff1,forw1', x, offsets, [lastdim]);\n\
		d = C * x\n\
		offsets is int32 array such as [1 nx nx-1 nx+1].\n\
		Usually lastdim is absent, so x is treated as single 'image',\n\
		and output d is of size [size(x) length(offsets)].\n\
		If lastdim is int32, then it specifies which dimension of x\n\
		is the last one associated with one 'image data'; other dims\n\
		that follow are further 'realizations', and output d has size\n\
		[size(x)(1:lastdim) length(offsets) size(s)(lastdim+1:end)].\n\
\n\
	d = function('diff1,forw2', ...);\n\
		d = C.^2 * x\n\
		like 'diff1,forw1' but with squared c_kj values.\n\
	d = function('diff1,forwA', ...);\n\
		trick: same as 'diff1,forw2' because |+/- 1| = (+/- 1)^2.\n\
\n\
	d = function('diff2,forw{1|2|A}', x, offsets, [lastdim]);\n\
		like 'diff1,forw{1|2|A}' but with 2nd differences.\n\
\n\
	x = function('diff1,back1', d, offsets, [lastdim]);\n\
		offsets is int32 array such as [-1 -nx -nx-1 -nx+1]\n\
		x = C' * d\n\
		output x is of size [length(d(:))/length(offsets)]\n\
\n\
	x = function('diff1,back2', wk, offsets, [lastdim]);\n\
		like diff1,back1 but with squared c_kj values.\n\
		useful for diagonal of Hessian of penalty.\n\
		x = C' .^2 * wk \n\
		same as 'diff1,backA' because |+/- 1| = (+/- 1)^2.\n\
\n\
	d = function('diff2,back{1|2}', x, offsets, [lastdim]);\n\
		like 'diff1,back{1|2}' but with 2nd differences.\n\
\n\
	d = function('diff2,backA', x, offsets, [lastdim]);\n\
		like 'diff2,back1' but with (1,2,1) not (-1,2,-1), for denom.\n\
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
\n\
	d = function('wk,leak,1', kappa, offsets, distance_power);\n\
		same but the 'leak' version that includes mask neighbors.\n\
\n\
	Use wk,*,2 for 2nd-order differences.\n\
\n\
	g = function('cgrad,offset',\n\
		x, kappa, offsets, beta_vector, pot_name, pot_param, \n\
			order, control, nthread)\n\
		gradient of R(x) = sum_k pot( [Cx]_k )\n\
			beta_vector can be [n_offset] or [n_offset nn]\n\
	d = function('denom,offset',\n\
		(same args as 'cgrad,offset' but provides SQS denominator)\n\
	R = function('penal,offset',\n\
		(same args as 'cgrad,offset' but provides penalty value)\n\
\n\
	g = function('cgrad,zxy',\n\
		x, kappa, offsets, beta_vector, pot_name, pot_param, \n\
			mask2, order, nthread)\n\
		gradient R(x) = sum_k pot( [Cx]_k )\n\
	[g w] = function('cgrad,denom,zxy',\n\
		x, kappa, offsets, beta_vector, pot_name, pot_param, \n\
			mask2, order, nthread)\n\
		gradient and denominator of R(x) = sum_k pot( [Cx]_k )\n\
\n\
	\n");
}


// penalty_wk_mex()
// function('wk,tight,?', kappa, offsets, distance_power)
// function( 'wk,leak,?', kappa, offsets, distance_power)
static sof penalty_wk_mex(
int nlhs, mxArray *plhs[],
int nrhs, Const mxArray *prhs[])
{
	char *arg;

	if (nlhs != 1 || nrhs != 4)
	{
		penalty_mex_help();
		Call(mxu_arg, (nrhs, prhs))
		Fail(Usage)
	}

	Const mxArray *mx_kappa = prhs[1];
	Const mxArray *mx_offsets = prhs[2];
	Const mxArray *mx_power = prhs[3];
	if (!mxIsRealSingle(mx_kappa)) Fail("mask must be real single")
	if (!mxIsInt32(mx_offsets)) Fail("offsets must be int32")
	if (!mxIsScalarDouble(mx_power)) Fail("power must be real scalar double")

	cSize_t nod_i = mxGetNumberOfDimensions(mx_kappa);
	const mwSize *dim_i;
	Call(dim_i = mxGetDimensions, (mx_kappa))

	cSize_t n_offset = mxGetM(mx_offsets) * mxGetN(mx_offsets);
	cSize_t nod_o = nod_i + 1;
	if (nod_i > Max_ndim) Fail("ndim limit exceeded!?")

	mwSize dim_o[Max_ndim+1];
	Bcopy(dim_i, dim_o, nod_i)
	dim_o[nod_o-1] = n_offset;

	Call(plhs[0] = mxCreateNumericArray,
		(nod_o, dim_o, mxSINGLE_CLASS, mxREAL))

	Call(arg = mxu_string, (prhs[0], "1st argument"))

	int int_dim_i[nod_i]; // trick=kludge
	for (size_t ii=0; ii < nod_i; ++ii)
		int_dim_i[ii] = (int) dim_i[ii];

	if (Streqn(arg, "wk,tight", 8))
	{
		int order;
		Call(1 == sscanf, (arg, "wk,tight,%d", &order))
		Call(penalty_diff_wk_tight, (
			(float *) mxGetData(plhs[0]),
			(cfloat *) mxGetData(mx_kappa),
			int_dim_i, nod_i,
			(cint *) mxGetData(mx_offsets), n_offset,
			*mxGetPr(mx_power), order))
	}

	else if (Streqn(arg, "wk,leak", 7))
	{
		int order;
		Call(1 == sscanf, (arg, "wk,leak,%d", &order))
		Call(penalty_diff_wk_leak, (
			(float *) mxGetData(plhs[0]),
			(cfloat *) mxGetData(mx_kappa),
			int_dim_i, nod_i,
			(cint *) mxGetData(mx_offsets), n_offset,
			*mxGetPr(mx_power), order))
	}

	else
		Fail1("unknown arg %s", arg)

	Call(mxu_string_free, (arg))

	Ok
}


// penalty_diff_forw_mex()
// function('diff1,forw1', x, offsets, [lastdim])
// also diff1,forw2 and diff2,forw1
static sof penalty_diff_forw_mex(
int nlhs, mxArray *plhs[],
int nrhs, Const mxArray *prhs[],
ctruf order, // 2 for 2nd differences
ctruf power)
{
	int nn = 1; // image size
	int nr = 1; // # of realizations

	if (nlhs != 1 || nrhs < 3 || nrhs > 4)
	{
		penalty_mex_help();
		Call(mxu_arg, (nrhs, prhs))
		Fail(Usage)
	}

	Cmx mx_x = prhs[1];
	Cmx mx_offsets = prhs[2];
	if (!mxIsRealSingle(mx_x)) Fail("x must be real single")
	if (!mxIsInt32(mx_offsets)) Fail("offsets must be int32")

	const mwSize *dim_i;
	Size_t nod_i = mxGetNumberOfDimensions(mx_x);
	Call(dim_i = mxGetDimensions, (mx_x))

	if (nod_i > Max_ndim) Fail("ndim limit exceeded!?")

	mwSize dim_o[Max_ndim+1];
	cSize_t n_offset = mxGetM(mx_offsets) * mxGetN(mx_offsets);
	cSize_t nod_o = nod_i + 1;

	// here x is a single realization (no "lastdim" arg)
	if (nrhs == 3)
	{
		Bcopy(dim_i, dim_o, nod_i)
		dim_o[nod_o-1] = n_offset;

		for (size_t ii=0; ii < nod_i; ++ii)
			nn *= dim_i[ii];
	}

	// here x is [n1,...,nl, nl+1,...,nd]
	// output d is [n1,...,nl, n_offset, nl+1,...,nd]
	else
	{
		Const mxArray *mx_lastdim = prhs[3];

		if (!mxIsInt32(mx_lastdim))
			Fail("lastdim must be int32")
		if (2 != mxGetNumberOfDimensions(mx_lastdim)
			|| 1 != mxGetM(mx_lastdim) * mxGetN(mx_lastdim))
			Fail("lastdim must be [1 1]")

		cSize_t lastdim = (cSize_t) *( (cint *) mxGetData(mx_lastdim) );
		if (lastdim < 1 || lastdim > nod_i)
			Fail1("bad lastdim %d", lastdim)

		Bcopy(dim_i, dim_o, lastdim)
		dim_o[lastdim] = n_offset;
		Bcopy(dim_i+lastdim, dim_o+lastdim+1, nod_i-lastdim)

		for (size_t ii=0; ii < lastdim; ++ii)
			nn *= dim_i[ii];

		for (size_t ii=lastdim; ii < nod_i; ++ii)
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


// penalty_diff_back_mex()
// function('diff1,back1', d, offsets, [lastdim])
// also diff1,back2 and diff2,back1
static sof penalty_diff_back_mex(
int nlhs, mxArray *plhs[],
int nrhs, Const mxArray *prhs[],
ctruf order,
ctruf power)
{
	Const mxArray *mx_d;
	Const mxArray *mx_offsets;
	int n_offset;
	int nod_i;
	const mwSize *dim_i;
	int nod_o;
	mwSize dim_o[Max_ndim];

	if (nlhs != 1 || nrhs < 3 || nrhs > 4)
	{
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

	// nod_i-1 would work since nod_i >= 2 says matlab, but be safe:
	nod_o = Max(nod_i - 1, 1);
	if (nod_o > Max_ndim) Fail("ndim limit exceeded!?")

	// single realization
	int nr = 1; // # of realizations
	int nn = 1; // image size
	if (nrhs == 3)
	{
		if (n_offset == 1)
			nod_o = nod_i; // trick: since matlab cuts final '1'
		else if ((int) dim_i[nod_i-1] != n_offset)
		{
			for (int ii=0; ii < nod_i; ++ii)
				Warn2("dim input [%d] = %d", ii, (int) dim_i[ii])
			Fail2("last input dim is %d != length(offsets)=%d",
				(int) dim_i[nod_i-1], n_offset)
		}

		Bcopy(dim_i, dim_o, nod_o)

		for (int ii=0; ii < nod_o; ++ii)
			nn *= dim_o[ii];
	}

	// multiple realizations (possibly)
	// input d is [n1,...,nl, n_offset, n_{l+1},...,nd]
	// output x is [n1,...,nl, n_{l+1},...,nd]
	// so nn = n1 * ... * nl and nr = n_{l+1} * ... * nd.
	// When n_offset = 1, d could be simply [n1,...,nl], not [n1,...,nl,1].
	else
	{
		Cmx mx_lastdim = prhs[3];

		Call(mxIsScalarInt32, (mx_lastdim))
		int lastdim = *( (cint *) mxGetData(mx_lastdim) );
		if (lastdim < 1) Fail1("bad lastdim=%d < 1", lastdim)

		if (n_offset == 1 && lastdim >= nod_i)
		{
			nod_o = Min(nod_i, lastdim);
			Bcopy(dim_i, dim_o, Min(nod_i, lastdim))
		}
		else
		{
			if (lastdim > nod_o)
				Fail3("bad lastdim=%d vs nod_o=%d, nod_i=%d",
					lastdim, nod_o, nod_i)

			if ((int) dim_i[lastdim] != n_offset)
				Fail3("size(d,%d)=%d != length(offset)=%d",
					lastdim, (int) dim_i[lastdim], n_offset)

			Bcopy(dim_i, dim_o, lastdim)
			Bcopy(dim_i+lastdim+1, dim_o+lastdim, nod_i-lastdim)

			for (int ii=lastdim+1; ii < nod_i; ++ii)
				nr *= dim_i[ii];
		}

		for (int ii=0; ii < Min(nod_i, lastdim); ++ii)
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


typedef
enum {p_error, p_cgrad, p_denom, p_penal} p_type_enum;


// penalty_cdp_offset_mex()
// function('{cgrad|denom|penal},offset', x, kappa,
//	offsets, beta_vector, pot_name, pot_param, order, control, nthread)
static sof penalty_cdp_offset_mex(
mxArray *plhs[],
Cmx mx_x, // [nn] single realization
Cmx mx_kappa, // [nn]
Cmx mx_offsets, // [n_offset]
Cmx mx_beta, // [n_offset]
Cmx mx_pot_name, // char
Cmx mx_pot_param, // [n_pot_param]
Cmx mx_order, // 1|2 int
Cmx mx_control,
Cmx mx_nthread, // int
p_type_enum p_type)
{
	char *pot_name;
	uint n_offset;
	const mwSize *dim_i;
	cint chat = 0;
	truf is_beta_array;

	Call(mxIsRealSingle, (mx_x))
	Call(mxIsRealSingle, (mx_kappa))
	Call(mxIsInt32, (mx_offsets))
	Call(mxIsRealSingle, (mx_beta))
	Call(mxIsChar, (mx_pot_name))
	Call(mxIsRealSingle, (mx_pot_param))
	Call(mxIsScalarInt32, (mx_order))
	Call(mxIsScalarInt32, (mx_control))
	Call(mxIsScalarInt32, (mx_nthread))

	Call(pot_name = mxu_string, (mx_pot_name,"pot_name"))

	cint nod_i = mxGetNumberOfDimensions(mx_x);
	Call(dim_i = mxGetDimensions, (mx_x))
	cuint nn = mxu_numel(mx_x);

	n_offset = mxGetM(mx_offsets) * mxGetN(mx_offsets);

	// trick: beta can be vector of length n_offset
	// or an array of size [n_offset nn]
	if (n_offset == mxGetM(mx_beta) && mxGetN(mx_beta) == 1)
		is_beta_array = 0;
	else if (n_offset == mxGetN(mx_beta) && mxGetM(mx_beta) == 1)
		is_beta_array = 0;
	else if (n_offset == mxGetM(mx_beta) && mxGetN(mx_beta) == nn)
		is_beta_array = 1;
	else
		Fail4("offsets(%d) nn=%d and beta size %d,%d mismatch",
			n_offset, nn, (int) mxGetM(mx_beta), (int) mxGetN(mx_beta))

#if 0 // not needed; it suffices to check that x and kappa have same numel
	cint nod_k = mxGetNumberOfDimensions(mx_kappa);
	if (nod_i != nod_k)
		Warn2("x %d kappa %d dimension mismatch", nod_i, nod_k)

	const mwSize *dim_k;
	Call(dim_k = mxGetDimensions, (mx_kappa))
	for (int ii=0; ii < nod_i; ii++)
		if (dim_k[ii] != dim_i[ii])
			Warn("x,kappa size mismatch")
#endif

	if (mxu_numel(mx_x) != mxu_numel(mx_kappa))
		Fail("x,kappa size mismatch")

#define the_args \
			nn, \
			mxGetPr_cfloat(mx_x), \
			mxGetPr_cfloat(mx_kappa), \
			mxGetPr_cint32(mx_offsets), \
			n_offset, \
			mxGetPr_cfloat(mx_beta), \
			is_beta_array, \
			pot_name, \
			mxGetPr_cfloat(mx_pot_param), \
			mxGetM(mx_pot_param) * mxGetN(mx_pot_param), \
			mxGetInt(mx_order), \
			mxGetInt(mx_control), \
			mxGetInt(mx_nthread), \
			Chat

	if (p_type == p_cgrad)
	{
		Call(plhs[0] = mxCreateNumericArray,
			(nod_i, dim_i, mxSINGLE_CLASS, mxREAL))
		Call(penalty_cgrad_offset, (
			(float *) mxGetData(plhs[0]), the_args))
	}

	else if (p_type == p_denom)
	{
		Call(plhs[0] = mxCreateNumericArray,
			(nod_i, dim_i, mxSINGLE_CLASS, mxREAL))
		Call(penalty_denom_offset, (
			(float *) mxGetData(plhs[0]), the_args))
	}

	else if (p_type == p_penal)
	{
		Call(plhs[0] = mxCreateNumericMatrix,
			(1, 1, mxDOUBLE_CLASS, mxREAL))
		Call(penalty_cost_offset, (
			mxGetPr(plhs[0]), /* scalar! */ the_args))
	}

	else
		Fail("bug")

#undef the_args

	Call(mxu_string_free, (pot_name))

	Ok
}


#if 0

// penalty_cgrad_mex()
// function('cgrad1', x, beta_array, kappa, offsets)
// single realization for now
static sof penalty_cgrad_mex(
int nlhs, mxArray *plhs[],
int nrhs, Const mxArray *prhs[],
ctruf order)
{
	Const mxArray *mx_x;
	Const mxArray *mx_beta;
	Const mxArray *mx_kappa;
	Const mxArray *mx_offsets;
	Const mxArray *mx_pot_name;
	Const mxArray *mx_pot_param;
	char *pot_name;
	int image_size = 1;
	int num_offsets = 1;
	int nod_i,ii,pot_name_size;
	const mwSize *dim_i;

/*	int nod_o;
	mwSize dim_o[Max_ndim]; */

	if (nlhs != 1 || nrhs != 7 )
	{
		penalty_mex_help();
		Call(mxu_arg, (nrhs, prhs))
		Fail(Usage)
	}
	mx_x = prhs[1];
	if (!mxIsRealSingle(mx_x)) Fail("x must be real single")
	mx_beta = prhs[2];
	if (!mxIsRealSingle(mx_beta)) Fail("beta must be real single")
	mx_kappa = prhs[3];
	if (!mxIsRealSingle(mx_kappa)) Fail("kappa must be real single")
	mx_offsets = prhs[4];
	if (!mxIsInt32(mx_offsets)) Fail("offsets must be int32")
	mx_pot_name = prhs[5];
	pot_name = mxu_string(mx_pot_name,"pot_name");
	pot_name_size = mxGetM(mx_pot_name) * mxGetN(mx_pot_name) + 1;
	mx_pot_param = prhs[6];
	if (!mxIsRealSingle(mx_pot_param)) Fail("pot_param must be single")

	nod_i = mxGetNumberOfDimensions(mx_x);
	Call(dim_i = mxGetDimensions, (mx_x))
	for(ii=0;ii<nod_i;ii++) image_size *= dim_i[ii];

	num_offsets = mxGetM(mx_offsets) * mxGetN(mx_offsets);

	Call(plhs[0] = mxCreateNumericArray,
		(nod_i, dim_i, mxSINGLE_CLASS, mxREAL))

	Call(penalty_cgrad, (
		image_size,
		num_offsets,
		(float *) mxGetData(plhs[0]),
		(cfloat *) mxGetData(mx_x),
		(cfloat *) mxGetData(mx_beta),
		(cfloat *) mxGetData(mx_kappa),
		(cint *) mxGetData(mx_offsets),
		(char *) pot_name,
		(cint) pot_name_size,
		(float *) mxGetData(mx_pot_param),
		(cint) mxGetM(mx_pot_param) * mxGetN(mx_pot_param),
		order) )

	Ok
}


// penalty_denom_mex()
// function('denom1', d, offsets, [lastdim])
static sof penalty_denom_mex(
int nlhs, mxArray *plhs[],
int nrhs, Const mxArray *prhs[],
ctruf order)
{
	Const mxArray *mx_x;
	Const mxArray *mx_beta;
	Const mxArray *mx_kappa;
	Const mxArray *mx_offsets;
	Const mxArray *mx_pot_name;
	Const mxArray *mx_pot_param;
	char *pot_name;
	int image_size = 1;
	int num_offsets = 1;
	int nod_i,ii,pot_name_size;
	const mwSize *dim_i;

/*	int nod_o;
	mwSize dim_o[Max_ndim]; */

	if (nlhs != 1 || nrhs != 7 )
	{
		penalty_mex_help();
		Call(mxu_arg, (nrhs, prhs))
		Fail(Usage)
	}
	mx_x = prhs[1];
	if (!mxIsRealSingle(mx_x)) Fail("x must be real single")
	mx_beta = prhs[2];
	if (!mxIsRealSingle(mx_beta)) Fail("beta must be real single")
	mx_kappa = prhs[3];
	if (!mxIsRealSingle(mx_kappa)) Fail("kappa must be real single")
	mx_offsets = prhs[4];
	if (!mxIsInt32(mx_offsets)) Fail("offsets must be int32")
	mx_pot_name = prhs[5];
	pot_name = mxu_string(mx_pot_name,"pot_name");
	pot_name_size = mxGetM(mx_pot_name) * mxGetN(mx_pot_name) + 1;
	mx_pot_param = prhs[6];
	if (!mxIsRealSingle(mx_pot_param)) Fail("pot_param must be single")

	nod_i = mxGetNumberOfDimensions(mx_x);
	Call(dim_i = mxGetDimensions, (mx_x))
	for(ii=0;ii<nod_i;ii++) image_size *= dim_i[ii];

	num_offsets = mxGetM(mx_offsets) * mxGetN(mx_offsets);

	Call(plhs[0] = mxCreateNumericArray,
		(nod_i, dim_i, mxSINGLE_CLASS, mxREAL))

	Call(penalty_denom, (
		image_size,
		num_offsets,
		(float *) mxGetData(plhs[0]),
		(cfloat *) mxGetData(mx_x),
		(cfloat *) mxGetData(mx_beta),
		(cfloat *) mxGetData(mx_kappa),
		(cint *) mxGetData(mx_offsets),
		(char *) pot_name,
		(cint) pot_name_size,
		(float *) mxGetData(mx_pot_param),
		(cint) mxGetM(mx_pot_param) * mxGetN(mx_pot_param),
		order) )

	Ok
}


// penalty_cgrad_denom_mex()
// function('cgrad_denom1', d, offsets, [lastdim])
static sof penalty_cgrad_denom_mex(
int nlhs, mxArray *plhs[],
int nrhs, Const mxArray *prhs[],
ctruf order)
{
	Const mxArray *mx_x;
	Const mxArray *mx_beta;
	Const mxArray *mx_kappa;
	Const mxArray *mx_offsets;
	Const mxArray *mx_pot_name;
	Const mxArray *mx_pot_param;
	Const mxArray *mx_numthreads;
	char *pot_name;
	int image_size = 1;
	int num_offsets = 1;
	int nod_i,ii,pot_name_size;
	double *numthreads;
	const mwSize *dim_i;

	int nod_o;
	mwSize dim_o[Max_ndim];

	if (nlhs != 1 || nrhs != 8 )
	{
		penalty_mex_help();
		Call(mxu_arg, (nrhs, prhs))
		Fail(Usage)
	}
	mx_x = prhs[1];
	if (!mxIsRealSingle(mx_x)) Fail("x must be real single")
	mx_beta = prhs[2];
	if (!mxIsRealSingle(mx_beta)) Fail("beta must be real single")
	mx_kappa = prhs[3];
	if (!mxIsRealSingle(mx_kappa)) Fail("kappa must be real single")
	mx_offsets = prhs[4];
	if (!mxIsInt32(mx_offsets)) Fail("offsets must be int32")
	mx_pot_name = prhs[5];
	pot_name = mxu_string(mx_pot_name,"pot_name");
	pot_name_size = mxGetM(mx_pot_name) * mxGetN(mx_pot_name) + 1;
	mx_pot_param = prhs[6];
	if (!mxIsRealSingle(mx_pot_param)) Fail("pot_param must be single")
	mx_numthreads = prhs[7];
	if (!mxIsDouble(mx_numthreads))
		Fail("numthreads must be double")

	nod_i = mxGetNumberOfDimensions(mx_x);
	nod_o = nod_i;

	Call(dim_i = mxGetDimensions, (mx_x))
	for(ii=0;ii<nod_i;ii++)
	{
		image_size *= dim_i[ii];
		if (ii==0)
			dim_o[ii] = 2*dim_i[ii];
		else
			dim_o[ii] = dim_i[ii];
	}

	num_offsets = mxGetM(mx_offsets) * mxGetN(mx_offsets);
	numthreads = (double *) mxGetData(mx_numthreads);

	Call(plhs[0] = mxCreateNumericArray,
		(nod_o, dim_o, mxSINGLE_CLASS, mxREAL))

	/*
	Call(penalty_cgrad_denom, (
		image_size,
		num_offsets,
		(float *) mxGetData(plhs[0]),
		(cfloat *) mxGetData(mx_x),
		(cfloat *) mxGetData(mx_beta),
		(cfloat *) mxGetData(mx_kappa),
		(cint *) mxGetData(mx_offsets),
		(char *) pot_name,
		(cint) pot_name_size,
		(float *) mxGetData(mx_pot_param),
		(cint) mxGetM(mx_pot_param) * mxGetN(mx_pot_param),
		order) )
	*/

	Call(penalty_cgrad_denom_par, (
		(int) *numthreads,
		image_size,
		num_offsets,
		(float *) mxGetData(plhs[0]),
		(float *) mxGetData(mx_x),
		(float *) mxGetData(mx_beta),
		(float *) mxGetData(mx_kappa),
		(int *) mxGetData(mx_offsets),
		(char *) pot_name,
		(int) pot_name_size,
		(float *) mxGetData(mx_pot_param),
		(int) mxGetM(mx_pot_param) * mxGetN(mx_pot_param),
		order) )

	Ok
}


// penalty_cgrad_denom_big_mex()
// function('cgrad_denom1', d, offsets, [lastdim])
static sof penalty_cgrad_denom_big_mex(
int nlhs, mxArray *plhs[],
int nrhs, Const mxArray *prhs[],
ctruf order)
{
	Const mxArray *mx_x;
	Const mxArray *mx_beta;
	Const mxArray *mx_kappa;
	Const mxArray *mx_offsets;
	Const mxArray *mx_pot_name;
	Const mxArray *mx_pot_param;
	Const mxArray *mx_numthreads;
	Const mxArray *mx_cgrad_denom;
	void *dummy;
	char *pot_name;
	int image_size = 1;
	int num_offsets = 1;
	int nod_i, pot_name_size;
	double *numthreads;
	const mwSize *dim_i;

	int nod_o;
	mwSize dim_o[Max_ndim];

	dummy = (void *) plhs;
	if (nlhs != 0 || nrhs != 9 )
	{
		penalty_mex_help();
		Call(mxu_arg, (nrhs, prhs))
		Fail(Usage)
	}
	mx_x = prhs[1];
	if (!mxIsRealSingle(mx_x)) Fail("x must be real single")
	mx_beta = prhs[2];
	if (!mxIsRealSingle(mx_beta)) Fail("beta must be real single")
	mx_kappa = prhs[3];
	if (!mxIsRealSingle(mx_kappa)) Fail("kappa must be real single")
	mx_offsets = prhs[4];
	if (!mxIsInt32(mx_offsets)) Fail("offsets must be int32")
	mx_pot_name = prhs[5];
	pot_name = mxu_string(mx_pot_name,"pot_name");
	pot_name_size = mxGetM(mx_pot_name) * mxGetN(mx_pot_name) + 1;
	mx_pot_param = prhs[6];
	if (!mxIsRealSingle(mx_pot_param)) Fail("pot_param must be single")
	mx_numthreads = prhs[7];
	if (!mxIsDouble(mx_numthreads))
		Fail("numthreads must be double")
	mx_cgrad_denom = prhs[8];
	if (!mxIsRealSingle(mx_cgrad_denom)) Fail("cgrad_denom must be real single")


	nod_i = mxGetNumberOfDimensions(mx_x);
	nod_o = nod_i;

	Call(dim_i = mxGetDimensions, (mx_x))
	for (int ii=0;ii<nod_i;ii++)
	{
		image_size *= dim_i[ii];
		if (ii==0)
			dim_o[ii] = 2*dim_i[ii];
		else
			dim_o[ii] = dim_i[ii];
	}

	num_offsets = mxGetM(mx_offsets) * mxGetN(mx_offsets);
	numthreads = (double *) mxGetData(mx_numthreads);
	/*
	Call(plhs[0] = mxCreateNumericArray,
		(nod_o, dim_o, mxSINGLE_CLASS, mxREAL))
	*/
	/*
	Call(penalty_cgrad_denom, (
		image_size,
		num_offsets,
		(float *) mxGetData(mx_cgrad_denom),
		(cfloat *) mxGetData(mx_x),
		(cfloat *) mxGetData(mx_beta),
		(cfloat *) mxGetData(mx_kappa),
		(cint *) mxGetData(mx_offsets),
		(char *) pot_name,
		(cint) pot_name_size,
		(float *) mxGetData(mx_pot_param),
		(cint) mxGetM(mx_pot_param) * mxGetN(mx_pot_param),
		order) )
	*/

	Call(penalty_cgrad_denom_par, (
		(int) *numthreads,
		image_size,
		num_offsets,
		(float *) mxGetData(mx_cgrad_denom),
		(float *) mxGetData(mx_x),
		(float *) mxGetData(mx_beta),
		(float *) mxGetData(mx_kappa),
		(int *) mxGetData(mx_offsets),
		(char *) pot_name,
		(int) pot_name_size,
		(float *) mxGetData(mx_pot_param),
		(int) mxGetM(mx_pot_param) * mxGetN(mx_pot_param),
		order) )

	Ok
}

#endif // if 0


// reg_cgrad_zxy_mex()
// g = function('cgrad,zxy',
//	x, kappa, offsets, betas, pot_name, pot_param,
//	mask2, order, nthread)
static sof reg_cgrad_zxy_mex(
mxArray *plhs[],
Cmx mx_x,
Cmx mx_kappa,
Cmx mx_offsets,
Cmx mx_betas,
Cmx mx_pot_name,
Cmx mx_pot_param,
Cmx mx_mask2,
Cmx mx_order,
Cmx mx_nthread)
{
	char *pot_name;
	int nod_i = mxGetNumberOfDimensions(mx_x);
	const mwSize *dim_i;
	cint chat = 0;
	int nx, ny, nz;
	reg_mask2 rm_, *rm = &rm_;

	if (!mxIsRealSingle(mx_x)) Fail("x must be real single")
	if (!mxIsRealSingle(mx_kappa)) Fail("kappa must be real single")
	if (!mxIsRealSingle(mx_betas)) Fail("beta must be real single")
	if (!mxIsInt32(mx_offsets)) Fail("offsets must be int32")
	pot_name = mxu_string(mx_pot_name, "pot_name");
	if (!mxIsRealSingle(mx_pot_param)) Fail("pot_param must be single")
	if (!mxIsUint8(mx_mask2)) Fail("mask2 must be uint8")
	if (!mxIsScalarInt32(mx_order)) Fail("order must be int32 scalar")
	if (!mxIsScalarInt32(mx_nthread)) Fail("nthread must be int32 scalar")
	cint nthread = mxGetInt(mx_nthread);

	Call(dim_i = mxGetDimensions, (mx_x))

	if (nod_i != 3) Fail("only 3D done now")
	if (3 != mxGetNumberOfDimensions(mx_kappa))
		Fail("kappa must be 3D")

	nz = dim_i[0];
	nx = dim_i[1];
	ny = dim_i[2];

	if (nz != (int) mxGetDimensions(mx_kappa)[0]
	 || nx != (int) mxGetDimensions(mx_kappa)[1]
	 || ny != (int) mxGetDimensions(mx_kappa)[2])
		Fail("x and kappa dimension mismatch")
	if (nx != (int) mxGetDimensions(mx_mask2)[0]
	 || ny != (int) mxGetDimensions(mx_mask2)[1])
		Fail("x and mask2 dimension mismatch")

	cbyte *mask2 = (cbyte *) mxGetData(mx_mask2); // binary
	byte *mask2_id; // 0, 1 ... nthread
	Mem0(mask2_id, nx*ny)
	Call(jf_thread1_mask_init, (mask2_id, mask2, nx*ny, nthread, Chat))

	Call(reg_mask2_set, (rm, mask2_id, nx, ny))

	Call(plhs[0] = mxCreateNumericArray,
		(nod_i, dim_i, mxSINGLE_CLASS, mxREAL))

#define reg_cgrad_zxy_mex_args \
	(float *) mxGetData(plhs[0]), /* cgrad */ \
	dim_i[1] /* nx */, dim_i[2] /* ny */, dim_i[0] /* nz */, \
	(cfloat *) mxGetData(mx_x), \
	(cfloat *) mxGetData(mx_kappa), \
	(cint *) mxGetData(mx_offsets), \
	mxGetM(mx_offsets) * mxGetN(mx_offsets), \
	(cfloat *) mxGetData(mx_betas), \
	pot_name, \
	(cfloat *) mxGetData(mx_pot_param), \
	mxGetM(mx_pot_param) * mxGetN(mx_pot_param), \
	rm, \
	mxGetInt(mx_order), \
	reg_end_replicate

	jf_timeval jt_, *jt = &jt_;
	Call(jf_time_tic, (jt))
	Call(reg_cgrad_zxy, (reg_cgrad_zxy_mex_args, nthread, chat))
	Call(jf_time_toc, (jt, "reg_cgrad_zxy"))
//	Note2("nthread=%d time=%g", nthread, jf_clock)

	Free0(mask2_id)
	Call(mxu_string_free, (pot_name))
	Ok
}


// reg_cgrad_denom_zxy_mex()
// [g w] = function('cgrad,denom,zxy',
//	x, kappa, offsets, betas, pot_name, pot_param,
//	mask2, order, nthread)
static sof reg_cgrad_denom_zxy_mex(
mxArray *plhs[],
Cmx mx_x,
Cmx mx_kappa,
Cmx mx_offsets,
Cmx mx_betas,
Cmx mx_pot_name,
Cmx mx_pot_param,
Cmx mx_mask2,
Cmx mx_order,
Cmx mx_nthread)
{
	char *pot_name;
	int nod_i = mxGetNumberOfDimensions(mx_x);
	const mwSize *dim_i;
	cint chat = 0;
	int nx, ny, nz;
	reg_mask2 rm_, *rm = &rm_;

	if (!mxIsRealSingle(mx_x)) Fail("x must be real single")
	if (!mxIsRealSingle(mx_kappa)) Fail("kappa must be real single")
	if (!mxIsRealSingle(mx_betas)) Fail("beta must be real single")
	if (!mxIsInt32(mx_offsets)) Fail("offsets must be int32")
	pot_name = mxu_string(mx_pot_name, "pot_name");
	if (!mxIsRealSingle(mx_pot_param)) Fail("pot_param must be single")
	if (!mxIsUint8(mx_mask2)) Fail("mask2 must be uint8")
	if (!mxIsScalarInt32(mx_order)) Fail("order must be int32 scalar")
	if (!mxIsScalarInt32(mx_nthread)) Fail("nthread must be int32 scalar")
	cint nthread = mxGetInt(mx_nthread);

	Call(dim_i = mxGetDimensions, (mx_x))

	if (nod_i != 3) Fail("only 3D done now")
	if (3 != mxGetNumberOfDimensions(mx_kappa))
		Fail("kappa must be 3D")

	nz = dim_i[0];
	nx = dim_i[1];
	ny = dim_i[2];

	if (nz != (int) mxGetDimensions(mx_kappa)[0]
	 || nx != (int) mxGetDimensions(mx_kappa)[1]
	 || ny != (int) mxGetDimensions(mx_kappa)[2])
		Fail("x and kappa dimension mismatch")
	if (nx != (int) mxGetDimensions(mx_mask2)[0]
	 || ny != (int) mxGetDimensions(mx_mask2)[1])
		Fail("x and mask2 dimension mismatch")

	cbyte *mask2 = (cbyte *) mxGetData(mx_mask2); // binary
	byte *mask2_id; // 0, 1 ... nthread
	Mem0(mask2_id, nx*ny)
	Call(jf_thread1_mask_init, (mask2_id, mask2, nx*ny, nthread, Chat))

	Call(reg_mask2_set, (rm, mask2_id, nx, ny))

	Call(plhs[0] = mxCreateNumericArray,
		(nod_i, dim_i, mxSINGLE_CLASS, mxREAL))
	Call(plhs[1] = mxCreateNumericArray,
		(nod_i, dim_i, mxSINGLE_CLASS, mxREAL))

#define reg_cgrad_denom_zxy_mex_args \
	(float *) mxGetData(plhs[0]), /* cgrad */ \
	(float *) mxGetData(plhs[1]), /* denom */ \
	dim_i[1] /* nx */, dim_i[2] /* ny */, dim_i[0] /* nz */, \
	(cfloat *) mxGetData(mx_x), \
	(cfloat *) mxGetData(mx_kappa), \
	(cint *) mxGetData(mx_offsets), \
	mxGetM(mx_offsets) * mxGetN(mx_offsets), \
	(cfloat *) mxGetData(mx_betas), \
	pot_name, \
	(cfloat *) mxGetData(mx_pot_param), \
	mxGetM(mx_pot_param) * mxGetN(mx_pot_param), \
	rm, \
	mxGetInt(mx_order), \
	reg_end_replicate

	jf_timeval jt_, *jt = &jt_;
	Call(jf_time_tic, (jt))
	Call(reg_cgrad_denom_zxy, (reg_cgrad_denom_zxy_mex_args, nthread, chat))
	Call(jf_time_toc, (jt, "reg_cgrad_denom_zxy"))
//	Note2("nthread=%d time=%g", nthread, jf_clock)

	Free0(mask2_id)
	Call(mxu_string_free, (pot_name))
	Ok
}


// power_parse()
static sof power_parse(int *ppow, cchar *str, cchar *arg)
{
	if (Streq(str, "1"))
		*ppow = 1;
	else if (Streq(str, "2"))
		*ppow = 2;
	else if (Streq(str, "A"))
		*ppow = -1; /* trick: |+/- 1| = (+/- 1)^2 */
	else
		Fail1("unknown (power?) '%s'", arg)
	Ok
}


// intermediate gateway routine
static sof penalty_mex(
int nlhs, mxArray *plhs[],
int nrhs, Const mxArray *prhs[])
{
	char *arg;

	if (nrhs < 1 || !mxIsChar(prhs[0]))
	{
		penalty_mex_help();
		Call(mxu_arg, (nrhs, prhs))
		Fail(Usage)
	}

	if (nlhs == 0 && nrhs == 1 && mxu_streq(prhs[0], "arg 1", "check"))
                Ok // check

	// 'help'
	if (nrhs == 1 && mxIsChar(prhs[0]))
	{
		penalty_mex_help();
		Ok
	}

	Call(arg = mxu_string, (prhs[0], "1st argument"))

	// diff forward: C * x
	if (Streqn(arg, "diff1,forw", 10))
	{
		int power;
		Call(power_parse, (&power, arg+10, arg))
		Call(penalty_diff_forw_mex, (nlhs, plhs, nrhs, prhs, 1, power))
	}

	else if (Streqn(arg, "diff2,forw", 10))
	{
		int power;
		Call(power_parse, (&power, arg+10, arg))
		Call(penalty_diff_forw_mex, (nlhs, plhs, nrhs, prhs, 2, power))
	}

	// diff back: C' * d
	else if (Streqn(arg, "diff1,back", 10))
	{
		int power;
		Call(power_parse, (&power, arg+10, arg))
		Call(penalty_diff_back_mex, (nlhs, plhs, nrhs, prhs, 1, power))
	}

	else if (Streqn(arg, "diff2,back", 10))
	{
		int power;
		Call(power_parse, (&power, arg+10, arg))
		Call(penalty_diff_back_mex, (nlhs, plhs, nrhs, prhs, 2, power))
	}

	// new low memory options
	// (power absorbed into beta vector)

	else if (Streq(arg, "cgrad,offset"))
	{
		if (nlhs != 1 || nrhs != 10) FailUse
		Call(penalty_cdp_offset_mex, (plhs,
			prhs[1], prhs[2], prhs[3], prhs[4],
			prhs[5], prhs[6], prhs[7], prhs[8], prhs[9], p_cgrad))
	}

	else if (Streq(arg, "denom,offset"))
	{
		if (nlhs != 1 || nrhs != 10) FailUse
		Call(penalty_cdp_offset_mex, (plhs,
			prhs[1], prhs[2], prhs[3], prhs[4],
			prhs[5], prhs[6], prhs[7], prhs[8], prhs[9], p_denom))
	}

	else if (Streq(arg, "penal,offset"))
	{
		if (nlhs != 1 || nrhs != 10) FailUse
		Call(penalty_cdp_offset_mex, (plhs,
			prhs[1], prhs[2], prhs[3], prhs[4],
			prhs[5], prhs[6], prhs[7], prhs[8], prhs[9], p_penal))
	}

#if 0
xxx
	else if (Streqn(arg, "cgrad1", 6)) /* power absorbed into beta */
		Call(penalty_cgrad_mex, (nlhs, plhs, nrhs, prhs, 1))
	else if (Streqn(arg, "cgrad2", 6))
		Call(penalty_cgrad_mex, (nlhs, plhs, nrhs, prhs, 2))
	else if (Streqn(arg, "denom1", 6))
		Call(penalty_denom_mex, (nlhs, plhs, nrhs, prhs, 1))
	else if (Streqn(arg, "denom2", 6))
		Call(penalty_denom_mex, (nlhs, plhs, nrhs, prhs, 2))
	else if (Streqn(arg, "cgrad_denom1", 6))
	{
		if (nlhs == 1)
			Call(penalty_cgrad_denom_mex, (nlhs, plhs, nrhs, prhs, 1))
		else
			Call(penalty_cgrad_denom_big_mex, (nlhs, plhs, nrhs, prhs, 1))
	}
#endif

	// wk factors for \sum_k w_k \pot([Cx]_k)
	else if (Streqn(arg, "wk,", 3))
		Call(penalty_wk_mex, (nlhs, plhs, nrhs, prhs))

	else if (Streq(arg, "cgrad,zxy"))
	{
		if (nlhs != 1 || nrhs != 10) FailUse
		Call(reg_cgrad_zxy_mex, (plhs,
			prhs[1], prhs[2], prhs[3], prhs[4],
			prhs[5], prhs[6], prhs[7], prhs[8], prhs[9]))
	}

	else if (Streq(arg, "cgrad,denom,zxy"))
	{
		if (nlhs != 2 || nrhs != 10) FailUse
		Call(reg_cgrad_denom_zxy_mex, (plhs,
			prhs[1], prhs[2], prhs[3], prhs[4],
			prhs[5], prhs[6], prhs[7], prhs[8], prhs[9]))
	}

	else
		Fail1("command '%s' unknown", arg)

	Call(mxu_string_free, (arg))
	Ok
}


// gateway routine
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, Const mxArray *prhs[])
{
	if (!nlhs && !nrhs)
	{
		penalty_mex_help();
		return;
	}
	if (!penalty_mex(nlhs, plhs, nrhs, prhs))
		mexErrMsgTxt("penalty_mex");
}
