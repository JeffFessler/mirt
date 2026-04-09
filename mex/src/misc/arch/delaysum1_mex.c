// delaysum1_mex.c
// Matlab mex gateway for delay-sum algorithm
//
// Copyright 2009-11-18, Jeff Fessler, University of Michigan

#include "jf,mex,def.h"

#if !defined(Need_delaysum1_mex_gateway)
#include "delaysum1,def.h"
#endif

#define Usage "usage error. see above"

static void delaysum1_mex_help(void)
{
	printf("\n\
\n\
y[n;m] = sum_{k=0}^{Nx-1} x[k] h[n - d[k;m]], n=0,...,Ny-1, m=0,...,Nm-1\n\
for 0 <= n - d[k;m] <= Nh - 1\n\
\n\
y = function('delaysum1,forw', h, d, nthread, x, Ny)\n\
	y: (single) [Ny Nm]\n\
\n\
	h: (single) [Nh 1] function to be delayed\n\
	d: (int32) [Nx Nm] delays (positive or negative)\n\
	x: (single) [Nx 1] coefficients\n\
	nthread: (int32) # of threads\n\
	Ny: (int32) output size\n\
\n\
x = function('delaysum1,back', h, d, nthread, y)\n\
	x: (single) [Nx 1]\n\
	y: (single) [Ny Nm] or [Ny*Nm 1]\n\
\n");
}


//
// delaysum1_back()
//
static sof delaysum1_back(
cfloat *yy, // [Ny Nm]
cint Ny,
cint Nm,
cfloat *hh, // [Nh]
cint Nh,
cfloat *dd, // [Nx Nm]
cint Nx,
float *p_xk, // [Nx]
cint nthread)
{
	if (nthread != 1) Warn("only 1 thread done")

#if 1
	for (int id=0; id < Nx; ++id) {
		double sum = 0;
		for (int im=0; im < Nm; ++im) {
			int delay = *(dd + Nx * im);
			cfloat *py = yy + Max(delay,0) + Ny * im;
			cfloat *ph = hh + Max(-delay,0);
			cint nh_new = Nh - Max(-delay,0);
			cint ny_new = Ny - Max(delay,0);
			cint npoint = Min(nh_new, ny_new);
			double sum1;
			VectInprod(float, sum1, py, ph, npoint)
			sum += sum1;
		}
		*p_xk++ = sum;
		dd++;
	}

#else
	for (int im=0; im < Nm; ++im, yy += Ny) {
		float *xk = p_xk;
		for (int id=0; id < Nx; ++id) {
			int delay = *dd++;
			cfloat *py = yy + Max(delay,0);
			cfloat *ph = hh + Max(-delay,0);
			cint nh_new = Nh - Max(-delay,0);
			cint ny_new = Ny - Max(delay,0);
			cint npoint = Min(nh_new, ny_new);
			float sum;
			VectInprod(float, sum, py, ph, npoint)
			*xk++ += sum;
		}
	}
#endif

	Ok
}


//
// delaysum1_mex_back()
//
static sof delaysum1_mex_back(
mxArray *plhs[], // [Nx]
Cmx mx_h, // [Nh]
Cmx mx_delay, // [Nx Nm]
Cmx mx_nthread,
Cmx mx_y) // [Ny Nm] or [Ny*Nm 1]
{
	Call(mxIsRealSingle, (mx_h))
	cint Nh = mxGetM(mx_h) * mxGetN(mx_h);

	Call(mxIsRealSingle, (mx_delay))
	cint Nx = mxGetM(mx_delay);
	cint Nm = mxGetN(mx_delay);

	Call(mxIsScalarInt32, (mx_nthread))
	cint nthread = mxGetInt(mx_nthread);

	Call(mxIsRealSingle, (mx_y))
	cint NyNm = mxGetM(mx_y) * mxGetN(mx_y);
	cint Ny = NyNm / Nm;
	if (Ny * Nm != NyNm)
		Fail3("d %d %d and y %d size mismatch", Nx, Nm, NyNm)

	// output
	Call(plhs[0] = mxCreateNumericMatrix, (Nx, 1, mxSINGLE_CLASS, mxREAL))

	Call(delaysum1_back, ( mxGetPr_cfloat(mx_y),
			Ny, Nm,
			mxGetPr_cfloat(mx_h), Nh,
			mxGetPr_cfloat(mx_delay), Nx,
			(float *) mxGetData(plhs[0]), // y
			nthread))
	Ok
}



//
// delaysum1_forw()
//
static sof delaysum1_forw(
float *yy, // [Ny Nm]
cint Ny,
cint Nm,
cfloat *hh, // [Nh]
cint Nh,
cfloat *dd, // [Nx Nm]
cint Nx,
cfloat *p_xk, // [Nx]
cint nthread)
{
	if (nthread != 1) Warn("only 1 thread done")

	for (int im=0; im < Nm; ++im, yy += Ny) {
		cfloat *xk = p_xk;
		for (int id=0; id < Nx; ++id) {
			int delay = *dd++;
			cfloat coef = *xk++;
			float *py = yy + Max(delay,0);
			cfloat *ph = hh + Max(-delay,0);
			cint nh_new = Nh - Max(-delay,0);
			cint ny_new = Ny - Max(delay,0);
			cint npoint = Min(nh_new, ny_new);
			VectAddScale(float, py, coef, ph, npoint)
		}
	}
	Ok
}


//
// delaysum1_mex_forw()
//
static sof delaysum1_mex_forw(
mxArray *plhs[], // [Ny Nm]
Cmx mx_h, // [Nh]
Cmx mx_delay, // [Nx Nm]
Cmx mx_nthread,
Cmx mx_x, // [Nx]
Cmx mx_Ny)
{
	Call(mxIsRealSingle, (mx_h))
	cint Nh = mxGetM(mx_h) * mxGetN(mx_h);

	Call(mxIsRealSingle, (mx_delay))
//	cint Nx = mxGetM(mx_delay);
	cint Nm = mxGetN(mx_delay);

	Call(mxIsScalarInt32, (mx_nthread))
	cint nthread = mxGetInt(mx_nthread);

	Call(mxIsRealSingle, (mx_x))
	cint Nx = mxGetM(mx_x) * mxGetN(mx_x);
	if (Nx != (int) mxGetM(mx_delay))
		Fail3("d %d %d and x %d size mismatch",
			(int) mxGetM(mx_delay), Nm, Nx)

	Call(mxIsScalarInt32, (mx_Ny))
	cint Ny = mxGetInt(mx_Ny);

	// output
	Call(plhs[0] = mxCreateNumericMatrix,
		(Ny*Nm, 1, mxSINGLE_CLASS, mxREAL))

	Call(delaysum1_forw, ( (float *) mxGetData(plhs[0]), Ny, Nm,
			mxGetPr_cfloat(mx_h), Nh,
			mxGetPr_cfloat(mx_delay), Nx,
			mxGetPr_cfloat(mx_x),
			nthread))
	Ok
}


// intermediate gateway routine
#if defined(Need_delaysum1_mex_gateway)
static
#endif
sof delaysum1_mex(
int nlhs, mxArray *plhs[],
int nrhs, Const mxArray *prhs[])
{
	char *arg;

	if (nrhs <= 1 || !mxIsChar(prhs[0])) {
		delaysum1_mex_help();
		Call(mxu_arg, (nrhs, prhs))
		Fail(Usage)
	}

	// help
	if (nrhs == 1 && mxIsChar(prhs[0])) {
		delaysum1_mex_help();
		Ok
	}

	Call(arg = mxu_string, (prhs[0], "1st argument"))

	// forw
	if (Streq(arg, "delaysum1,forw")) {
		if (nrhs != 6 || nlhs != 1)
			Fail(Usage)
		Call(delaysum1_mex_forw, (plhs,
			prhs[1], prhs[2], prhs[3], prhs[4], prhs[5]))
	}

	// back
	else if (Streq(arg, "delaysum1,back")) {
		if (nrhs != 5 || nlhs != 1)
			Fail(Usage)
		Call(delaysum1_mex_back, (plhs,
			prhs[1], prhs[2], prhs[3], prhs[4]))
	}

	else
		Fail1("command '%s' unknown", arg)

	Call(mxu_string_free, (arg))
	Ok
}


#if defined(Need_delaysum1_mex_gateway)
// gateway routine
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, Const mxArray *prhs[])
{
	if (!nlhs && !nrhs) {
		delaysum1_mex_help();
		return;
	}
	if (!delaysum1_mex(nlhs, plhs, nrhs, prhs))
		mexErrMsgTxt("delaysum1_mex");
}
#endif
