// delaysum1_mex.c
// Matlab mex gateway for delay-sum algorithm
//
// Copyright 2009-11-18, Jeff Fessler, University of Michigan
//
// Added threading, 2012-01-26, Steve Schmitt

#include "jf,mex,def.h"
#include "jf,thread1.h"
// #include <unistd.h>

#if !defined(Need_delaysum1_mex_gateway)
#include "delaysum1,def.h"
#endif

#define Usage "usage error. see above"

static void delaysum1_mex_help(void)
{
	// int nc = jf_thread1_ncore(-1);
	printf("\n\
\n\
delaysum1_mex usage:\n\
\n\
y[n;m] = sum_{k=0}^{Nx-1} x[k] h[n - d[k;m]], n=0,...,Ny-1, m=0,...,Nm-1\n\
for 0 <= n - d[k;m] <= Nh - 1\n\
\n\
y = function('delaysum1,forw[,thr]', h, d, nthread, x, Ny, chat)\n\
	y: (single) [Ny Nm]\n\
\n\
	h: (single) [Nh 1] function to be delayed\n\
	d: (single) [Nx Nm] delays (positive or negative)\n\
	nthread: (int32) # of threads\n\
	x: (single) [Nx 1] coefficients\n\
	Ny: (int32) output size\n\
	chat: (int32) be verbose with threading\n\
\n\
x = function('delaysum1,back[,thr]', h, d, nthread, y, chat)\n\
	x: (single) [Nx 1]\n\
\n\
	y: (single) [Ny Nm] or [Ny*Nm 1]\n\
\n\
	This is the threaded 2012 version.\n\
\n");
}


// delaysum1_back()
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
	for (int id=0; id < Nx; ++id)
	{
		double sum = 0;
		for (int im=0; im < Nm; ++im)
		{
			int delay = (int) *(dd + Nx * im);
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
	for (int im=0; im < Nm; ++im, yy += Ny)
	{
		float *xk = p_xk;
		for (int id=0; id < Nx; ++id)
		{
			int delay = (int) *dd++;
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


// delaysum1_mex_back()
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


typedef struct
{
	float *yy;
	int Ny;
	int Nm;
	float *hh;
	int Nh;
	float *dd;
	int Nx;
	float *p_xk;
	int chat;
	int* xstart;
	int* xcount;
} delaysum_back_worker_args;


// delaysum1_back_worker()
static sof delaysum1_back_worker(void* p, cint tid, cint nthread)
{
	delaysum_back_worker_args* args = (delaysum_back_worker_args*) p;
	
	if (args->chat)
		printf("delaysum1_back_worker(%x,%d,%d)\n",p,tid,nthread);
	// try to fast-forward pointers as they should be after m_loops outer loops
	// already ran.
	float *xk = args->p_xk + args->xstart[tid];
	float *dd = args->dd + args->xstart[tid];

	for (int id=0; id < args->xcount[tid]; ++id)
	{
		double sum = 0;
		for (int im=0; im < args->Nm; ++im)
		{
			int delay = (int) *(dd + args->Nx * im);
			cfloat *py = args->yy + Max(delay,0) + args->Ny * im;
			cfloat *ph = args->hh + Max(-delay,0);
			cint nh_new = args->Nh - Max(-delay,0);
			cint ny_new = args->Ny - Max(delay,0);
			cint npoint = Min(nh_new, ny_new);
			double sum1;
			VectInprod(float, sum1, py, ph, npoint)
			sum += sum1;
		}
		*xk++ = sum;
		dd++;
	}

	Ok
}


// delaysum1_back_thr()
static sof delaysum1_back_thr(
cfloat *yy, // [Ny Nm]
cint Ny,
cint Nm,
cfloat *hh, // [Nh]
cint Nh,
cfloat *dd, // [Nx Nm]
cint Nx,
float *p_xk, // [Nx]
cint nthread,
cint chat)
{
	delaysum_back_worker_args args;
	args.yy = (float*) yy;
	args.Ny = Ny;
	args.Nm = Nm;
	args.hh = (float*) hh;
	args.Nh = Nh;
	args.dd = (float*) dd;
	args.Nx = Nx;
	args.p_xk = p_xk;
	args.chat = chat;

	Mem0(args.xstart, nthread)
	Mem0(args.xcount, nthread)

	cint x_per_thread = Nx / nthread;
	cint x_leftover = Nx % nthread;
	int x_claimed = 0;

	for (int i = 0; i < nthread; i++)
	{
		args.xstart[i] = x_claimed;
		int this_claim = x_per_thread + (i < x_leftover ? 1 : 0);
		args.xcount[i] = this_claim;
		x_claimed += this_claim;
		if (chat)
			printf("delaysum1_back_thr: claiming %d loops for thread %d.\n",
				this_claim, i);
	}

	if (x_claimed != Nx)
	{
		Fail2("delaysum1_back_thr claimed %d but needed to claim %d.",
			x_claimed, Nx);
	}
	
	jf_thread1_top(delaysum1_back_worker, NULL, &args, nthread, 0);
	
	Free0(args.xstart)
	Free0(args.xcount)
	Ok
}


// delaysum1_mex_back_thr()
static sof delaysum1_mex_back_thr(
mxArray *plhs[], // [Nx]
Cmx mx_h, // [Nh]
Cmx mx_delay, // [Nx Nm]
Cmx mx_nthread,
Cmx mx_y, // [Ny Nm] or [Ny*Nm 1]
Cmx mx_chat)
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

	Call(mxIsScalarInt32, (mx_chat))
	cint chat = mxGetInt(mx_chat);

	// output
	Call(plhs[0] = mxCreateNumericMatrix, (Nx, 1, mxSINGLE_CLASS, mxREAL))

	Call(delaysum1_back_thr, ( mxGetPr_cfloat(mx_y),
			Ny, Nm,
			mxGetPr_cfloat(mx_h), Nh,
			mxGetPr_cfloat(mx_delay), Nx,
			(float *) mxGetData(plhs[0]), // y
			nthread, chat))
	Ok
}


// delaysum1_forw()
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

	for (int im=0; im < Nm; ++im, yy += Ny)
	{
		cfloat *xk = p_xk;
		for (int id=0; id < Nx; ++id)
		{
			int delay = (int) *dd++;
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
// Threaded delaysum1_forw:
//

typedef struct
{
	float *yy;
	int Ny;
	int Nm;
	float *hh;
	int Nh;
	float *dd;
	int Nx;
	float *p_xk;
	int *mstart;
	int *mcount;
	int chat;
} delaysum1_forw_worker_args;

sof delaysum1_forw_worker(void* p, cint id, cint nthread);


// delaysum1_forw_thr()
static sof delaysum1_forw_thr(
float *yy, // [Ny Nm]
cint Ny,
cint Nm,
cfloat *hh, // [Nh]
cint Nh,
cfloat *dd, // [Nx Nm]
cint Nx,
cfloat *p_xk, // [Nx]
cint nthread,
cint chat)
{
	delaysum1_forw_worker_args args;
	args.yy = yy;
	args.Ny = Ny;
	args.Nm = Nm;
	args.hh = (float*) hh;
	args.Nh = Nh;
	args.dd = (float*) dd;
	args.Nx = Nx;
	args.p_xk = (float*) p_xk;
	args.chat = chat;

	Mem0(args.mstart, nthread)
	Mem0(args.mcount, nthread)

	cint m_per_thread = Nm / nthread;
	cint m_leftover = Nm % nthread;
	int m_claimed = 0;

	for (int i = 0; i < nthread; i++)
	{
		args.mstart[i] = m_claimed;
		int this_claim = m_per_thread + (i < m_leftover ? 1 : 0);
		args.mcount[i] = this_claim;
		m_claimed += this_claim;
		if (chat)
			printf("delaysum1_forw_thr: claiming %d loops for thread %d.\n",
				this_claim, i);
	}

	if (m_claimed != Nm)
	{
		Fail2("delaysum1_forw_thr claimed %d but needed to claim %d.",
			m_claimed, Nm);
	}
	
	jf_thread1_top(delaysum1_forw_worker, NULL, &args, nthread, 0);
	
	free(args.mstart);
	free(args.mcount);
	Ok
}


// delaysum1_forw_worker()
sof delaysum1_forw_worker(void* p, cint tid, cint nthread)
{
	delaysum1_forw_worker_args * args = p;
	if (args->chat)
		printf("delaysum1_forw_worker thread %d/%d started.\n",tid,nthread);
	
	int m_loops = args->mstart[tid];
	// try to fast-forward pointers as they should be after m_loops outer loops
	// already ran.
	float *yy = args->yy + args->Ny*m_loops;
	cfloat *dd = args->dd + args->Nx*m_loops;

	for (int im = 0; im < args->mcount[tid]; ++im, yy += args->Ny)
	{
		cfloat *xk = args->p_xk;
		for (int id=0; id < args->Nx; ++id)
		{
			int delay = (int) *dd++;
			cfloat coef = *xk++;
			float *py = yy + Max(delay,0);
			cfloat *ph = args->hh + Max(-delay,0);
			cint nh_new = args->Nh - Max(-delay,0);
			cint ny_new = args->Ny - Max(delay,0);
			cint npoint = Min(nh_new, ny_new);
			VectAddScale(float, py, coef, ph, npoint)
		}
	}
	Ok
}


// delaysum1_mex_forw()
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


// delaysum1_mex_forw_thr()
static sof delaysum1_mex_forw_thr(
mxArray *plhs[], // [Ny Nm]
Cmx mx_h, // [Nh]
Cmx mx_delay, // [Nx Nm]
Cmx mx_nthread,
Cmx mx_x, // [Nx]
Cmx mx_Ny,
Cmx mx_chat)
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

	Call(mxIsScalarInt32, (mx_chat))
	cint chat = mxGetInt(mx_chat);
	
	// output
	Call(plhs[0] = mxCreateNumericMatrix,
		(Ny*Nm, 1, mxSINGLE_CLASS, mxREAL))

	Call(delaysum1_forw_thr, ( (float *) mxGetData(plhs[0]), Ny, Nm,
			mxGetPr_cfloat(mx_h), Nh,
			mxGetPr_cfloat(mx_delay), Nx,
			mxGetPr_cfloat(mx_x),
			nthread, chat))
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
	// help (if no args)
	if (nrhs < 1 || !mxIsChar(prhs[0]))
	{
		delaysum1_mex_help();
		Call(mxu_arg, (nrhs, prhs))
		Fail(Usage)
	}

	// for "check" option needed for test_all_mex
	if (nrhs == 1 && mxIsChar(prhs[0]))
		Ok

	char *arg;
	Call(arg = mxu_string, (prhs[0], "1st argument"))

	// forw
	if (Streq(arg, "delaysum1,forw"))
	{
		if (nrhs != 7 || nlhs != 1) // require but ignore chat
			Fail(Usage)
		Call(delaysum1_mex_forw, (plhs,
			prhs[1], prhs[2], prhs[3], prhs[4], prhs[5]))
	}

	// forw, threaded
	else if (Streq(arg, "delaysum1,forw,thr"))
	{
		if (nrhs != 7 || nlhs != 1)
			Fail(Usage)
		Call(delaysum1_mex_forw_thr, (plhs,
			prhs[1], prhs[2], prhs[3], prhs[4], prhs[5], prhs[6]))
	}

	// back
	else if (Streq(arg, "delaysum1,back"))
	{
		if (nrhs != 6 || nlhs != 1) // require but ignore chat
			Fail(Usage)
		Call(delaysum1_mex_back, (plhs,
			prhs[1], prhs[2], prhs[3], prhs[4]))
	}

	// back, threaded
	else if (Streq(arg, "delaysum1,back,thr"))
	{
		if (nrhs != 6 || nlhs != 1)
			Fail(Usage)
		Call(delaysum1_mex_back_thr, (plhs,
			prhs[1], prhs[2], prhs[3], prhs[4], prhs[5]))
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
	if (!nlhs && !nrhs)
	{
		delaysum1_mex_help();
		return;
	}
	if (!delaysum1_mex(nlhs, plhs, nrhs, prhs))
		mexErrMsgTxt("delaysum1_mex");
}
#endif
