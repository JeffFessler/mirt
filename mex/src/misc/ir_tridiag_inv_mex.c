// ir_tridiag_inv_mex.c
// Matlab mex file for tridiagonal solver
// T * x = y
// assumes real tridiagonal matrix T, rhs y is formulated as NxM matrix,
// computes M columns of x in parallel
// user chooses nthreads
//
// 2013-11-09 Mai Le, initial version
// 2016-01-25 Mai Le and JF, updates

// #include "mex.h"
#include "def,mexarg.h"
#include "defs-env.h"
#include "jf,mex,def.h"
#include "pthread.h"

#define Usage "usage error. see above"
// #define MAX_THREADS 32


static void tridiag_inv_mex_help(void)
{
	printf("\n\
	Usage for ir_tridiag_inv_mex: \n\
		output = ir_tridiag_inv_mex(subdiag, diagvals, supdiag, rhs, nthread) \n\
	\n\
		T * x = y \n\
	\n\
		subdiag: (single, real) [N-1 1] -1st subdiagonal values of T \n\
		diagvals: (single, real) [N 1] diagonal values of T \n\
		supdiag: (single, real) [N-1 1] 1st diagonal values of T \n\
		rhs: (single) [N M] rhs of inverse problem, y \n\
		nthread: (int32) # of threads \n\
		output: (single) [N M] \n\
	\n");
}

struct thread_data
{
	int thread_id;
	int block_size;
	int num_blocks_for_me;
	cfloat *subdiag_ptr;
	cfloat *diagvals_ptr;
	cfloat *supdiag_ptr;
	cfloat *rhsr_ptr;
	cfloat *rhsi_ptr;
	float *outr_ptr;
	float *outi_ptr;
};


// tridiag_inv()
// tridiag solver over one N x N block of real values
static void tridiag_inv(
cfloat *a, // subdiag
cfloat *b, // diag
cfloat *c, // supdiag
cfloat *d, // rhs
cint N,
float *x, // out
float *new_c, // work
float *new_d) // work
{
	new_c[0] = c[0] / b[0];
	new_d[0] = d[0] / b[0];

	for (int ii = 1; ii <= N-2; ii++)
	{
		cfloat a_prev = a[ii - 1];
		cfloat new_c_prev = new_c[ii - 1];
		new_c[ii] = c[ii] / (b[ii] - new_c_prev * a_prev);
		new_d[ii] = (d[ii] - new_d[ii - 1] * a_prev)
				/ (b[ii] - new_c_prev * a_prev);
	}
	new_d[N - 1] = (d[N - 1] - new_d[N - 2] * (a[N - 2]))
			/ (b[N - 1] - new_c[N - 2] * (a[N - 2]));

	x[N - 1] = new_d[N - 1];
	for (int ii = N-2; ii >= 0; ii--)
		x[ii] = new_d[ii] - new_c[ii] * x[ii + 1];
}


// tridiag_inv_loop_thr()
// each thread loops over tridiag_inv
static sof tridiag_inv_loop_thr(void *threadarg)
{
	struct thread_data *my_data = (struct thread_data *) threadarg;
	cint N = my_data -> block_size;

	float new_c[N-1]; // work space
	float new_d[N];
//	float *new_c;
//	float *new_d;
//	Note1("N = %d", N)
//	Note1("num_blocks_for_me = %d", num_runs)
//	Mem0(new_c, N - 1)
//	Mem0(new_d, N)

	cfloat *a = my_data -> subdiag_ptr;
	cfloat *b = my_data -> diagvals_ptr;
	cfloat *c = my_data -> supdiag_ptr;
	cfloat *dr = my_data -> rhsr_ptr;
	cfloat *di = my_data -> rhsi_ptr;
	float *xr = my_data -> outr_ptr; // output
	float *xi = my_data -> outi_ptr;
	cint num_runs = my_data -> num_blocks_for_me;

	for (int ii = 0; ii < num_runs; ii++)
	{
		tridiag_inv(a, b, c, dr, N, xr, new_c, new_d);
		dr += N;
		xr += N;
		if (xi != NULL)
		{
			tridiag_inv(a, b, c, di, N, xi, new_c, new_d);
			di += N;
			xi += N;
		}
	}

//	Free0(new_c)
//	Free0(new_d)

	Ok
}


// tridiag_inv_loop_thr_void()
static void *tridiag_inv_loop_thr_void(void *threadarg)
{
	if (!tridiag_inv_loop_thr(threadarg))
		Warn("Failure with tridiag_inv_thread()")
	return NULL;
}


// check_types_and_sizes()
static sof check_types_and_sizes(Const mxArray *prhs[])
{
	cint Nsub = mxGetM(prhs[0]);
	cint Msub = mxGetN(prhs[0]);
	cint Ndiag = mxGetM(prhs[1]);
	cint Mdiag = mxGetN(prhs[1]);
	cint Nsup = mxGetM(prhs[2]);
	cint Msup = mxGetN(prhs[2]);
	cint Nrhs = mxGetM(prhs[3]);

	if (!((Nsub == Nrhs - 1) && (Msub == 1)) && !((Nsub == 1) && (Msub == Nrhs - 1)))
		Fail("subdiag size [%d %d] does not match rhs length of %d", Nsub, Msub, Nrhs)
	if (!((Ndiag == Nrhs) && (Mdiag == 1)) && !((Ndiag == 1) && (Mdiag == Nrhs)))
		Fail("diag size [%d %d] does not match rhs length of %d", Ndiag, Mdiag, Nrhs)
	if (!((Nsup == Nrhs - 1) && (Msup == 1)) && !((Nsup == 1) && (Msup == Nrhs - 1)))
		Fail("supdiag size [%d %d] does not match rhs length of %d", Nsup, Msup, Nrhs)

	if (mxIsComplex(prhs[0]))
		Fail("subdiag cannot be complex")
	if (!mxIsClass(prhs[0], "single"))
		Fail("subdiag must be single")
	if (mxIsComplex(prhs[1]))
		Fail("diag cannot be complex")
	if (!mxIsClass(prhs[1], "single"))
		Fail("diag must be single")
	if (mxIsComplex(prhs[2]))
		Fail("supdiag cannot be complex")
	if (!mxIsClass(prhs[2], "single"))
		Fail("supdiag must be single")
	if (!mxIsClass(prhs[3], "single"))
		Fail("rhs must be single")
	if (!mxIsScalarInt32(prhs[4]))
		Fail("nthread must be int32")

	Ok
}


// currently assume all inputs real
// wrapper function for thread
static sof tridiag_inv_mex_thr(
cfloat *subdiag_ptr, cfloat *diagvals_ptr, cfloat *supdiag_ptr,
cfloat *rhs_real_ptr, cfloat *rhs_imag_ptr,
cint block_size, cint nblocks,
float *out_real_ptr,
float *out_imag_ptr,
cint nthreads)
{
	struct thread_data *thread_data_array;
	Mem0(thread_data_array, nthreads)
//	struct thread_data thread_data_array[nthreads];
	int blocks_per_thread[nthreads];
	int cum_blocks[nthreads + 1];

	pthread_attr_t attr;
	pthread_t threads[nthreads];
	pthread_attr_init(&attr);
	pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);

	int remainder_blocks = nblocks - (nblocks/nthreads)*nthreads;
	cum_blocks[0] = 0;
	for (int th_ndx = 0; th_ndx < nthreads; th_ndx++)
	{
		blocks_per_thread[th_ndx] = nblocks/nthreads;
		if (th_ndx < remainder_blocks)
			blocks_per_thread[th_ndx]++;
		cum_blocks[th_ndx + 1] = cum_blocks[th_ndx] + blocks_per_thread[th_ndx];
	}

	for (int th_id = 0; th_id < nthreads; th_id++)
	{
		thread_data_array[th_id].thread_id = th_id;
		thread_data_array[th_id].block_size = block_size;
//		Note("block_size[%d] = %d", th_id, thread_data_array[th_id].block_size)
		thread_data_array[th_id].subdiag_ptr = subdiag_ptr;
		thread_data_array[th_id].diagvals_ptr = diagvals_ptr;
		thread_data_array[th_id].supdiag_ptr = supdiag_ptr;
		thread_data_array[th_id].rhsr_ptr = rhs_real_ptr + cum_blocks[th_id] * block_size;
		thread_data_array[th_id].outr_ptr = out_real_ptr + cum_blocks[th_id] * block_size;
		if (rhs_imag_ptr != NULL) {
			thread_data_array[th_id].rhsi_ptr =
			rhs_imag_ptr + cum_blocks[th_id] * block_size;
			thread_data_array[th_id].outi_ptr =
			out_imag_ptr + cum_blocks[th_id] * block_size;
		} else {
			thread_data_array[th_id].rhsi_ptr = NULL;
			thread_data_array[th_id].outi_ptr = NULL;
		}
		thread_data_array[th_id].num_blocks_for_me = blocks_per_thread[th_id];

		cint rc = pthread_create(&threads[th_id], &attr,
			tridiag_inv_loop_thr_void,
			(void *) (thread_data_array + th_id));
		if (rc) Fail("pthread_create")
	}

	for (int it = 0; it < nthreads; it++)
	{
		if (pthread_join(threads[it], NULL))
			Fail1("pthread_join %d failed", it)

	}

	Free0(thread_data_array)
	Ok
}


// intermediate GateWay routine
static sof tridiag_inv_mex_gw(
int nlhs, mxArray *plhs[],
int nrhs, Const mxArray *prhs[])
{
	// "check" for test_all_mex
	if (nlhs == 0 && nrhs == 1 && mxIsChar(prhs[0]))
		Ok

	if (nlhs != 1 || nrhs != 5)
	{
		Warn("Incorrect number of inputs or outputs.")
		tridiag_inv_mex_help();
		Fail(Usage)
	}

	Call(check_types_and_sizes, (prhs))

	cfloat *sub = (cfloat *) mxGetData(prhs[0]); // sub-diagonal
	cfloat *diag = (cfloat *) mxGetData(prhs[1]); // diagonal 1xN
	cfloat *sup = (cfloat *) mxGetData(prhs[2]); // sup-diagonal
	cfloat *rhs = (cfloat *) mxGetData(prhs[3]);
	cint nthread = ((cint *) mxGetData(prhs[4]))[0];
	size_t N = mxGetM(prhs[3]);	// size of tridiag matrix
	size_t M = mxGetN(prhs[3]);	// numcols of rhs matrix

#if 0
	if (nthread > MAX_THREADS)
	{
		printf("nthread:%d > MAX_THREADS: %d, truncated to %d \n",
			nthread, MAX_THREADS, MAX_THREADS);
		nthread = MAX_THREADS;
	}
#endif

	float *rhs_imag;
	// output
	if (mxIsComplex(prhs[3])) {
		rhs_imag = (float *) mxGetImagData(prhs[3]);
		plhs[0] = mxCreateNumericMatrix(N, M, mxSINGLE_CLASS, mxCOMPLEX);
	} else {
		rhs_imag = NULL;
		plhs[0] = mxCreateNumericMatrix(N, M, mxSINGLE_CLASS, mxREAL);
	}

	float *x_real = (float *) mxGetData(plhs[0]);
	float *x_imag = (float *) mxGetImagData(plhs[0]); // NULL when rhs is real

	if ((rhs_imag == NULL) ^ (x_imag == NULL)) // debug check
		Fail("problem: only one NULL vec, should be both or neither")

	Call(tridiag_inv_mex_thr,
		(sub, diag, sup, rhs, rhs_imag, N, M, x_real, x_imag, nthread))
	Ok
}


// gateway routine
#if defined(Need_ir_tridiag_inv_mex_gateway)
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	if (!nlhs && !nrhs) {
		tridiag_inv_mex_help();
		return;
	}
	if (!tridiag_inv_mex_gw(nlhs, plhs, nrhs, prhs))
		mexErrMsgTxt("ir_tridiag_inv_mex");
}
#endif
