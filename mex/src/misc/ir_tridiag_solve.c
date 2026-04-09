// ir_tridiag_solve.c
//
// UNDER DEVELOPMENT!
// Tridiagonal solver: T * x = y so x = T \ y
//
// Assumes real NxN tridiagonal matrix T.
// RHS y is NxM matrix.
// Computes M columns of x in parallel; user chooses nthreads
//
// 2013-11-09 Mai Le, initial version
// 2016-01-25 Mai Le and JF, updates
// 2016-02-03 Mai Le and JF, handle M different tridiagonal matrices

#include "defs-env.h"
#include "pthread.h"


// tridiag_solve1()
// tridiag solver for one NxN block
// https://en.wikipedia.org/wiki/Tridiagonal_matrix_algorithm
static void tridiag_solve1(
cfloat *a, // [N-1] subdiag
cfloat *b, // [N] diag
cfloat *c, // [N-1] supdiag
cfloat *d, // [N] rhs
cint N,
float *x, // [N] out
float *new_c, // [N-1] work
float *new_d) // [N] work
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



// tridiag_solve_many_rhs()
// same matrix but M RHS vectors
void tridiag_solve_many_rhs(
cfloat *a, // [N-1] subdiag
cfloat *b, // [N] diag
cfloat *c, // [N-1] supdiag
cfloat *dd, // [N M] rhs
cint N,
float *xx, // [N M] out
float *new_c, // [N-1] work
float *new_d, // [N] work
cint M)
{
	for (int mm = 0; mm < M; mm++, dd += N, xx += N)
		tridiag_solve1(a, b, c, dd, N, xx, new_c, new_d);
}


// tridiag_solve_many_sys()
// M different matrices and RHS vectors
void tridiag_solve_many_sys(
cfloat *aa, // [N-1 M] subdiag
cfloat *bb, // [N M] diag
cfloat *cc, // [N-1 M] supdiag
cfloat *dd, // [N M] rhs
cint N,
float *xx, // [N M] out
float *new_c, // [N-1] work
float *new_d, // [N] work
cint M)
{
	for (int mm = 0; mm < M; mm++,
		aa += N, bb += N, cc += N, dd += N, xx += N)
		tridiag_solve1(aa, bb, cc, dd, N, xx, new_c, new_d);
}



// tridiag_solve_loop_thr_void()
// todo: cut this?  seems useless
static void *tridiag_solve_loop_thr_void(void *threadarg)
{
	if (!tridiag_solve_loop_thr(threadarg))
		Warn("Failure with tridiag_solve_thread()")
	return NULL;
}


// todo: try openMP instead of all this?

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


// tridiag_solve_loop_thr()
// each thread loops over tridiag_solve
static sof tridiag_solve_loop_thr(void *threadarg)
{
	struct thread_data *my_data = (struct thread_data *) threadarg;
	cint N = my_data -> block_size;

	float *new_c;
	float *new_d;
	Mem0(new_c, N - 1)
	Mem0(new_d, N)

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
		tridiag_solve(a, b, c, dr, N, xr, new_c, new_d);
		dr += N;
		xr += N;
		if (xi != NULL)
		{
			tridiag_solve(a, b, c, di, N, xi, new_c, new_d);
			di += N;
			xi += N;
		}
		a += N - 1;
		b += N;
		c += N - 1;
	}

	Free0(new_c)
	Free0(new_d)
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
		thread_data_array[th_id].subdiag_ptr = subdiag_ptr + cum_blocks[th_id] * (block_size - 1);
		thread_data_array[th_id].diagvals_ptr = diagvals_ptr + cum_blocks[th_id] * block_size;
		thread_data_array[th_id].supdiag_ptr = supdiag_ptr + cum_blocks[th_id] * (block_size - 1);\
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

	for (int th_id = 0; th_id < nthreads; th_id++) {
		cint rc = pthread_join(threads[th_id], NULL);
		if (rc) {
			printf("ERROR; return code from pthread_join() is %d\n", rc);
			exit(-1);
		}
	}

	Free0(thread_data_array)
	Ok
}


// intermediate GateWay routine
static sof tridiag_inv_mex_gw(
int nlhs, mxArray *plhs[],
int nrhs, Const mxArray *prhs[])
{
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

	if (nthread > M)
		nthread = M;

	float *rhs_imag;
	// output
	if (mxIsComplex(prhs[3])) {
		rhs_imag = (float *) mxGetImagData(prhs[3]);
		plhs[0] = mxCreateNumericMatrix(N * M, 1, mxSINGLE_CLASS, mxCOMPLEX);
	} else {
		rhs_imag = NULL;
		plhs[0] = mxCreateNumericMatrix(N * M, 1, mxSINGLE_CLASS, mxREAL);
	}

	float *x_real = (float *) mxGetData(plhs[0]);
	float *x_imag = (float *) mxGetImagData(plhs[0]); // NULL when rhs is real

	if ((rhs_imag == NULL) ^ (x_imag == NULL)) // debug check
		Fail("problem: only one NULL vec, should be both or neither")

	Call(tridiag_inv_mex_thr,
		(sub, diag, sup, rhs, rhs_imag, N, M, x_real, x_imag, nthread))
	Ok
}
