// ir_shrink1.c
// Iterative 1D shrinkage algorithm
//
// min_x 1/2 (y - x)^2 + reg * pot(x)
// user must provide tabulated dpot to be used with linear interpolation.
//
// Copyright 2012-06-25, Jeff Fessler, University of Michigan

#include "ir_shrink1.h"
#include "jf,thread1.h"


// ir_shrink1()
// returns solution to
// min_x 1/2 (y - x)^2 + reg * pot(x)
static float ir_shrink1(
cfloat y, // data value
cfloat reg, // regularization parameter
cfloat *table, // [K] dpot([1:K] * dt)
cfloat dt, // table spacing
cint K, // # of points in table
cint niter) // max # of iterations
{
	cfloat ay = Abs(y);
	cfloat sy = (y > 0) ? 1. : -1.;

	cfloat d1 = table[0]; // dpot(dt)
	cfloat b1 = dt + reg * d1;
	if (ay < b1) // in quadratic regime for |y| near 0
		return y / (1 + reg * d1 / dt);

	cfloat dK = table[K-1]; // dpot(K dt)
	cfloat bK = K * dt + reg * dK;
	if (ay >= bK)
		return sy * (ay - reg * dK);

	float x = ay; // trick: work with y > 0 until return

//	Note4("y=%g b1=%g bK=%g niter=%d", y, b1, bK, niter)

	for (int ii=0; ii < niter; ++ii)
	{
		cint k = floor(x / dt);
		if (k == 0 || k >= K)
		{
			Note1("y=%g bug", y)
			return y;
		}

		cfloat slope = (table[k] - table[k-1]) / dt;
		cfloat dpot = table[k-1] + (x - k * dt) * slope;
		x = ay / (1. + reg * dpot / x);
//		Note2("y=%g x=%g", ay, x)
	}

	return sy * x;
}


typedef struct
{
	float *x;
	cfloat *y;
	int N;
	cfloat *reg;
	cfloat *table;
	float dt;
	int K;
	int niter;
	float thresh;
	int chat;
} ir_shrink1_s;


// ir_shrink1_work()
static sof ir_shrink1_work(void *p, cint id, cint nthread)
{
	ir_shrink1_s* st = (ir_shrink1_s *) p;

	cint N = st->N;
	cint n_per_thread = (int) ceil(N / (double) nthread);
	cint n_full = N - (n_per_thread - 1) * nthread;
	int jj0, jj1;

	if (id < n_full) // fully loaded threads
        {
		jj0 = id * n_per_thread;
                jj1 = jj0 + n_per_thread;
        }
        else // partially loaded threads each have one less job
        {
		jj0 = n_full * n_per_thread
                        + (id - n_full) * (n_per_thread - 1);
 		jj1 = jj0 + (n_per_thread - 1);
        }

	if (st->chat)
		Note4("id=%d of %d: jj0=%d jj1=%d", id, nthread, jj0, jj1)

	float *x = st->x + jj0;
	cfloat *y = st->y + jj0;
	cfloat *reg = st->reg + jj0;
	cfloat *table = st->table;
	cfloat dt = st->dt;
	cint K = st->K;
	cint niter = st->niter;

	for (int jj=jj0; jj < jj1; ++jj)
		*x++ = ir_shrink1(*y++, *reg++, table, dt, K, niter);

	Ok
}


// ir_shrink1_p()
sof ir_shrink1_p(
float *x, // [N]
cfloat *y, // [N]
cint N,
cfloat *reg, // [N]
cfloat *table, // [K]
cfloat dt,
cint K,
cfloat thresh,
cint niter,
cint nthread,
cint chat)
{
	ir_shrink1_s st;

	st.x = x;
	st.y = y;
	st.N = N;
	st.reg = reg;
	st.table = table;
	st.dt = dt;
	st.K = K;
	st.thresh = thresh;
	st.niter = niter;
	st.chat = chat;

	if (chat)
		Note5("N=%d K=%d dt=%g thresh=%g niter=%d",
				N, K, dt, thresh, niter)

	Call(jf_thread1_top,
		(ir_shrink1_work, NULL, &st, nthread, Chat))
	Ok
}
