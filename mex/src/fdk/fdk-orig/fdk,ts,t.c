// fdk,ts,t.c
// Threaded versions of FDK back-projection
// For detector index (t,s).
// Copyright 2008-10-09, Jeff Fessler, University of Michigan

#include "defs-env.h"
#include "def,fdk.h"
#include "jf,thread1.h"
#include "jf,time.h"

// fdk_ts_put_view()
// put projection view into center of view_work which has 1 pixel all around
void fdk_ts_put_view(
float *view_work, // [(nt+2) (ns+2)], assume outer border already 0
cfloat *proj, // [nt ns]
cint ns,
cint nt)
{
	for (int is=0; is < ns; ++is)
		Bcopy(proj + is*nt, view_work + (is+1) * (nt+2) + 1, nt)
}


typedef struct
{
	float *image; // [nz nx ny] <- trick!
	const cbct_ig *ig; // image geometry
	const cbct_cg *cg; // cone-beam CT system geometry
	int na; // # of views
	cfloat *proj; // [nt ns na] <- trick! projection views
	cdouble *beta; // [na] source angles [radians]
	int view_size; // (nt+2) * (ns+2)
	float *view_work; // [nt+2 ns+2 nthread] work space for views
} fdk_ts_s;


// fdk_ts_back_init()
// interface routine for threaded versions
static sof fdk_ts_back_init(void *in, cint id, cint nthread)
{
	fdk_ts_s *pa = (fdk_ts_s *) in;
	const cbct_ig *ig = pa->ig;
	const cbct_cg *cg = pa->cg;
	cint na = pa->na;
	cfloat *proj = pa->proj;
	cdouble *beta = pa->beta;
	(void) nthread;

	float *view_work = pa->view_work + id * pa->view_size;

	cint ns = cg->ns;
	cint nt = cg->nt;
	for (int ia=0; ia < na; ++ia, proj += ns * nt) // each view
	{
		fdk_ts_put_view(view_work, proj, ns, nt);
		Call(fdk_ts_back1, (pa->image,
			ig->nx, ig->ny, ig->nz,
			ig->dx, ig->dy, ig->dz,
			ig->offset_x, ig->offset_y, ig->offset_z,
			ig->mask2, id + 1, // each thread does some voxels only
			cg->dso, cg->dsd, cg->dfs,
			ns, nt,
			cg->ds, cg->dt, cg->offset_s, cg->offset_t,
			beta[ia], view_work))
	}
	Ok
}


// fdk_ts_back_t()
// entry point for threaded FDK back-projector
sof fdk_ts_back_t(
float *image, // [nz nx ny] <- trick!
const cbct_ig *ig,
const cbct_cg *cg,
cint na, // # of views
cfloat *proj, // [nt ns na] <- trick! projection views
cdouble *beta, // [na] source angles [radians]
cint nthread, // # of threads
cint chat)
{
	fdk_ts_s st;

	cint view_size = (cg->nt+2) * (cg->ns+2);
	float *view_work;
	Mem0(view_work, view_size * nthread)

#define put(arg) st.arg = arg;
	put(image)
	put(ig)
	put(cg)
	put(na)
	put(proj)
	put(beta)
	put(view_size)
	put(view_work)
#undef put

	Bzero(image, ig->nx * ig->ny * ig->nz) // initialize image volume to 0

	Call(jf_time, ("start"))

	Call(jf_thread1_top, (fdk_ts_back_init,
                NULL /* wrap up */, &st, nthread, Chat))

	if (chat)
		Call(jf_time, ("report fdk,ts,back"))

        Ok
}
