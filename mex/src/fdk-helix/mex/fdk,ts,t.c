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
	cdouble *source_zs; // [na] source z-location by DK
	int na1_half; // [1] na1_half by DK
	int cone_par; // [1] cone_par by DK
	int w3d; // [1] w3d by DK
	double source_dz_per_rad; // [1] source_dz_per_rad by DK
	float pitch; // [1] pitch by DK
	float *pow_tab; // [nt/2] by DK
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
	cdouble *source_zs = pa->source_zs; // DK
	cint na1_half = pa->na1_half; // DK
	cint cone_par = pa->cone_par; // DK
	cint w3d = pa->w3d; // DK
	cdouble source_dz_per_rad = pa->source_dz_per_rad; // DK
	cfloat pitch = pa->pitch; // DK
	cfloat *pow_tab = pa->pow_tab; // DK
	(void) nthread;

	float *view_work = pa->view_work + id * pa->view_size;

	cint ns = cg->ns;
	cint nt = cg->nt;
	cfloat dz = ig->dz; // DK
	cfloat wz = (ig->nz-1)/2. + ig->offset_z; // DK

	cfloat source_zs_min = source_zs[0]; // DK
	cfloat source_zs_max = source_zs[na-1]; // DK

	// todo na -> na/2
	for (int ia=0; ia < na; ++ia, proj += ns * nt) // each view
	{
		// DK: todo!!: z-slices to be updated for [beta - pi, beta + pi]
		int iz_min_beta = 0;
		int iz_max_beta = ig->nz;
		if (pitch)
		{
			iz_min_beta = Max(Ceilf(wz + source_zs[Max(ia - na1_half, 0)]/dz), 0); 
			iz_max_beta = Min(Ceilf(wz + source_zs[Min(ia + na1_half, na-1)]/dz), ig->nz);
		}

		fdk_ts_put_view(view_work, proj, ns, nt);
		Call(fdk_ts_back1, (pa->image,
			ig->nx, ig->ny, ig->nz,
			ig->dx, ig->dy, ig->dz,
			ig->offset_x, ig->offset_y, ig->offset_z,
			ig->mask2, id + 1, // each thread does some voxels only
			cg->dso, cg->dsd, cg->dfs,
			ns, nt,
			cg->ds, cg->dt, cg->offset_s, cg->offset_t,
			beta[ia], 
			source_zs[ia], // DK
			source_zs_min, source_zs_max, // DK
			iz_min_beta, iz_max_beta, // DK
			cone_par, w3d, source_dz_per_rad, pitch, pow_tab, // DK
			view_work))
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
cdouble *source_zs, // [na] source z-location by DK
cint na1_half, // [1] na1_half by DK
cint cone_par, // [1] cone_par by DK
cint w3d, // [1] w3d by DK
cdouble source_dz_per_rad, // [1] source_dz_per_rad by DK
cfloat pitch, // [1] pitch by DK
cint nthread, // # of threads
cint chat)
{
	fdk_ts_s st;

	cint view_size = (cg->nt+2) * (cg->ns+2);
	float *view_work;
	Mem0(view_work, view_size * nthread)

	// DK
	float *pow_tab;
	cint ntab = pitch ? Ceilf(cg->nt/2) : (Ceilf(cg->nt/2)*Ceilf(ig->nz/2));
	Mem0(pow_tab, ntab);
	cfloat wz = (ig->nz-1)/2. + ig->offset_z;
	Call(fdk_w3d_pow_tab_init, (pow_tab,  
			cg->nt, ig->nz, wz, ig->dz, cg->dsd, 
			w3d, pitch));

#define put(arg) st.arg = arg;
	put(image)
	put(ig)
	put(cg)
	put(na)
	put(proj)
	put(beta)
	put(source_zs) // DK
	put(na1_half) // DK
	put(cone_par) // DK
	put(w3d) // DK
	put(source_dz_per_rad) // DK
	put(pitch) // DK
	put(pow_tab) // DK
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
