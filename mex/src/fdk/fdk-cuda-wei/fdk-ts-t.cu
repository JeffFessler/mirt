// fdk-ts-t.cu
// Threaded versions of FDK back-projection
// For detector index (t,s).
// Copyright 2008-10-09, Jeff Fessler, University of Michigan

#include "jf-cuda.h"
#include "def,fdk.h"
#include "jf,thread1.h"
#include "fdk-gpu.h"


typedef struct {
	float *image; // [nz nx ny] <- trick!
	const cbct_ig *ig; // image geometry
	const cbct_cg *cg; // cone-beam CT system geometry
	int na; // # of views
	cfloat *proj; // [nt ns na] <- trick! projection views
	cdouble *beta; // [na] source angles [radians]
} fdk_ts_s;


//
// fdk_ts_back_init()
// interface routine for threaded versions
//
static sof fdk_ts_back_init(void *in, cint id, cint nthread)
{
	fdk_ts_s *pa = (fdk_ts_s *) in;
	const cbct_ig *ig = pa->ig;
	const cbct_cg *cg = pa->cg;
	cint na = pa->na;
	cfloat *proj = pa->proj;
	cdouble *beta = pa->beta;
	cint nst = cg->ns * cg->nt;

	if (nthread > 1) Fail("only nthread=1 done")

#ifdef fdk_gpu
	cint nxyz = ig->nx * ig->ny * ig->nz;
	float *dev_img;
	jf_gpu_malloc(dev_img, nxyz) // image memory on device
	jf_gpu_memset(dev_img, 0, nxyz) // initialize device image to 0

	float *dev_proj;
	jf_gpu_malloc(dev_proj, nst) // one projection view on device
 

	byte *dev_mask2;
	cint nxy = ig->nx * ig->ny;
	jf_gpu_malloc(dev_mask2, nxy) // 2D mask
	jf_gpu_put(dev_mask2, ig->mask2, nxy)

// NEW: values of s for each (x,y) set
	float *dev_sval;			// s values for each x,y pair on device
	jf_gpu_malloc(dev_sval, nxy*3)	// 2D img
	jf_gpu_memset(dev_sval, 0, nxy)	// initialize values to 0

//TEST:
/*	FILE * fptr;
	fptr = fopen("test.txt","w");
	for (int i=0; i<(ig->nx); i++){
		for (int j=0; j<(ig->ny); j++){
			fprintf(fptr,"%x ", ig->mask2[i*(ig->nx) + j]);
		}
		fprintf(fptr,"\n");
	}
*/


#endif

	for (int ia=0; ia < na; ++ia, proj += nst) { // each view

#ifdef fdk_gpu
		// copy this view to gpu
		jf_gpu_put(dev_proj, proj, nst)

#else
		float *dev_img = pa->image; // already zeroed
		cfloat *dev_proj = proj;
		cbyte *dev_mask2 = ig->mask2;
		
#endif

		if (!fdk_ts_back1_gpu(dev_img,
#ifdef fdk_gpu		// NEW
			dev_sval,
#endif
			ig->nx, ig->ny, ig->nz,
			ig->dx, ig->dy, ig->dz,
			ig->offset_x, ig->offset_y, ig->offset_z,
			dev_mask2, id + 1, // each thread does some voxels only
			cg->dso, cg->dsd, cg->dfs,
			cg->ns, cg->nt,
			cg->ds, cg->dt, cg->offset_s, cg->offset_t,
			dev_proj, beta[ia]))
			Fail("fdk_ts_back1_gpu()")
	}

#ifdef fdk_gpu
//	Note("Copying image to host")
	jf_gpu_get(pa->image, dev_img, nxyz) // caution: works only for 1 thread

//	Note("freeing dev_img memory")
	jf_gpu_free(dev_img)

//	Note("freeing dev_proj memory\n")
	jf_gpu_free(dev_proj)

// NEW: freeing dev_sval
// Note("freeing dev_sval memory\n")
	jf_gpu_free(dev_sval)

// NEW: freeing the mask
	jf_gpu_free(dev_mask2)

#endif

	Ok
}


//
// fdk_ts_back_t()
// entry point for threaded FDK back-projector
//
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
#define put(arg) st.arg = arg;
	put(image)
	put(ig)
	put(cg)
	put(na)
	put(proj)
	put(beta)
#undef put

// TESTING
// printf("MY FILE COOLNESS! \n");

	Bzero(image, ig->nx * ig->ny * ig->nz) // initialize image volume to 0

	Call(jf_thread1_top, (fdk_ts_back_init,
                NULL /* wrap up */, &st, nthread, Chat))
        Ok
}
