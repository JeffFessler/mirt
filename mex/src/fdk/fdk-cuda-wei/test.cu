// fdk-mex-cuda.cu
// Matlab mex gateway for CUDA/GPU version FDK.
//
// Copyright 2005-6-27, Jeff Fessler, University of Michigan

#include "jf,mex,def.h"
#include "def,fdk.h"

#define Usage	"usage error. see above"

static void fdk_mex_help(void)
{
	printf("Usage for FDK mex routines for\n\
cone-beam weighted 3d backprojection of 2d projection view(s).\n\
");

	fdk_st_help();
	fdk_ts_help();
}


static sof fdk_mex_back(
mxArray *plhs[],
Cmx mx_nxyz, // [3] image dimensions
Cmx mx_dxyz, // [3] voxel size
Cmx mx_cxyz, // [3] image center offset
Cmx mx_mask, // [nx ny] 2D mask
Cmx mx_dso, // distance source to origin (iso-center)
Cmx mx_dsd, // distance source to detector
Cmx mx_dfs, // 0 for 3rd gen CT or inf for flat
Cmx mx_dst, // [2] ds dt
Cmx mx_ost, // [2] offset_s offset_t
Cmx mx_proj, // [ns nt] or [nt ns] or [ns nt na] or [nt ns na]
Cmx mx_beta, // [1] or [na]
Cmx mx_nthread,
ctruf is_st)
{
	int ns, nt, na;

	Call(mxIsInt32, (mx_nxyz))
	if (3 != mxGetM(mx_nxyz) * mxGetN(mx_nxyz))
		Fail("mx_nxyz must have 3 elements")
	cint nx = ((cint *) mxGetPr(mx_nxyz))[0];
	cint ny = ((cint *) mxGetPr(mx_nxyz))[1];
	cint nz = ((cint *) mxGetPr(mx_nxyz))[2];

	Call(mxIsRealDouble, (mx_dxyz))
	if (3 != mxGetM(mx_dxyz) * mxGetN(mx_dxyz))
		Fail("mx_dxyz must have 3 elements")

	Call(mxIsRealDouble, (mx_cxyz))
	if (3 != mxGetM(mx_cxyz) * mxGetN(mx_cxyz))
		Fail("mx_cxyz must have 3 elements")

	Call(mxIsUint8, (mx_mask))
	if (2 != mxGetNumberOfDimensions(mx_mask))
		Fail("mask must be 2d")
	if (nx != (int) mxGetM(mx_mask) || ny != (int) mxGetN(mx_mask))
		Fail("mask size not nx by ny?")

	Call(mxIsScalarDouble, (mx_dso))
	Call(mxIsScalarDouble, (mx_dsd))
	Call(mxIsScalarDouble, (mx_dfs))

	Call(mxIsRealDouble, (mx_dst))
	if (2 != mxGetM(mx_dst) * mxGetN(mx_dst))
		Fail("mx_dst must have 2 elements")

	Call(mxIsRealDouble, (mx_ost))
	if (2 != mxGetM(mx_ost) * mxGetN(mx_ost))
		Fail("mx_ost must have 2 elements")

	Call(mxIsRealSingle, (mx_proj))
	if (2 == mxGetNumberOfDimensions(mx_proj)) {
		na = 1;
		Call(mxIsScalarDouble, (mx_beta))
	}
	else if (3 == mxGetNumberOfDimensions(mx_proj)) {
		na = mxGetDimensions(mx_proj)[2];
		Call(mxIsRealDouble, (mx_beta))
		int na_beta = (int) mxGetM(mx_beta) * mxGetN(mx_beta);
		if (na != na_beta)
			Fail2("beta length=%d vs na=size(proj,3)=%d", na_beta, na)
	}
	else
		Fail("mx_proj must be 2d or 3d")

	Call(mxIsScalarInt32, (mx_nthread))
	int nthread = mxGetInt(mx_nthread);
	if (nthread > 1) {
		Warn1("ignoring nthread=%d, using 1 instead", nthread)
		nthread = 1;
	}

	if (is_st) {
		mwSize odim[3] = {nx, ny, nz}; // usual order
		Call(plhs[0] = mxCreateNumericArray, (3, odim, mxSINGLE_CLASS, mxREAL))
		ns = mxGetDimensions(mx_proj)[0];
		nt = mxGetDimensions(mx_proj)[1];

	}
	else {
		mwSize odim[3] = {nz, nx, ny}; // trick: zxy order!
		Call(plhs[0] = mxCreateNumericArray, (3, odim, mxSINGLE_CLASS, mxREAL))
		nt = mxGetDimensions(mx_proj)[0];
		ns = mxGetDimensions(mx_proj)[1];
		// Note2("nt=%d ns=%d", nt, ns)
	}

	// populate system and image geometry structures
	cbct_ig ig_, *ig = &ig_;
	cbct_cg cg_, *cg = &cg_;
	cint chat = 0;

	Mem0(ig->mask2, nx * ny)
	Call(cbct_mask_init, (ig->mask2, (cbyte *) mxGetData(mx_mask),
		nx, ny, nthread, chat))

	ig->nx = nx;
	ig->ny = ny;
	ig->nz = nz;
	ig->dx = mxGetPr_cdouble(mx_dxyz)[0];
	ig->dy = mxGetPr_cdouble(mx_dxyz)[1];
	ig->dz = mxGetPr_cdouble(mx_dxyz)[2];
	ig->offset_x = mxGetPr_cdouble(mx_cxyz)[0];
	ig->offset_y = mxGetPr_cdouble(mx_cxyz)[1];
	ig->offset_z = mxGetPr_cdouble(mx_cxyz)[2];

	cg->dso = mxGetDouble(mx_dso);
	cg->dsd = mxGetDouble(mx_dsd);
	cg->dfs = mxGetDouble(mx_dfs);
	cg->ns = ns;
	cg->nt = nt;
	cg->ds = mxGetPr_cdouble(mx_dst)[0];
	cg->dt = mxGetPr_cdouble(mx_dst)[1];
	cg->offset_s = mxGetPr_cdouble(mx_ost)[0];
	cg->offset_t = mxGetPr_cdouble(mx_ost)[1];

	if (is_st) {
		if (nthread != 1)
			Warn("only nthread=1 supported, ignoring")

		for (int ia=0; ia < na; ++ia)
			fdk_st_back1(
				(float *) mxGetData(plhs[0]),
				ig->nx, ig->ny, ig->nz,
				ig->dx, ig->dy, ig->dz,
				ig->offset_x, ig->offset_y, ig->offset_z,
				ig->mask2,
				cg->dso, cg->dsd, /* cg->dfs, */
				cg->ns, cg->nt,
				cg->ds, cg->dt, cg->offset_s, cg->offset_t,
				mxGetPr_cfloat(mx_proj) + ia * nt * ns,
				mxGetPr_cdouble(mx_beta)[ia],
				cg->dfs == 0 /* is_arc */);

	}
	else { // new faster way
			//TEST:
//			printf("YAY\n");

		Call(fdk_ts_back_t, (
			(float *) mxGetData(plhs[0]),
			ig, cg, na,
			mxGetPr_cfloat(mx_proj),
			mxGetPr_cdouble(mx_beta),
			nthread, chat))
	}

	Free0(ig->mask2)
	Ok
}


// intermediate gateway routine
sof fdk_mex(
int nlhs, mxArray *plhs[],
int nrhs, Const mxArray *prhs[])
{
	char *arg;

	if (nrhs <= 1 || !mxIsChar(prhs[0])) { // expect 'fdk,ts,back', ...
		fdk_mex_help();
		Call(mxu_arg, (nrhs, prhs))
		Fail(Usage)
	}

	if (nrhs == 1 && mxIsChar(prhs[0])) { // 'help'
		fdk_mex_help();
		Ok
	}

	Call(arg = mxu_string, (prhs[0], "1st argument"))

	// fdk,{st|ts},back
	if (Streq(arg, "fdk,st,back") || Streq(arg, "fdk,ts,back")) {
		ctruf is_st = Streqn(arg, "fdk,st", 6);
		if (nrhs != 13 || nlhs != 1)
			Fail(Usage)
		Call(fdk_mex_back, (plhs,
			prhs[1], prhs[2], prhs[3], prhs[4], prhs[5],
			prhs[6], prhs[7], prhs[8], prhs[9], prhs[10],
			prhs[11], prhs[12], is_st))
	}

	else
		Fail1("command '%s' unknown", arg)

	Call(mxu_string_free, (arg))
	Ok
}

#ifdef Use_fdk_mex
// gateway routine
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, Const mxArray *prhs[])
{
	if (!nlhs && !nrhs) {
		fdk_mex_help();
		return;
	}
	if (!fdk_mex(nlhs, plhs, nrhs, prhs))
		mexErrMsgTxt("fdk_mex");
}
#endif
// fdk-ts.cu
// Feldkamp aka FDK backprojection for arc/flat detector.
// For detector index (t,s).
// CUDA/GPU version
// Copyright 2005-6-27, Jeff Fessler, University of Michigan

// VERSION 5: shared memory 

#include "jf-cuda.h"
#include "def,fdk.h"
#include "fdk-gpu.h"



/* **************************************************************
BEGIN KERNEL
*************************************************************** */

//
// fdk_ts_back1_kernel()
// The FDK backprojection is *added* to the image, so the user must zero it!
//
static
#ifdef fdk_gpu
__global__
#endif
void fdk_ts_back1_kernel(
#ifndef fdk_gpu
cint ix, // ix, iy : index of x and y 
cint iy,
#endif
// NEW:
#ifdef fdk_gpu
float *s_val, // [nx ny 4] s values, w2, mag for x,y pairs (2D)
#endif
float *image, // [nz nx ny] <- trick!
int nx,
int ny,
int nz,
float dx, // voxel size
float dy, // can be negative to cause flip
float dz,
float offset_x, // image volume center offset in pixels (usually 0)
float offset_y,
float offset_z,
cbyte *mask2, // [nx ny] 2D support mask: 0, 1, ..., nthread
byte mask_id, // 1 ... nthread
float dso, // distance from source to isocenter
float dsd, // distance from source to detector
truf is_arc,
int ns, // projection view dimensions
int nt,
float ds, // horizontal ray spacing (view sample spacing)
float dt, // vertical ray spacing (view sample spacing)
float offset_s, // channel offset [pixels]
float offset_t, // vertical offset on detector [pixels]
cfloat *proj, // [nt ns] <- trick! projection view at angle beta
float beta) // source angle [radians]
{
#ifdef fdk_gpu
	// index into image array
	cint ix = blockIdx.x * blockDim.x + threadIdx.x;	// determine the index of x, y
	cint iy = blockIdx.y * blockDim.y + threadIdx.y;
#endif

// if index is out of bound
	if (ix >= nx || iy >= ny)
		return;

//NEW: ATTN: not figuring out the image location if using gpu
#ifndef fdk_gpu
	image += (ix + iy * nx) * nz;
#endif

// bound issues
	if (mask2[ix + iy*nx] != mask_id) // each thread does its part only
		return;


#if 0
	if (iy == 40 || ix == 7) return;
	for (int iz=0; iz < nz; ++iz)
		image[iz] = 7+0*iz;
	return;
#endif

	cfloat wx = (nx-1)/2. + offset_x;
	cfloat wy = (ny-1)/2. + offset_y;
	cfloat ws = (ns-1)/2. + offset_s;
// NEW:
#ifndef fdk_gpu
	cfloat wz = (nz-1)/2. + offset_z;
	cfloat wt = (nt-1)/2. + offset_t;
#endif
	cfloat sinb = sinf(beta);
	cfloat cosb = cosf(beta);

	cfloat yy = dy * (iy - wy);
	cfloat xx = dx * (ix - wx);
	cfloat xbeta = xx * cosb + yy * sinb;
	cfloat ybetas = dso - (-xx * sinb + yy * cosb);
	cfloat mag = dsd / ybetas;
	cfloat ss = is_arc ? (dsd * atan2f(xbeta, ybetas))
				: (mag * xbeta);
	cfloat ss_bin = ss / ds + ws;

// NEW:
#ifdef fdk_gpu
	cfloat w2 = is_arc ? // fan-beam image domain weighting
		(Sqr(dsd) / (Sqr(ybetas) + Sqr(xbeta))) : Sqr(mag);

	s_val[ix + iy*nx] = (float)ss_bin;
	s_val[ix + iy*nx + nx*ny] = (float)w2;
	s_val[ix + iy*nx + 2*nx*ny] = (float)mag;
#else
// index of s is "is"
	cint is = floorf(ss_bin); // index of nearest neighbor in "s"

// Check if index out of bound
	if (is < 0 || is >= ns-1) // out of FOV
		return;



	cfloat w2 = is_arc ? // fan-beam image domain weighting
		(Sqr(dsd) / (Sqr(ybetas) + Sqr(xbeta))) : Sqr(mag);

	cfloat wr = ss_bin - is; // horizontal bilinear
	cfloat wl = 1. - wr; // interpolation factors
	float *pi = image; //IMAGE!
	cfloat *pp1 = proj + is * nt;
	cfloat *pp2 = proj + (is+1) * nt;

// with ix, iy fixed, loop through all possible iz values
	for (int iz = 0; iz < nz; ++iz, ++pi) { // slice loop
		cfloat zz = dz * (iz - wz);
		cfloat tt = mag * zz;
		cfloat tt_bin = tt / dt + wt;
	// z value is used to determine index of t "it"
		cint it = floorf(tt_bin); // nearest nbr in "t"

		if (it < 0 || it >= nt-1) // out of FOV
			continue;
		else { // reconstructing the image
			cfloat wu = tt_bin - it;
			cfloat wd = 1. - wu;
			cfloat p1 = wl * pp1[it]
				+ wr * pp2[it]; // interpolate
			cfloat p2 = wl * pp1[it+1]
				+ wr * pp2[it+1]; // horizontal

			// final vertical interpolation:
			*pi += w2 * (wu * p1 + wd * p2);
		}
	}
#endif
}


/* *******************************************************************
END OF KERNEL
******************************************************************** */



/* ******************************************************************
NEW: SECOND KERNEL for z coordinate
******************************************************************** */
#ifdef fdk_gpu
__global__ void fdk_ts_back1_kernel_2(
float *s_val, // [nx ny 3] s values, w2, mag for x,y pairs (2D)
float *image, // [nz nx ny] <- trick!
int nx,
int ny,
int nz,
//float dx, // voxel size
//float dy, // can be negative to cause flip
float dz,
//float offset_x, // image volume center offset in pixels (usually 0)
//float offset_y,
float offset_z,
cbyte *mask2, // [nx ny] 2D support mask: 0, 1, ..., nthread
byte mask_id, // 1 ... nthread
//float dso, // distance from source to isocenter
//float dsd, // distance from source to detector
//truf is_arc,
int ns, // projection view dimensions
int nt,
//float ds, // horizontal ray spacing (view sample spacing)
float dt, // vertical ray spacing (view sample spacing)
//float offset_s, // channel offset [pixels]
float offset_t, // vertical offset on detector [pixels]
cfloat *proj) // [nt ns] <- trick! projection view at angle beta
//float beta) // source angle [radians]
{

/*	cint ix = fmodf((blockIdx.x*blockDim.x + threadIdx.x), nx);
	cint iy = blockIdx.y*blockDim.y + threadIdx.y;
	cint iz = floorf((blockIdx.x*blockDim.x+threadIdx.x)/nx)*blockDim.z + threadIdx.z;
*/
	// NEW 5:
	// index into image array
	cint ix = fmodf(blockIdx.x, nx);
	cint iy = blockIdx.y; 
	cint iz = threadIdx.x+512*(nz/512);	// if more than 512, iz=threadIdx.x+512
													// else iz=threadIdx.x 

	cfloat wz = (nz-1)/2. + offset_z;
	cfloat wt = (nt-1)/2. + offset_t;

	//NEW 5: use shared mem
	__shared__ float ss_bin;
	__shared__ float w2;
	__shared__ float mag;

	if (iz==1){
		ss_bin = s_val[ix + iy*nx];
	   w2 = s_val[ix + iy*nx + nx*ny];
	   mag = s_val[ix + iy*nx + 2*nx*ny];
	}

	__syncthreads();
		

/*	cfloat ss_bin = s_val[ix + iy*nx];
	cfloat w2 = s_val[ix + iy*nx + nx*ny];
	cfloat mag = s_val[ix + iy*nx + 2*nx*ny];
*/	

// index of s is "is"
	cint is = floorf(ss_bin); // index of nearest neighbor in "s"

// Check if index out of bound
	if (is < 0 || is >= ns-1 || ix>=nx || iy>=ny ||(mask2[ix + iy*nx] != mask_id)) // each thread does its part only
		return;



	cfloat wr = ss_bin - is; // horizontal bilinear
	cfloat wl = 1. - wr; // interpolation factors
	//float *pi = image; //IMAGE!
	image += (ix + iy * nx) * nz + iz;	// NEW: don't need pi

	cfloat *pp1 = proj + is * nt;
	cfloat *pp2 = proj + (is+1) * nt;


		cfloat zz = dz * (iz - wz);
		cfloat tt = mag * zz;
		cfloat tt_bin = tt / dt + wt;
	// z value is used to determine index of t "it"
		cint it = floorf(tt_bin); // nearest nbr in "t"

		if (it < 0 || it >= nt-1) // out of FOV
			return;
		else { // reconstructing the image

			cfloat wu = tt_bin - it;
			cfloat wd = 1. - wu;
			cfloat p1 = wl * pp1[it]
				+ wr * pp2[it]; // interpolate
			cfloat p2 = wl * pp1[it+1]
				+ wr * pp2[it+1]; // horizontal

			// final vertical interpolation:
			*image += w2 * (wu * p1 + wd * p2);
		}
	
}
#endif
/* ******************************************************************
NEW: END OF SECOND KERNEL for z coordinate
******************************************************************** */











#ifdef fdk_gpu
static int iDivUp(int a, int b) {
	return (a % b != 0) ? (a / b + 1) : (a / b);
}
#endif

//
// fdk_ts_back1()
// The FDK backprojection is *added* to the image, so the user must zero it!
//
sof fdk_ts_back1_gpu(
float *image, // [nz nx ny] <- trick!
#ifdef fdk_gpu
float *s_val, // NEW: [nx ny] s values for each x,y pair
#endif
cint nx,
cint ny,
cint nz,
cfloat dx, // voxel size
cfloat dy, // can be negative to cause flip
cfloat dz,
cfloat offset_x, // image volume center offset in pixels (usually 0)
cfloat offset_y,
cfloat offset_z,
cbyte *mask2, // [nx ny] 2D support mask: 0, 1, ..., nthread
cbyte mask_id, // 1 ... nthread
cfloat dso, // distance from source to isocenter
cfloat dsd, // distance from source to detector
cfloat dfs, // distance from focal point to source (0 or inf)
cint ns, // projection view dimensions
cint nt,
cfloat ds, // horizontal ray spacing (view sample spacing)
cfloat dt, // vertical ray spacing (view sample spacing)
cfloat offset_s, // channel offset [pixels]
cfloat offset_t, // vertical offset on detector [pixels]
cfloat *proj, // [nt ns] <- trick! projection view at angle beta
cfloat beta) // source angle [radians]
{
	truf is_arc = 0;
	if (dfs == 0)
		is_arc = 1;
	else if (!Isinf(dfs))
		Warn("dfs not done - junk!")

	static truf told = False;
	if (!told) {
		Note2("nx=%d ny=%d", nx, ny)
		told = True;
	}

#ifdef fdk_gpu
//	jf_gpu_memset(image, 0, nx*ny*nz) // initialize device image to 0
#endif

#if 0 // the lines below crash. cpu should not address via device pointers!
	for (int ii=0; ii < nx*ny*nz; ++ii)
		image[ii] = 700.;
	Ok
#endif

#ifdef fdk_gpu
//	dim3 dimBlock(nx, ny);
//	dim3 dimBlock(nx/16, ny/16);
//	int tmp2 = 16;
//	int tmp = 16;
	dim3 dimBlock(16, 16);
//	dim3 dimBlock(32, 16);
//	dim3 dimBlock(1, 1);
	dim3 dimGrid(iDivUp(nx,dimBlock.x), iDivUp(ny,dimBlock.y));

// testing
//	printf("MY FILE COOLNESS!: yay GPU \n");
//	printf("nx:%d ny:%d nz%d ns%d nt%d \n", nx, ny, nz, ns, nt);

	fdk_ts_back1_kernel<<<dimGrid, dimBlock>>>(
// NEW
	s_val,
#else
// testing
// printf("MY FILE COOLNESS!: NOooo no GPU :( \n");

	for (int iy=0; iy < ny; ++iy)
	for (int ix=0; ix < nx; ++ix)
		fdk_ts_back1_kernel(ix, iy,
#endif
		image,
		nx,
		ny,
		nz,
		dx,
		dy,
		dz,
		offset_x,
		offset_y,
		offset_z,
		mask2,
		mask_id,
		dso,
		dsd,
		is_arc,
		ns,
		nt,
		ds,
		dt,
		offset_s,
		offset_t,
		proj,
		beta);

// NEW:
#ifdef fdk_gpu

// Testing only
//	int a = iDivUp(10,4);
//	printf("a: %d\n",a);

/*	int dBlock_x = 8;
	int dBlock_y = 8;
	int dBlock_z = 8;
	int nzLayers = iDivUp(nz, dBlock_z);	// number of z layers

	dim3 dimBlock2(dBlock_x, dBlock_y, dBlock_z);
	dim3 dimGrid2(iDivUp(nx*nzLayers,dimBlock2.x), iDivUp(ny,dimBlock2.y));
*/

// NEW 5: case where nz <=512
	dim3 dimBlock2(nz, 1, 1);		// all the z's in the same block for a given (x,y)
	int numBlock_x = ((nz/512)+1)*nx;
   dim3 dimGrid2(numBlock_x, ny);


	// Testing:
//	printf("dimGrid2: %d %d\n",dimGrid2.x, dimGrid2.y);

	fdk_ts_back1_kernel_2<<<dimGrid2, dimBlock2>>>(
		s_val, // [nx ny 4] s values, w2, mag for x,y pairs (2D)
		image, // [nz nx ny] <- trick!
		nx,
		ny,
		nz,
		//float dx, // voxel size
		//float dy, // can be negative to cause flip
		dz,
		//float offset_x, // image volume center offset in pixels (usually 0)
		//float offset_y,
		offset_z,
		mask2, // [nx ny] 2D support mask: 0, 1, ..., nthread
		mask_id, // 1 ... nthread
		//float dso, // distance from source to isocenter
		//float dsd, // distance from source to detector
		//truf is_arc,
		ns, // projection view dimensions
		nt,
		//float ds, // horizontal ray spacing (view sample spacing)
		dt, // vertical ray spacing (view sample spacing)
		//float offset_s, // channel offset [pixels]
		offset_t, // vertical offset on detector [pixels]
		proj);

#endif

	Ok
}
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
// fdk-ts-h.cu
#include <stdio.h>

void fdk_ts_help(void)
{
	printf("\n\
\n\
	image = function('fdk,ts,back', nx,ny,nz, dx,dy,dz, \n\
		offset_x, offset_y, offset_z, mask2, \n\
		dso, dsd, ds, dt, offset_s, offset_t, proj, beta, nthread)\n\
\n\
		image output is single [nz nx ny] <- trick!\n\
		nx,ny,nz: (int32) image size\n\
		dx,dy,dz: (double) voxel size\n\
		offset_x,_y,_z: (double) center offset in pixels (usually 0)\n\
		mask2: (uint8) [nx ny] 2D support mask\n\
		dso: (double) distance from source to isocenter\n\
		dsd: (double) distance from source to detector\n\
		dfs: (double) distance from focal point to source (0 or inf)\n\
		ds: (double) horizontal ray spacing\n\
		dt: (double) vertical ray spacing\n\
		offset_s: (double) channel offset [pixels]\n\
		offset_t: (double) vertical offset on detector [pixels]\n\
		nthread: (int32) # of processors\n\
		proj: (single) [nt ns na] (trick!) projection view for each beta\n\
		beta: (double) [na] source angle(s) [radians]\n\
	(CUDA version)\n\
\n");
}
// cbct,mask2.c
// mask operations

#include "cbct,def.h"

//
// cbct_mask_init()
// initialize mask with values for each thread in it
// input and output masks can be the same pointer!
//
sof cbct_mask_init(
byte *mask_int,	// [nx ny] out: integers 0 or 1 ... nthread
cbyte *mask_bin, // [nx ny] in: binary, possibly null
cint nx, cint ny,
cint nthread,
ctruf chat)
{
	if (nthread > 250) Fail("too many threads")

	long sum = 0;
	if (!mask_bin)
		sum = nx * ny;
	else
		for (int jj=0; jj < nx * ny; ++jj)
			sum += mask_bin[jj] != 0;

	const int nj_per_thread = (int) Ceilf(sum / (float) nthread);
	if (chat > 99)
		Note3("%ld voxels / %d threads = %d each",
			sum, nthread, nj_per_thread)

	int jj=0;
	for (int it=1; it <= nthread; ++it) {
		int nj_here = Min(sum, it * nj_per_thread) - (it-1) * nj_per_thread;
		if (chat > 99) Note2("thread %d nj %d", it-1, nj_here)

		while (nj_here) {
			if (!mask_bin || mask_bin[jj]) {
				mask_int[jj] = it;
				--nj_here;
			}
			else
				mask_int[jj] = 0;
			++jj;
		}
	}
	// if (jj != nx * ny) Fail("bug")
	// Iwrite2byte("mask-int.fld", mask_int, nx, ny)
	Ok
}
// fdk,st.c
// Feldkamp aka FDK back-projection for arc/flat detector.
// For detector index (s,t).
// Copyright 2005-6-27, Jeff Fessler, University of Michigan

#include "defs-env.h"
#include "def,fdk.h"

void fdk_st_help(void)
{
	printf("\n\
\n\
	image = function('fdk,st,{arc|flat},back', nx,ny,nz, dx,dy,dz, \n\
		offset_x, offset_y, offset_z, mask2, \n\
		dso, dsd, ds, dt, offset_s, offset_t, proj, beta, nthread)\n\
\n\
		image output is single [nx ny nz]\n\
		nx,ny,nz: (int32) image size\n\
		dx,dy,dz: (double) voxel size\n\
		offset_x,_y,_z: (double) center offset in pixels (usually 0)\n\
		mask2: (uint8) [nx,ny] 2D support mask\n\
		dso: (double) distance from source to isocenter\n\
		dsd: (double) distance from source to detector\n\
		ds: (double) horizontal ray spacing\n\
		dt: (double) vertical ray spacing\n\
		offset_s: (double) channel offset [pixels]\n\
		offset_t: (double) vertical offset on detector [pixels]\n\
		proj: (single) [ns nt na] projection view at angle beta\n\
		beta: (double) [na] source angle(s) [radians]\n\
\n");
}


//
// fdk_st_back1()
// The FDK backprojection is *added* to the image, so the user must zero it!
//
void fdk_st_back1(
float *image, // [nx ny nz]
cint nx,
cint ny,
cint nz,
cfloat dx,
cfloat dy, // can be negative to cause flip
cfloat dz,
cfloat offset_x, // center offset in pixels (usually 0)
cfloat offset_y,
cfloat offset_z,
cbyte *mask2, // [nx ny] 2D support mask
#if 0
cbyte mask_id, // 1 + thread id
#endif
cfloat dso, // distance from source to isocenter
cfloat dsd, // distance from source to detector
cint ns,
cint nt,
cfloat ds, // horizontal ray spacing
cfloat dt, // vertical ray spacing
cfloat offset_s, // channel offset [pixels]
cfloat offset_t, // vertical offset on detector [pixels]
cfloat *proj, // [ns nt] projection view at angle beta
cfloat beta, // source angle [radians]
ctruf is_arc)
{
	cint nxy = nx * ny;
	cfloat wx = (nx-1)/2. + offset_x;
	cfloat wy = (ny-1)/2. + offset_y;
	cfloat wz = (nz-1)/2. + offset_z;
	cfloat ws = (ns-1)/2. + offset_s;
	cfloat wt = (nt-1)/2. + offset_t;
	cfloat sinb = sin(beta);
	cfloat cosb = cos(beta);

#if 0
	printf("nxyz=%d,%d,%d dxyz=%g,%g,%g offset_xyz=%g,%g,%g\n"
		"beta=%g dso=%g dsd=%g ds=%g dt=%g offset_st=%g,%g\n",
		nx, ny, nz, dx, dy, dz,
		offset_x, offset_y, offset_z,
		beta, dso, dsd, ds, dt, offset_s, offset_t);
#endif

	// loop over pixels
	for (int iy = 0; iy < ny; ++iy) {
		cfloat yy = dy * (iy - wy);
	 for (int ix = 0; ix < nx; ++ix, ++image) {
		// if (mask2[ix + iy*nx] != mask_id)
		if (!mask2[ix + iy*nx])
			continue;

		cfloat xx = dx * (ix - wx);
		cfloat xbeta = xx * cosb + yy * sinb;
		cfloat ybetas = dso - (-xx * sinb + yy * cosb);
		cfloat mag = dsd / ybetas;
		cfloat ss = is_arc ? (dsd * atan2(xbeta, ybetas))
				: (mag * xbeta);
		cfloat ss_bin = ss / ds + ws;
		cint is = (int) Floorf(ss_bin);

		if (is < 0 || is >= ns-1) // out of FOV
			continue;

		cfloat wr = ss_bin - is;
		cfloat wl = 1. - wr;
		cfloat w2 = is_arc ?
			(Sqr(dsd) / (Sqr(ybetas) + Sqr(xbeta))) : Sqr(mag);

		for (int iz = 0; iz < nz; ++iz) {
			cfloat zz = dz * (iz - wz);
			cfloat tt = mag * zz;
			cfloat tt_bin = tt / dt + wt;
			cint it = (int) Floorf(tt_bin);

			if (it < 0 || it >= nt-1) // out of FOV
				continue;

			cfloat wu = tt_bin - it;
			cfloat wd = 1. - wu;
			cfloat p1 = wl * proj[is + it*ns]
					+ wr * proj[is+1 + it*ns];
			cfloat p2 = wl * proj[is + (it+1)*ns]
					+ wr * proj[is+1 + (it+1)*ns];

			image[iz*nxy] += w2 * (wu * p1 + wd * p2);
		}
	 }
	}
}
/*
* mexarg.c
* Show matlab mex arguments (for debugging and error messages)
*
* Copyright 01-04-23, Jeff Fessler, University of Michigan
*/
#include "def,mexarg.h"

#ifdef Mmex // needed for wt.c -> wtfmex.c

/*
* mxu_string()
* caller must free using mxu_string_free()
*/
char *mxu_string(Const mxArray *mx, cchar *arg)
{
	char	*string;
	int	n = mxGetM(mx) * mxGetN(mx) + 1;

	if (!mxIsChar(mx))
		Fail1("%s must be char array", arg)

#if 1
	Call(string = (char *) mxCalloc, (n, sizeof(char)))
#else
	Mem0(string, n)
#endif
	if (mxGetString(mx, string, n))
		Warn("bug with mxGetString")
	return string;
}

sof mxu_string_free(char *s)
{
#if 1
	mxFree(s);
#else
	Free0(s)
#endif
	Ok
}

/*
* mxu_showdim()
*/
static sof mxu_showdim(Const mxArray *mx)
{
	int id, ndim;
	Const mwSize *dims;

	Call(ndim = mxGetNumberOfDimensions, (mx))
	Call(dims = mxGetDimensions, (mx))

	printf("dims");
	for (id=0; id < ndim; ++id)
		printf(" %d", (int) dims[id]);
#if 0
	if (ndim == 1)	Note1("dims %d", dims[0])
	if (ndim == 2)	Note2("dims %d %d", dims[0], dims[1])
	if (ndim == 3)	Note3("dims %d %d %d", dims[0], dims[1], dims[2])
	if (ndim == 4)	Note4("dims %d %d %d %d", dims[0], dims[1], dims[2], dims[3])
	if (ndim > 4)	Note4("dims %d %d %d ... %d", dims[0], dims[1], dims[2], dims[ndim-1])
#endif
	Ok
}


/*
* mxu_arg()
* Show arguments
*/
sof mxu_arg(cint nmx, Const mxArray *pmx[])
{
	int ii;

	Note1("narg=%d", nmx)

	for (ii=0; ii < nmx; ++ii) {
		Const mxArray *mx = pmx[ii];

		printf("arg %d, ", ii);
		Call(mxu_showdim, (mx))
		printf(", ");

		if (mxIsChar(mx)) {
			char *arg;
			Call(arg = mxu_string, (mx, ""))
			printf("char, '%s'", arg);
			Call(mxu_string_free, (arg))
		}

		else if (mxIsUint8(mx)) {
			byte val = *((cbyte *) mxGetData(mx));
			printf("uint8, val[0] = %d", (int) val);
		}

		else if (mxIsInt32(mx)) {
			int val = mxGetInt(mx);
			printf("int32, val[0] = %d", val);
		}

		else if (mxIsSingle(mx)) {
			float val = *((cfloat *) mxGetData(mx));
			printf("single, val[0] = %g", val);
		}

		else if (mxIsDouble(mx)) {
			double val = *((cdouble *) mxGetData(mx));
			printf("double, val[0] = %g", val);
		}

		else if (mxIsCell(mx))
			printf("cell");

		else if (mxIsStruct(mx))
			printf("struct");

		else
			printf("UNKNOWN!?");

		printf("\n");
	}

	Ok
}


/*
* mxu_numel()
* # of elements in a matlab array
*/
int mxu_numel(Cmx mx)
{
	int ndim;
	Const mwSize *dims;

	Call(ndim = mxGetNumberOfDimensions, (mx))
	Call(dims = mxGetDimensions, (mx))

	int numel = 1;
	int ii; // for cuda
	for (ii=0; ii < ndim; ++ii)
		numel *= dims[ii];

	return numel;
}


#endif // Mmex
// jf,thread1.c
// generic routines for multi-threading
// for linux, now using "Portable Linux Processor Affinity (PLPA)"
// http://www.open-mpi.org/software/plpa/overview.php
//
// Copyright 2006-5-3, Jeff Fessler, University of Michigan


#define jf_nthread_big 32 // avoid alloc/free if nthread <= this value

#ifdef Use_thread

# ifdef Use_nptl // linux: native posix thread library
#  ifndef _GNU_SOURCE
#   define _GNU_SOURCE
#  endif
#  include <sched.h>
#  include <nptl/pthread.h>

# elif Use_plpa
#  include <plpa.h>
#  include <pthread.h>

# else // mac os
#  include <pthread.h>
# endif

#endif // Use_thread

#include "jf,thread1.h"

#ifdef Use_ncore_sysctl // old way for mac
#include <sys/sysctl.h>
#endif

#include <unistd.h> // for sysconf() to get ncore

#ifdef Use_aff_mac1 // affinity for mac leopard (and above?)
#include <mach/thread_policy.h>
#include <mach/mach.h>
#endif

/*
* jf_thread1_ncore()
* return number of available cores, if possible
* caution: returns 0 on failure
*/
int jf_thread1_ncore(cint nwant)
{
	int ncore = 0;

#ifdef Use_ncore_sysctl // old way for mac
        int mib[2] = {CTL_HW, HW_NCPU};
        size_t len;
        len = sizeof(ncore);
        sysctl(mib, 2, &ncore, &len, NULL, 0);
	if (ncore <= 0)
		Fail1("sysctl returned %d", ncore)
#endif

#if 1 // this works on mac osx and linux
	ncore = sysconf( _SC_NPROCESSORS_ONLN );
	if (ncore == -1)
		Fail("sysconf() failed")
#endif

	if (nwant == -1) { // user wants as many as possible
		if (ncore)
			return ncore;
		else
			Fail("cannot determine # of cores")
	}

	else if (nwant == 0) { // user would like many, but will accept 1
		if (ncore)
			return ncore;
		else {
			Warn("cannot determine # of cores, defaulting to 1")
			return 1;
		}
	}

	else {
		if (ncore && ncore < nwant)
			Fail2("want %d cores but have only %d", nwant, ncore)
	}

        return ncore;
}


/*
* fake temporary routine for setting affinity.
* needed only if gcc compiler cannot find the real one.
* probably superceded by plpa library
*/
#ifdef Provide_setaffinity
void pthread_attr_setaffinity_np(void)
{
	static int warned = 0;
	if (!warned) {
		Note("calling dummy setaffinity")
		warned = 1;
	}
}
#endif


#ifdef Use_aff_mac1

/*
* jf_thread1_setaffinity_mac()
* for mac (leopard and above?) we set affinity after starting the thread
* but before doing any work!?
* apple thread affinity help, but does not mention pthread:
* http://developer.apple.com/releasenotes/Performance/RN-AffinityAPI/
* for example see:
* http://www.opensource.apple.com/darwinsource/projects/other/xnu-1228.3.13/tools/tests/affinity/sets.c
*/
static sof jf_thread1_setaffinity_mac(cint ithread)
// Const jf_thread1_affinity *aff, // affinity control
{
	thread_extended_policy_data_t epolicy;
	epolicy.timeshare = FALSE;
	kern_return_t ret = thread_policy_set(
		mach_thread_self(), THREAD_EXTENDED_POLICY,
		(thread_policy_t) &epolicy, THREAD_EXTENDED_POLICY_COUNT);
	if (ret != KERN_SUCCESS)
		Fail1("thread_policy_set returned %d", ret)

	thread_affinity_policy_data_t apolicy;
	apolicy.affinity_tag = ithread + 1; // set affinity tag

	ret = thread_policy_set(
		mach_thread_self(), THREAD_EXTENDED_POLICY,
		(thread_policy_t) &apolicy, THREAD_EXTENDED_POLICY_COUNT);
	if (ret != KERN_SUCCESS)
		Fail1("thread_policy_set returned %d", ret)

	Ok
}

#endif // Use_aff_mac1


/*
* jf_thread1_affinity_check()
*/
sof jf_thread1_affinity_check(cint chat)
{
#if Use_nptl
	if (chat)
		Note("using nptl, so affinity should work")
#elif Use_plpa
	Call(PLPA_PROBE_OK == plpa_api_probe, ())
	if (chat)
		Note("using plpa, and affinity probe ok")
#elif Use_aff_mac1
	if (chat)
		Note("using mac affinity sets, which i hope works")
#else
	if (chat)
		Warn("affinity check called without support, disregarding")
#endif
	Ok
}


/*
* jf_thread1_setaffinity_attr()
* match thread to a given cpu
*/
static sof jf_thread1_setaffinity_attr(
pthread_attr_t *attr,
cint ithread,
Const jf_thread1_affinity *aff, // affinity control
cint chat)
{
	if (!aff || aff->type == jf_thread1_affinity_none) Ok

#if Use_nptl || Use_plpa
	int affinity = ithread; // usual affinity
	if (aff->type == jf_thread1_affinity_mod && aff->nmod)
		affinity = ithread % aff->nmod;
	if (aff->type == jf_thread1_affinity_list && aff->list)
		affinity = aff->list[ithread];
#endif

#if Use_nptl
	{
	cpu_set_t cs;
	size_t cpu_set_size = sizeof(cs);
	__CPU_ZERO(&cs);
	__CPU_SET(affinity, &cs);
	pthread_attr_setaffinity_np(attr, cpu_set_size, &cs);
	if (chat) Note2("set affinity for thread %d to %d", ithread, affinity)
	}

#elif Use_plpa
	(void) attr;
	(void) chat;
	{
	plpa_cpu_set_t cs;
	size_t cpu_set_size = sizeof(cs);
	int ret;
	PLPA_CPU_ZERO(&cs);
	PLPA_CPU_SET(affinity, &cs);
	ret = plpa_sched_setaffinity(0, cpu_set_size, &cs);
	if (ret)
		Fail2("plpa_sched_setaffinity(affinity=%d) returned %d\n"
			"Perhaps you tried to use more threads than cores??",
			affinity, ret)
	}

#elif Use_aff_mac1
	// mac doesn't use attr to set affinity
	(void) attr;
	(void) chat;
	if (ithread == 0 && aff->type != jf_thread1_affinity_try)
		Fail("mac version supports only basic affinity support")

#else
	(void) attr;
	(void) chat;
	if (ithread == 0 && aff->type != jf_thread1_affinity_try)
		Warn("affinity support requested but not enabled!")
#endif
	Ok
}


/*
* jf_thread1_glue()
* interface routine for threads
*/
static void *jf_thread1_glue(void *in)
{
	jf_thread1_s *pt = (jf_thread1_s *) in;

// 2009-6-2 found that PRTS increases every time we call this!?
#ifdef Use_aff_mac1
	if (!jf_thread1_setaffinity_mac(pt->id)) {
		pt->ok = sof_failure;
		return NULL;
	}
#endif // Use_aff_mac1
//	(void) jf_thread1_setaffinity_mac;

	pt->ok = (pt->init)(pt->ps, pt->id, pt->nthread);

//	pthread_exit((void*) in); // 2009-6-2 per llnl example
	return NULL; // ""
}


/*
* jf_thread1_tops()
* top-level interface to threaded operations
* trick: only one of "ps" or "pps" should be used!
*/
sof jf_thread1_tops(
jf_thread1_init_t fun_init, // required user function
jf_thread1_wrap_t fun_wrap, // optional user function
void *ps, // pointer to data structure used by threads
void **pps, // [nthread] pointers to structures ""
cint nthread, // # threads
Const jf_thread1_affinity *aff, // affinity control
cint chat)
{
	int it;
	jf_thread1_s *pt;
	jf_thread1_s pt_pre[jf_nthread_big];

	if (nthread > jf_nthread_big) {
		Warn1("allocating space for %d threads", nthread)
		Mem0pure(pt, nthread)
	}
	else
		pt = pt_pre;

	if (ps && pps) Fail("only one of 'ps' and 'pps' may be non-null")
	if (!ps && !pps) Fail("one of 'ps' and 'pps' must be non-null")

	for (it=0; it < nthread; ++it) {
		pt[it].init = fun_init;
		pt[it].ok = sof_failure;
		pt[it].id = it;
		pt[it].nthread = nthread;
		if (ps)
			pt[it].ps = ps; // all threads get same structure!
		else
			pt[it].ps = pps[it]; // each thread gets its own
	}

	if (nthread == 1) { // to support non-threaded compiles
		jf_thread1_glue(pt+0);
		if (!pt[0].ok)
			Fail("single thread failed")

		if (fun_wrap)
			Warn("fun_wrap unused for nthread=1")
	}

#ifdef Use_thread
	else {
		pthread_t *pid;
		pthread_t pid_pre[jf_nthread_big];

		if (nthread > jf_nthread_big)
			Mem0pure(pid, nthread)
		else
			pid = pid_pre;

		pthread_attr_t attr_, *p_attr = &attr_;
		pthread_attr_init(p_attr);
		pthread_attr_setdetachstate(p_attr, PTHREAD_CREATE_JOINABLE);

		for (it=0; it < nthread; ++it) {
			// match thread to a given cpu, if requested
			Call(jf_thread1_setaffinity_attr, (p_attr, it, aff, chat))

			if (pthread_create(pid+it, p_attr, jf_thread1_glue,
				(void *) (pt+it)))
				Fail1("error creating thread %d", it)
		}

		if (pthread_attr_destroy(p_attr))
			Fail("pthread_attr_destroy()")

		for (it=0; it < nthread; ++it) {
			if (pthread_join(pid[it], NULL))
				Fail1("pthread_join %d failed", it)
			if (!pt[it].ok)
				Fail1("thread %d failed", it)
		}

		if (fun_wrap)
			Call(fun_wrap, (pt, nthread))

		if (nthread > jf_nthread_big)
			Free0pure(pid)
	}
#else

	else
		Fail1("threads %d not done", nthread)

	(void) fun_wrap;
	(void) aff;
	(void) chat;
	(void) jf_thread1_setaffinity_attr;

#endif

	if (nthread > jf_nthread_big)
		Free0pure(pt)

	Ok
}


/*
* jf_thread1_top()
* simpler top-level interface to threaded operations
*/
sof jf_thread1_top(
jf_thread1_init_t fun_init, // required user function
jf_thread1_wrap_t fun_wrap, // optional user function
void *ps, // pointer to data structure passed to threads
cint nthread, // # threads
cint chat)
{
	jf_thread1_affinity aff_, *aff = &aff_;
	aff->type = jf_thread1_affinity_try;
	Call(jf_thread1_tops,
		(fun_init, fun_wrap, ps, NULL, nthread, aff, Chat))
	Ok
}
