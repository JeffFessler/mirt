// Version 3.0
// proj in shared memory
// failed


// fdk-ts.cu
// Feldkamp aka FDK backprojection for arc/flat detector.
// For detector index (t,s).
// CUDA/GPU version
// Copyright 2005-6-27, Jeff Fessler, University of Michigan

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
float *s_val, // [nx ny 3] s values, w2, mag for x,y pairs (2D)
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
float *s_val, // [nx ny 4] s values, w2, mag for x,y pairs (2D)
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

//	cint ix = fmodf((blockIdx.x*blockDim.x + threadIdx.x), nx);
//	cint iy = blockIdx.y*blockDim.y + threadIdx.y;
//	cint iz = floorf((blockIdx.x*blockDim.x+threadIdx.x)/nx)*blockDim.z + threadIdx.z;
	cint ix = blockIdx.x*2 + threadIdx.y;
	cint iy = blockIdx.y;
	cint iz = threadIdx.x;

	cfloat ss_bin = s_val[ix + iy*nx]; //NEW2: position important
	cfloat wz = (nz-1)/2. + offset_z;
	cfloat wt = (nt-1)/2. + offset_t;

	// NEW3:
//	extern __shared__ float pp1_s[blockDim.x*blockDim.y][blockDim.z];
//	extern __shared__ float pp2_s[blockDim.x*blockDim.y][blockDim.z];
//	extern __shared__ float pp1_s[];
//	extern __shared__ float pp2_s[];
//	__shared__ float pp1_s[512];
//	__shared__ float pp2_s[512]
	__shared__ float proj_s[1024];
	

// index of s is "is"
	cint is = floorf(ss_bin); // index of nearest neighbor in "s"

// Check if index out of bound
	if (is < 0 || is >= ns-1 || ix>=nx || iy>=ny ||(mask2[ix + iy*nx] != mask_id)) // each thread does its part only
		return;


	cfloat mag = s_val[ix + iy*nx + 2*nx*ny];	// NEW2: position


	cfloat wr = ss_bin - is; // horizontal bilinear
	cfloat wl = 1. - wr; // interpolation factors
	//float *pi = image; //IMAGE!
	image += (ix + iy * nx) * nz + iz;	// NEW: don't need pi

// NEW3:
//	cfloat *pp1 = proj + is * nt;
//	cfloat *pp2 = proj + (is+1) * nt;


	cfloat zz = dz * (iz - wz);
	cfloat tt = mag * zz;
	cfloat tt_bin = tt / dt + wt;
	// z value is used to determine index of t "it"
	cint it = floorf(tt_bin); // nearest nbr in "t"

	if (it < 0 || it >= nt-1) // out of FOV
		return;

	//load into shared memory
//	pp1_s[(blockDim.x*threadIdx.y+threadIdx.x)*blockDim.z+threadIdx.z] = *(proj+is*nt+it);
//	pp2_s[(blockDim.x*threadIdx.y+threadIdx.x)*blockDim.z+threadIdx.z] = *(proj+(is+1)*nt+it);

//	cint index = (fmodf(ix,4)+fmodf(iy,2)*4)*64+fmodf(iz,64);
//	pp1_s[index] = *(proj+is*nt+it);
//	pp2_s[index] = *(proj+(is+1)*nt+it);
	proj_s[threadIdx.y*2*nz+iz]= *(proj+is*nt+iz);		//2*nz>=nt
	proj_s[threadIdx.y*2*nz+nz+iz]= *(proj+is*nt+nz+iz);

	__syncthreads();


 // reconstructing the image
	cfloat w2 = s_val[ix + iy*nx + nx*ny];	// NEW 2: position			


	cfloat wu = tt_bin - it;
	cfloat wd = 1. - wu;
/*	NEW3:
	cfloat p1 = wl * pp1[it]
		+ wr * pp2[it]; // interpolate
	cfloat p2 = wl * pp1[it+1]
		+ wr * pp2[it+1]; // horizontal */
	cint neg = (~(threadIdx.y&1))&1;
	cfloat p1 = wl * proj_s[threadIdx.y*2*nz + it] + wr * proj_s[neg*2*nz+it]; // interpolate

	cfloat p2 = wl * proj_s[threadIdx.y*2*nz + it+1]+ wr * proj_s[neg*2*nz+it+1]; // horizontal */

	// final vertical interpolation:
	*image += w2 * (wu * p1 + wd * p2);
	
	__syncthreads();	
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

	//ATTN: switch x and z
	int dBlock_x = nz; // NEW: nz <= 512
	int dBlock_y = 2;
	int dBlock_z = 1;
	//int nzLayers = iDivUp(nz, dBlock_z);	// number of z layers

	dim3 dimBlock2(dBlock_x, dBlock_y, dBlock_z);
//	dim3 dimGrid2(iDivUp(nx*nzLayers,dimBlock2.x), iDivUp(ny,dimBlock2.y));
	dim3 dimGrid2(nx/2,ny);

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
