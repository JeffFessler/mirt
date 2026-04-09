// fdk-ts.cu
// Feldkamp aka FDK backprojection for arc/flat detector.
// For detector index (t,s).
// CUDA/GPU version
// Copyright 2005-6-27, Jeff Fessler, University of Michigan

#include "jf-cuda.h"
#include "def,fdk.h"
#include "fdk-gpu.h"




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
cint ix,
cint iy,
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
	cint ix = blockIdx.x * blockDim.x + threadIdx.x;
	cint iy = blockIdx.y * blockDim.y + threadIdx.y;
	
#endif

	if (ix >= nx || iy >= ny)
		return;
	
	image += (ix + iy * nx) * nz;;

	if (mask2[ix + iy*nx] != mask_id) // each thread does its part only
		return;

#ifdef fdk_gpu
	__shared__ float shared_img[2][2][120];
	float *temp_img = shared_img[threadIdx.x][threadIdx.y];
	
	//extern __shared__ float shared_img[];
	//float *temp_img = shared_img + (threadIdx.y + (threadIdx.x)*5)*nz;
	
	
	#pragma unroll
	for(int i=0; i < nz; i++){
		temp_img[i] = image[i];
	}
	
#endif


	cfloat wx = (nx-1)/2. + offset_x;
	cfloat wy = (ny-1)/2. + offset_y;
	cfloat wz = (nz-1)/2. + offset_z;
	cfloat ws = (ns-1)/2. + offset_s;
	cfloat wt = (nt-1)/2. + offset_t;
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
	cint is = floorf(ss_bin); // index of nearest neighbor in "s"

	if (is < 0 || is >= ns-1) // out of FOV
		return;

	cfloat w2 = is_arc ? // fan-beam image domain weighting
		(Sqr(dsd) / (Sqr(ybetas) + Sqr(xbeta))) : Sqr(mag);

	cfloat wr = ss_bin - is; // horizontal bilinear
	cfloat wl = 1. - wr; // interpolation factors
	
#ifdef fdk_gpu
	float *pi = temp_img;
#else
	float *pi = image;
#endif
	
	float *pt = image;
	cfloat *pp1 = proj + is * nt;
	cfloat *pp2 = proj + (is+1) * nt;
	
	#pragma unroll 120
	for (int iz = 0; iz < nz; ++iz, ++pi, ++pt) { // slice loop
		cfloat zz = dz * (iz - wz);
		cfloat tt = mag * zz;
		cfloat tt_bin = tt / dt + wt;
		cint it = floorf(tt_bin); // nearest nbr in "t"

		if (it < 0 || it >= nt-1) // out of FOV
			continue;
		else {
			cfloat wu = tt_bin - it;
			cfloat wd = 1. - wu;
			cfloat p1 = wl * pp1[it]
				+ wr * pp2[it]; // interpolate
			cfloat p2 = wl * pp1[it+1]
				+ wr * pp2[it+1]; // horizontal

			// final vertical interpolation:
			//*pi += w2 * (wu * p1 + wd * p2);
			*pt = *pi + w2 * (wu * p1 + wd * p2);
			//*pt += + w2 * (wu * p1 + wd * p2);
		}
	}

/*
#ifdef fdk_gpu
	__syncthreads();
	#pragma unroll
	for(int i=0; i<nz; i++){
		image[i] = temp_img[i];
	}
#endif
*/
	
}


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
//	dim3 dimBlock(nx, ny);
	dim3 dimBlock(2, 2);
	dim3 dimGrid(iDivUp(nx,dimBlock.x), iDivUp(ny,dimBlock.y));
	
	//call when using internal shared memory
	fdk_ts_back1_kernel<<<dimGrid, dimBlock>>>(
	
	//call when using extern shared memory
	//fdk_ts_back1_kernel<<<dimGrid, dimBlock, 6*5*nz*sizeof(float)>>>(
#else
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
	Ok
}








