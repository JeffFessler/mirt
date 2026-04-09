// fdk-ts.cu
// Feldkamp aka FDK backprojection for arc/flat detector.
// For detector index (t,s).
// CUDA/GPU version
// Copyright 2005-6-27, Jeff Fessler, University of Michigan

#include "jf-cuda.h"
#include "def,fdk.h"
#include "fdk-gpu.h"

#define kernel_loop



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
cfloat *proj, // [nt ns] <- trick! projection view at angle beta
cbyte *mask2, // [nx ny] 2D support mask: 0, 1, ..., nthread
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
float beta) // source angle [radians]
{
#ifdef fdk_gpu
	// index into image array
	cint ix = blockIdx.x * blockDim.x + threadIdx.x;
	cint iy = blockIdx.y * blockDim.y + threadIdx.y;
#endif

	if (ix >= nx || iy >= ny)
		return;
	
	int img_index = (ix + iy * nx) * nz;		
	image += img_index;

#ifdef fdk_gpu
	if (tex1Dfetch(tex_mask2, ix + iy*nx) != mask_id) // each thread does its part only
		return;
#else
	if (mask2[ix + iy*nx] != mask_id) // each thread does its part only
		return;
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
	
	float *pi = image;
	int pp1_index = is * nt;
	int pp2_index = (is+1) * nt;

#ifndef fdk_gpu
	cfloat *pp1 = proj + pp1_index;
	cfloat *pp2 = proj + pp2_index;
#endif
	
	int iz_min, iz_max;
#ifdef kernel_loop
	if (dz / dt < 0) {
		iz_min = Ceilf(wz - wt * abs(dt) / (mag * dz));
		iz_max = Ceilf(wz + (nt-1 - wt) * abs(dt) / (mag * dz));
	}
	else{
		iz_min = Ceilf(wz - wt * dt / (mag * dz));
		iz_max = Ceilf(wz + (nt-1 - wt) * dt / (mag * dz));
	}
	
	iz_min = Max(iz_min, 0);
	iz_max = Min(iz_max, nz);

#else
	iz_min = 0;
	iz_max = nz;
#endif

	pi += iz_min;
	#pragma unroll
	for (int iz = iz_min; iz < iz_max; ++iz, ++pi) { // slice loop
		cfloat zz = dz * (iz - wz);
		cfloat tt = mag * zz;
		cfloat tt_bin = tt / dt + wt;
		cint it = floorf(tt_bin); // nearest nbr in "t"

		cfloat wu = tt_bin - it;
		cfloat wd = 1. - wu;
		
		#ifndef kernel_loop
		if (it < 0 || it >= nt-1) // out of FOV
			continue;
		#endif

		#ifdef fdk_gpu		
			
			#ifdef tex_1d
			cfloat p1 = wl * tex1Dfetch(tex_proj, pp1_index + it)
				+ wr * tex1Dfetch(tex_proj, pp2_index + it); // interpolate		
			
			cfloat p2 = wl * tex1Dfetch(tex_proj, pp1_index + it+1)
				+ wr * tex1Dfetch(tex_proj, pp2_index + it+1); // horizontal
			// final vertical interpolation:			
			*pi = tex1Dfetch(tex_img, img_index + iz) + w2 * (wu * p1 + wd * p2);			
			#else
			
			*pi = tex1Dfetch(tex_img, img_index + iz) + tex2D(tex_proj2d, tt_bin, ss_bin) * w2;			
			
			#endif

		#else	//CPU version
			cfloat p1 = wl * pp1[it]
				+ wr * pp2[it]; // interpolate
			cfloat p2 = wl * pp1[it+1]
				+ wr * pp2[it+1]; // horizontal
			// final vertical interpolation:
			*pi += w2 * (wu * p1 + wd * p2);
		#endif		
	}
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
#else
	for (int iy=0; iy < ny; ++iy)
	for (int ix=0; ix < nx; ++ix)
		fdk_ts_back1_kernel(ix, iy, proj, mask2,
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
		beta);
	Ok
}








