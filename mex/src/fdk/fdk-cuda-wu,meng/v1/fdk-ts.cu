// fdk-ts.cu
// Feldkamp aka FDK backprojection for arc/flat detector.
// For detector index (t,s).
// CUDA/GPU version
// Copyright 2005-6-27, Jeff Fessler, University of Michigan

// VERSION 5: shared memory 

#include "jf-cuda.h"
#include "def,fdk.h"
#include "fdk-gpu.h"
#include "omp.h"

#ifndef fdk_gpu

sof fdk_ts_back(
float *image, // [nz nx ny] <- trick!
cint nx,
cint ny,
cint nz,
cfloat dx, // voxel size
cfloat dy, // can be negative to cause flip
cfloat dz,
cfloat wx, // image volume center offset in pixels (usually 0)
cfloat wy,
cfloat wz,
//cbyte *mask2, // [nx ny] 2D support mask: 0, 1, ..., nthread
//cbyte mask_id, // 1 ... nthread
cfloat dso, // distance from source to isocenter
cfloat dsd, // distance from source to detector
cfloat dfs, // distance from focal point to source (0 or inf)
cint ns, // projection view dimensions
cint nt,
cfloat ds, // horizontal ray spacing (view sample spacing)
cfloat dt, // vertical ray spacing (view sample spacing)
cfloat ws, // channel offset [pixels]
cfloat wt, // vertical offset on detector [pixels]
cfloat *proj, // [nt ns] <- trick! projection view at angle beta
cfloat beta,
truf is_arc) // source angle [radians]
{
	float sinb = sinf(beta);
	float cosb = cosf(beta);
	image += (ix + iy * nx) * nz;
	for (int iy=0; iy < ny; ++iy){
		for (int ix=0; ix < nx; ++ix){
			if (ix >= nx || iy >= ny)
				return 0;
//			if (mask2[ix + iy*nx] != mask_id) 
				// each thread does its part only
//				return 0;
			cfloat yy = dy * iy - wy;
			cfloat xx = dx * ix - wx;
			cfloat xbeta = xx * cosb + yy * sinb;
			cfloat ybetas = dso - (-xx * sinb + yy * cosb);
			cfloat mag = dsd / ybetas;
			cfloat ss = is_arc ? (dsd * atan2f(xbeta, ybetas)): (mag * xbeta);
			cfloat ss_bin = ss / ds + ws;
			// index of s is "is"
			cint is = floorf(ss_bin); // index of nearest neighbor in "s"

			// Check if index out of bound
			if (is < 0 || is >= ns-1) // out of FOV
				return 0;
			cfloat w2 = is_arc ? // fan-beam image domain weighting
				(Sqr(dsd) / (Sqr(ybetas) + Sqr(xbeta))) : Sqr(mag);

			cfloat wr = ss_bin - is; // horizontal bilinear
			cfloat wl = 1. - wr; // interpolation factors
			float *pi = image; //IMAGE!
			cfloat *pp1 = proj + is * nt;
			cfloat *pp2 = proj + (is+1) * nt;

			// with ix, iy fixed, loop through all possible iz values
			for (int iz = 0; iz < nz; ++iz, ++pi) { // slice loop
				cfloat zz = dz * iz - wz;
				cfloat tt = mag * zz;
				cfloat tt_bin = tt / dt + wt;
			// z value is used to determine index of t "it"
				cint it = floorf(tt_bin); // nearest nbr in "t"

				if (it < 0 || it >= nt-1) // out of FOV
					continue;
				else { // reconstructing the image
					cfloat wu = tt_bin - it;
					cfloat wd = 1. - wu;
					cfloat p1 = wl * pp1[it] + wr * pp2[it]; // interpolate
					cfloat p2 = wl * pp1[it+1] + wr * pp2[it+1]; // horizontal

					// final vertical interpolation:
					*pi += w2 * (wu * p1 + wd * p2);
				}
			}
		}
	}
	Ok
}

#endif

