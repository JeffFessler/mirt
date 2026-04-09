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
cfloat dfs, // distance from focal point to source (0 or inf)
cint ns,
cint nt,
cfloat ds, // horizontal ray spacing
cfloat dt, // vertical ray spacing
cfloat offset_s, // channel offset [pixels]
cfloat offset_t, // vertical offset on detector [pixels]
cfloat *proj, // [ns nt] projection view at angle beta
cfloat beta) // source angle [radians]
{
	cint nxy = nx * ny;
	cfloat wx = (nx-1)/2. + offset_x;
	cfloat wy = (ny-1)/2. + offset_y;
	cfloat wz = (nz-1)/2. + offset_z;
	cfloat ws = (ns-1)/2. + offset_s;
	cfloat wt = (nt-1)/2. + offset_t;
	cfloat sinb = sinf(beta);
	cfloat cosb = cosf(beta);

#if 0
	printf("nxyz=%d,%d,%d dxyz=%g,%g,%g offset_xyz=%g,%g,%g\n"
		"beta=%g dso=%g dsd=%g ds=%g dt=%g offset_st=%g,%g\n",
		nx, ny, nz, dx, dy, dz,
		offset_x, offset_y, offset_z,
		beta, dso, dsd, ds, dt, offset_s, offset_t);
#endif

	truf is_par = Isinf(dsd) || Isinf(dso); // parallel-beam?
	truf is_arc = False; // flat detector
	if (!is_par)
	{
		if (dfs == 0)
			is_arc = True; // 3rd-gen CT arc detector
		else if (!Isinf(dfs))
			Warn("dfs not done - junk!")
	}

	// loop over pixels
	for (int iy = 0; iy < ny; ++iy)
	{
		cfloat yy = dy * (iy - wy);
	 for (int ix = 0; ix < nx; ++ix, ++image)
	 {
		// if (mask2[ix + iy*nx] != mask_id)
		if (!mask2[ix + iy*nx])
			continue;

		cfloat xx = dx * (ix - wx);
		cfloat xbeta = xx * cosb + yy * sinb;
		cfloat ybetas = dso - (-xx * sinb + yy * cosb);
		cfloat mag = is_par ? 1. : (dsd / ybetas);
		cfloat ss = is_arc ? (dsd * atan2(xbeta, ybetas)) // arc
				: (mag * xbeta); // flat or par
		cfloat ss_bin = ss / ds + ws;
		cint is = (int) Floorf(ss_bin);

		if (is < 0 || is >= ns-1) // out of FOV
			continue;

		cfloat wr = ss_bin - is;
		cfloat wl = 1. - wr;
		cfloat w2 = is_arc ? // fan-beam image domain weighting
			(Sqr(dsd) / (Sqr(ybetas) + Sqr(xbeta)))
			: Sqr(mag); // flat or par

		for (int iz = 0; iz < nz; ++iz)
		{
			cfloat zz = dz * (iz - wz);
			cfloat tt = mag * zz; // todo: is this correct for flat?
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
