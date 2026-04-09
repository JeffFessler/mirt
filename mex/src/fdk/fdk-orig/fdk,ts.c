// fdk,ts.c
// Feldkamp (aka FDK) cone-beam back-projection for arc/flat detector.
// For detector ordering (t,s): detector row (along axis) varies fastest.
// Copyright 2005-6-27, Jeff Fessler, University of Michigan
//
// 2010-05-07
// added 1-pixel padding to all sides of projection views
// so that data in all detector rows and columns can be used.
// this also (should) help with making results less sensitive to
// floating point variations like -ffast-math by ensuring continuity.
//
// 2010-05-10
// is: 0 1 ... ns-1 with proj[is] at each of those ns values
// with triangular interpolation, possibly nonzero for (-1 < s < ns)
// since is0 = floor(s) and is1 = is0 + 1, we care if -1 <= is0 and is0 <= ns-1

#include "defs-env.h"
#include "def,fdk.h"

// fdk_ts_help()
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
		dx,dy,dz: (single) voxel size\n\
		offset_x,_y,_z: (single) center offset in pixels (usually 0)\n\
		mask2: (uint8) [nx ny] 2D support mask\n\
		dso: (single) distance from source to isocenter\n\
		dsd: (single) distance from source to detector\n\
		dfs: (single) distance from focal point to source (0 or inf)\n\
		ds: (single) horizontal ray spacing\n\
		dt: (single) vertical ray spacing\n\
		offset_s: (single) channel offset [pixels]\n\
		offset_t: (single) vertical offset on detector [pixels]\n\
		nthread: (int32) # of processors\n\
		proj: (single) [nt ns na] (trick!) projection view for each beta\n\
		beta: (single) [na] source angle(s) [radians]\n\
\n");
}


// fdk_ts_back1()
// The back-projection is *added* to the image, so the user must zero it!
sof fdk_ts_back1(
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
cfloat beta, // source angle [radians]
cfloat *view_work) // [(nt+2) (ns+2)] <- trick: padded projection view for beta
{
	cfloat wx = (nx-1)/2. + offset_x;
	cfloat wy = (ny-1)/2. + offset_y;
	cfloat wz = (nz-1)/2. + offset_z;
	cfloat ws = (ns-1)/2. + offset_s;
	cfloat wt = (nt-1)/2. + offset_t;
	cfloat sinb = sinf(beta);
	cfloat cosb = cosf(beta);

	truf is_par = Isinf(dsd) || Isinf(dso); // parallel-beam?
	truf is_arc = False; // flat detector
	if (!is_par)
	{
		if (dfs == 0)
			is_arc = True; // 3rd-gen CT arc detector
		else if (!Isinf(dfs))
			Warn("dfs not done - junk!")
	}

	// loop over pixels in top slice (x-y plane)
	// gpu version probably should parallelize these two loops
	for (int iy = 0; iy < ny; ++iy)
	{
		cfloat yy = dy * (iy - wy);
	for (int ix = 0; ix < nx; ++ix, image += nz)
	{
		if (mask2[ix + iy*nx] != mask_id) // each pthread does its part
			continue;

		cfloat xx = dx * (ix - wx);
		cfloat xbeta = xx * cosb + yy * sinb;
		cfloat ybetas = dso - (-xx * sinb + yy * cosb);
		cfloat mag = is_par ? 1. : (dsd / ybetas);
		cfloat ss = is_arc ? (dsd * atan2(xbeta, ybetas)) // arc
					: (mag * xbeta); // flat or par
		cfloat ss_bin = ss / ds + ws;
		cint is = (int) Floorf(ss_bin); // nearest neighbor in "s"

		// view is padded by a single zero pixel on each side
		if (is < -1 || is >= ns) // out of FOV; see note at top
			continue;

		cfloat w2 = is_arc ? // fan-beam image domain weighting
			(Sqr(dsd) / (Sqr(ybetas) + Sqr(xbeta)))
			: Sqr(mag); // flat or par
		cfloat wr = ss_bin - is; // horizontal bilinear
		cfloat wl = 1. - wr; // interpolation factors

		// the extra "+1" factors below are because of zero padding
		register cfloat *pp1 = view_work + (is+1) * (nt+2) + 1;
		register cfloat *pp2 = view_work + (is+2) * (nt+2) + 1;

		// todo: is this correct for flat detector case?
		cfloat t_inc = mag * dz / dt; // t' increment for each iz++
		int iz_min, iz_max;

		if (dz / dt < 0)
			Fail("todo: dz/dt < 0 not done")
		else
		{
			// design iz_min so that first t' > -1
			// where t' = t_inc * (iz - wz) + wt is unitless
			iz_min = Floorf(wz - (wt+1) / t_inc) + 1;

			// design iz_max so that last t' < nt
			// caution: the range is [iz_min, iz_max), hence "+1"
			iz_max = Ceilf(wz + (nt - wt) / t_inc) - 1 + 1;
		}

		iz_min = Max(iz_min, 0);
		iz_max = Min(iz_max, nz);

		float t_bin = t_inc * (iz_min - wz) + wt; // initial t'
		if (t_bin < -1)
			Warn1("bug t_bin=%g < -1", t_bin)
		else if (t_bin == -1)
		{
//			Warn1("t_bin=%g == -1 (inefficient)", t_bin) // rare
			iz_min++;
			t_bin += t_inc;
		}

		register float *pi = image + iz_min;

		for (int iz = iz_min; iz < iz_max; ) // slice loop
		{
			cint it = (int) Floorf(t_bin); // nearest nbr in "t"

#if 0 // change ranges for debug
//			if (it < 0 || it >= nt-1) // out of FOV
			if (it >= nt) // out of FOV
				Warn4("bug iz_max=%d iz=%d it=%d >= nt=%d",
					iz_max, iz, it, nt)
#endif

			cfloat wu = t_bin - it;
			cfloat wd = 1. - wu;
			cfloat p1 = wl * pp1[it] + wr * pp2[it]; // interpolate
			cfloat p2 = wl * pp1[it+1] + wr * pp2[it+1]; // horizontal

			*pi++ += w2 * // fan-beam weighting:
				(wu * p1 + wd * p2); // vertical interpolation
			++iz;
			t_bin += t_inc;
		}

	} // ix
	} // iy

	Ok
}
