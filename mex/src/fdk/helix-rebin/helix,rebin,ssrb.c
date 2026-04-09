// helix_rebin_ssrb.c
// Single slice rebinning method for reconstructing cone-beam tomography data
// collected with a helical source trajectory.
// 2010-07-15 initial version by Gregory Handy
// 2010-09-12 Jeff Fessler: numerous modifications and corrections

#include "helix,rebin,def.h"


// helix_rebin_dim()
sof helix_rebin_dim(
int *p_na1,
float *p_orbit_short,
cfloat dso,
cfloat dsd,
cfloat dfs,
cfloat ds,
cfloat offset_s,
cfloat orbit,
cint ns,
cint na)
{
	(void) dso;

	cfloat na1_float = na / orbit * 360.; // # views in 1 turn
	if (Abs(round(na1_float) - na1_float) > 1e-4)
		Fail("bad na/orbit; need integral # per rotation")

	// find gamma_max
	cfloat ws = (ns-1)/2. + offset_s;
	cfloat s0 = ds * (0 - ws);
	cfloat s1 = ds * (ns-1 - ws);
	cfloat smax = Max(Abs(s0), Abs(s1));

	float gamma_max;
	if (dfs == 0) // arc
		gamma_max = smax / dsd;
	else if (Isinf(dfs)) // flat
		gamma_max = atan(smax / dsd);
	else
		Fail("bad dfs")

	cfloat fan_angle = gamma_max * 180. / M_PI;
	cfloat orbit_short_ideal = 180. + 2 * fan_angle;
	int na1 = 2 * ceilf(orbit_short_ideal / (orbit/na) / 2); // make even

	if (na1 == na + 2)
	{
		Warn2("na1=%d na=%d is off by 2; assuming short scan", na1, na)
		na1 = na;
	}

	if (na1 > na)
		Fail4("na1 %d > na %d? osi=%g fan=%g",
			na1, na, orbit_short_ideal, fan_angle)

	*p_na1 = na1;

	cfloat orbit_short_actual = na1 * (orbit / na);
	*p_orbit_short = orbit_short_actual;

	Ok
}


// rebin_helix_scale()
static float rebin_helix_scale(cfloat dfs, cfloat dsd, cfloat ss, cfloat tt)
{
	if (Isinf(dfs)) // flat noo:99:ssr (1)
		return sqrtf(Sqr(ss) + Sqr(tt) + Sqr(dsd))
			/ sqrtf(Sqr(tt) + Sqr(dsd));
	// case 0 % arc
		return dsd / sqrtf(Sqr(tt) + Sqr(dsd));
}



// helix_rebin_ssrb_z()
// make short-scan sinogram for a single slice
sof helix_rebin_ssrb_z(
float *sino, // [ns na_rebin]
float *p_z_orbit_start,
cint na_rebin,
cfloat zmid,
cfloat dso,
cfloat dsd,
cfloat dfs,
cint ns,
cint nt,
cint na,
cfloat ds,
cfloat dt,
cfloat offset_s,
cfloat offset_t,
cfloat pitch,
cfloat source_z0, // first z position of source
cfloat orbit,
cfloat orbit_start,
cfloat *proj) // [ns nt na] I: projection
{
	cfloat ws = (ns - 1) / 2. + offset_s;
	cfloat wt = (nt - 1) / 2. + offset_t;

	cfloat na_per_360 = na * (360. / orbit); // # views per turn
	cfloat zfov = Isinf(dsd) || Isinf(dso) ? nt * dt : dso / dsd * nt * dt;
	cfloat source_dz = pitch * zfov / na_per_360;

	// find view index (ia) for closest source position 
	// source_zs[ia] = source_z0 + ia * source_dz
	int ia_middle;
	if (!source_dz) // axial
		ia_middle = ceilf(na / 2.);
	else
		ia_middle = roundf((zmid - source_z0) / source_dz);
//	Note3("zmid=%g ia_middle=%d source_z=%g",
//		zmid, ia_middle, source_z0 + ia_middle * source_dz)

	cint na1_half = na_rebin / 2; // because even
	int ia_start = ia_middle - na1_half;
	if (ia_start < 0)
		ia_start = 0; // resort to first views
	if (ia_start + na_rebin > na)
		ia_start = na - na_rebin; // resort to last views

	*p_z_orbit_start = orbit_start + ia_start * orbit / na;

	for (int ia_new=0; ia_new < na_rebin; ++ia_new)
	{
		cint ia = ia_new + ia_start;
		cfloat *view = proj + ia * ns * nt; // [ns nt]
		cfloat source_z = source_z0 + ia * source_dz;
		cfloat zdiff = zmid - source_z;

		// which rows (t coordinate) of helical projection views
		// based on point where CB ray hits fan-beam plane
		for (int is=0; is < ns; ++is)
		{
			cfloat ss = (is - ws) * ds;

			// linear interpolation

			float tt;
			if (Isinf(dfs)) // flat noo:99:ssr (2)
				tt = zdiff * (Sqr(dsd) + Sqr(ss)) / (dsd * dso);
			else // arc
			{
				cfloat gamma = ss / dsd;
				tt = zdiff * (dsd / dso) / cosf(gamma);
			}

			float itr = tt/dt + wt; // float
			itr = Max(itr, 0); // extrapolate using first row
			itr = Min(itr, nt-1); // extrapolate using last row

			tt = (itr - wt) * dt; // trick: actual t value!
			int it0 = floorf(itr); // nearest neighbor
			it0 = Min(it0, nt-2); // so that it0+1 is ok
			cfloat frac = itr - it0; // for linear interpolation

			cfloat scale = rebin_helix_scale(dfs, dsd, ss, tt);

			cint ii0 = is + it0*ns;
			cfloat tmp0 = view[ii0];
			cfloat tmp1 = view[ii0 + ns];
//			view = (1 - frac) .* tmp0 + frac .* tmp1;
			cfloat val = tmp0 + frac * (tmp1 - tmp0);

			sino[is + ia_new * ns] = scale * val;
//			used([it0+1; it0+2], ia1) = 1;
		}
	}

	Ok

#if 0
	// greg below here

	double num_turns = orbit/360;
	int num_betas = floor(na/num_turns);
	double myPitch = pitch * nt * (dso / dsd) * dt;
	double currentZ;
	double upper_limit, lower_limit, deltaZ, spoints[ns], tpoints;
	double scale, alpha;
	double t_index, y_0, y_1;
	int ia, iz, is, zindex, x_0, x_1;
	double betas[2];
	double z_slices[nz];
	double *fan_beam_proj;
	double new_orbit;
	double change_betas;
	int orbit_found;
	double zloc;
	int angle, nextAngle;
	double smax=0, rmax;
	double image_z0;

	Mem0(fan_beam_proj, ns*num_betas*nz);

	for (is = 0; is < ns; ++is)
	{
		spoints[is]= ds *( is - ((ns-1)/2+offset_s));

		if (smax < abs(spoints[is]))
			smax = abs(spoints[is]);
	}

	if (dfs == 0)
		rmax = dso * sin(smax/dsd);
	else
		rmax = dso * sin(atan(smax / dsd));

	double delta = asin(rmax/dso);
	double dist = .5 * myPitch * (M_PI + 2 * delta)/(2*M_PI);
	int newNA = na_rebin;

	// used to calculate change_betas farther down
	betas[0] = orbit_start * (M_PI/180);
	betas[1] = (orbit_start + orbit * 1/na)*(M_PI/180);

	// z locations of the z-slices
	image_z0 = (0-((nz-1)/2+offset_z))*dz;
	for (iz=0; iz < nz; ++iz)
		{
			z_slices[iz] = image_z0 + iz*dz;
		}

	// loop over the different view angles
	for (ia = 0; ia < na; ++ia)
		{
			currentZ = source_z0 + ia * (myPitch/na) * num_turns;
			upper_limit = currentZ + dist;
			lower_limit = currentZ - dist;

			zindex = 0;

			// enter into the acceptable range for the z-slices
			while ((zindex < nz) && (z_slices[zindex] < lower_limit))
	{
		zindex = zindex + 1;
	}

			// loop over the acceptable z-slices for the current view angle
			while ((zindex < nz) && (z_slices[zindex] < upper_limit))
	{
		deltaZ = z_slices[zindex]-currentZ;

		// Calculate the values of t to be used for each s

		for (is = 0; is<ns; ++is)
			{
				tpoints = ((spoints[is]*spoints[is]+dsd*dsd)/ (dso*dsd)) * deltaZ;
				t_index = (tpoints / dt) + ((nt-1)/ 2+ offset_t);
				x_0 = floor(t_index);
				x_1 = 1 + x_0;

				// Scaling factor due to different ray lenghts
				scale = sqrt((spoints[is]*spoints[is]+dsd*dsd))/sqrt(spoints[is]*spoints[is]+tpoints*tpoints+dsd*dsd);

				alpha = t_index - x_0;
				y_1 = proj[is+ns*x_1+ns*nt*ia];
				y_0 = proj[is+ns*x_0+ns*nt*ia];
				fan_beam_proj[is+ns*(ia%num_betas)+ns*num_betas*zindex] = scale*(alpha*y_1+(1-alpha)*y_0);
			}

		zindex = zindex + 1;
	}
		}

	// Next portion of the code prepares the sinogram for FB2

	// decide on beginning orbit for each z-slice
	new_orbit = orbit_start*(M_PI/180);

	zloc = source_z0;

	change_betas = betas[1]-betas[0];
	for (iz = 0; iz < nz; ++iz)
		{
			orbit_found = 0;
			while (orbit_found == 0)
	{
		if (z_slices[iz] < zloc + dist)
			{
				z_orbit_start[iz] = new_orbit * (180/M_PI);
				orbit_found = 1;
			}
		else
			{
				new_orbit = new_orbit + change_betas;
				zloc = zloc + (myPitch/na) * num_turns;
			}
	}
		}

	// deletes the empty rows of the sinograms and reorganizes the rows
	for(iz = 0; iz<nz; ++iz)
		{
			angle = round(((z_orbit_start[iz]-orbit_start)*(M_PI/180)) / change_betas);
			for (ia = 0; ia<newNA; ++ia)
	{
		nextAngle = angle + ia;
		for (is=0; is<ns; ++is)
			{
				sino[is+ns*ia+ns*newNA*iz] = fan_beam_proj[is+ns*(nextAngle%num_betas)+ns*num_betas*iz];
			}
	}
		}

	// should apply Parker Weight to each sinogram

	Free0(fan_beam_proj);
	Ok
#endif
}


// helix_rebin_ssrb()
sof helix_rebin_ssrb(
float *sino, // [ns na_rebin nz]
float *z_orbit_start, // [nz]
cint na_rebin,
cint nz,
cfloat dz,
cfloat offset_z,
cfloat dso,
cfloat dsd,
cfloat dfs,
cint ns,
cint nt,
cint na,
cfloat ds, // horizontal ray spacing
cfloat dt, // vertical ray spacing)
cfloat offset_s, // channel offset [pixels]
cfloat offset_t, // vertical offset on detector [pixels]
cfloat pitch,
cfloat source_z0, // first z position of source
cfloat orbit, // total rotation angle [degrees]
cfloat orbit_start, // starting angle [degrees]
cfloat *proj) // [ns nt na] I: projection
{
	cfloat wz = (nz - 1) / 2. + offset_z;

	if (na_rebin > na)
		Fail("na_rebin > na")

	for (int iz=0; iz < nz; ++iz)
	{
		cfloat zmid = (iz - wz) * dz; // z coordinate of slice center

		Call(helix_rebin_ssrb_z, (
			sino + iz * ns * na_rebin, // [ns na_rebin]
			z_orbit_start + iz,
			na_rebin,
			zmid,
			dso, dsd, dfs,
			ns, nt, na,
			ds, dt, offset_s, offset_t,
			pitch, source_z0,
			orbit, orbit_start,
			proj))
	}

	Ok
}
