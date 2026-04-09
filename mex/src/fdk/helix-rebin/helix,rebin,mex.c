// helix,rebin,mex.c
// Matlab mex gateway for helical CT rebinning to fan-beam sinograms
// Copyright 2010-07-26, Jeff Fessler, University of Michigan

#include "jf,mex,def.h"
#include "helix,rebin,def.h"

#define Usage	"usage error. see above"


// helix_rebin_mex_help()
static void helix_rebin_mex_help(void)
{
	printf("Usage for helical CT rebin routines.\n");

	printf("\n\
\n\
	[sino orbits] = function('helix,rebin,ssrb',\n\
		nz, dz, offset_z, \n\
		dso, dsd, dfs, [ds dt], [offset_s offset_t], \n\
		pitch, source_z0, orbit, orbit_start, proj, nthread)\n\
\n\
		output sinogram output is single [ns na_rebin nz], half scan\n\
		output orbits is single [nz 2] (orbit, orbit_start)\
\n\
		nz (int32) # of image slices desired\n\
		dz (single) slice spacing\n\
		offset_z: (single) center offset in pixels (usually 0)\n\
			z = dz * (iz - wz), wz = (nz-1)/2 + offset_z\n\
		dso: (single) distance from source to isocenter\n\
		dsd: (single) distance from source to detector\n\
		ds: (single) horizontal ray spacing\n\
		dt: (single) vertical ray spacing\n\
		offset_s: (single) horizontal (channel) offset [pixels]\n\
		offset_t: (single) vertical (axial) offset on detector [pixels]\n\
		pitch: (single) translation of source per 360 rotation\n\
			in units of fraction of axial FOV at isocenter\n\
		source_z0: (single) source z location for 1st projection view\n\
		orbit: (single) [degrees]\n\
		orbit_start: (single) [degrees]\n\
		proj: (single) [ns nt na] projection views\n\
		nthread: (int32) (ignored)\n\
\n");
}


// helix_rebin_ssrb_mex()
static sof helix_rebin_ssrb_mex(
mxArray *plhs[],
Cmx mx_nz, // # of slices
Cmx mx_dz, // slice spacing
Cmx mx_cz, // image slice offset [unitless]
Cmx mx_dso, // distance source to origin (iso-center)
Cmx mx_dsd, // distance source to detector
Cmx mx_dfs, // 0 for 3rd gen CT or inf for flat
Cmx mx_dst, // [2] ds dt
Cmx mx_ost, // [2] offset_s offset_t
Cmx mx_pitch,
Cmx mx_source_z0,
Cmx mx_orbit,
Cmx mx_orbit_start,
Cmx mx_proj, // [ns nt] or [nt ns] or [ns nt na] or [nt ns na]
Cmx mx_nthread)
{
	Call(mxIsScalarInt32, (mx_nz))
	cint nz = mxGetInt(mx_nz);

	Call(mxIsScalarSingle, (mx_dz))
	Call(mxIsScalarSingle, (mx_cz))

	Call(mxIsScalarSingle, (mx_dso))
	cfloat dso = mxGetSingle(mx_dso);

	Call(mxIsScalarSingle, (mx_dsd))
	cfloat dsd = mxGetSingle(mx_dsd);

	Call(mxIsScalarSingle, (mx_dfs))
	cfloat dfs = mxGetSingle(mx_dfs);

	Call(mxIsRealSingle, (mx_dst))
	if (2 != mxGetM(mx_dst) * mxGetN(mx_dst))
		Fail("mx_dst must have 2 elements")
	cfloat ds = mxGetPr_cfloat(mx_dst)[0];
	cfloat dt = mxGetPr_cfloat(mx_dst)[1];

	Call(mxIsRealSingle, (mx_ost))
	if (2 != mxGetM(mx_ost) * mxGetN(mx_ost))
		Fail("mx_ost must have 2 elements")
	cfloat offset_s = mxGetPr_cfloat(mx_ost)[0];
	cfloat offset_t = mxGetPr_cfloat(mx_ost)[1];

	Call(mxIsScalarSingle, (mx_pitch))
	cfloat pitch = mxGetSingle(mx_pitch);

	Call(mxIsScalarSingle, (mx_source_z0))
	cfloat source_z0 = mxGetSingle(mx_source_z0);

	Call(mxIsScalarSingle, (mx_orbit))
	cfloat orbit = mxGetSingle(mx_orbit);

	Call(mxIsScalarSingle, (mx_orbit_start))
	cfloat orbit_start = mxGetSingle(mx_orbit_start);

	Call(mxIsRealSingle, (mx_proj))
	if (3 != mxGetNumberOfDimensions(mx_proj))
		Fail("mx_proj must be 3d")
	cint ns = mxGetDimensions(mx_proj)[0];
	cint nt = mxGetDimensions(mx_proj)[1];
	cint na = mxGetDimensions(mx_proj)[2];

	Call(mxIsScalarInt32, (mx_nthread))
	cint nthread = mxGetInt(mx_nthread);

	float orbit_short = 0;
	int na_rebin;
	Call(helix_rebin_dim, (&na_rebin, &orbit_short,
		dso, dsd, dfs, ds, offset_s, orbit, ns, na))

	mwSize odim[3] = {ns, na_rebin, nz}; // usual order
	Call(plhs[0] = mxCreateNumericArray, (3, odim, mxSINGLE_CLASS, mxREAL))

	Call(plhs[1] = mxCreateNumericMatrix, (nz, 2, mxSINGLE_CLASS, mxREAL));
	float *orbits = (float *) mxGetPr(plhs[1]); // [nz 2]
	for (int iz=0; iz < nz; ++iz)
		orbits[iz] = orbit_short;

	if (nthread != 1)
		Warn("only nthread=1 supported, ignoring")

	Call(helix_rebin_ssrb, (
			(float *) mxGetPr(plhs[0]), // output [ns na_rebin nz]
			orbits + nz, // orbits [nz 2]
			na_rebin,
			nz,
			mxGetSingle(mx_dz),
			mxGetSingle(mx_cz), // offset_z
			dso, dsd, dfs,
			ns, nt, na,
			ds, dt,
			offset_s, offset_t,
			pitch, source_z0,
			orbit, orbit_start,
			mxGetPr_cfloat(mx_proj)))
	Ok
}


// intermediate gateway routine
static sof helix_rebin_mex(
int nlhs, mxArray *plhs[],
int nrhs, Const mxArray *prhs[])
{
	char *arg;

	if (nrhs <= 1 || !mxIsChar(prhs[0])) // expect 'helix,rebin,?', ...
	{
		helix_rebin_mex_help();
		Call(mxu_arg, (nrhs, prhs))
		Fail(Usage)
	}

	if (nrhs == 1 && mxIsChar(prhs[0])) // 'help'
	{
		helix_rebin_mex_help();
		Ok
	}

	Call(arg = mxu_string, (prhs[0], "1st argument"))

	if (Streq(arg, "helix,rebin,ssrb"))
	{
		if (nrhs != 15 || nlhs != 2)
			Fail(Usage)
		Call(helix_rebin_ssrb_mex, (plhs,
			prhs[1], prhs[2], prhs[3], prhs[4], prhs[5],
			prhs[6], prhs[7], prhs[8], prhs[9], prhs[10],
			prhs[11], prhs[12], prhs[13], prhs[14]))
	}

	else
		Fail1("command '%s' unknown", arg)

	Call(mxu_string_free, (arg))
	Ok
}


#ifdef Use_helix_rebin_mex
// gateway routine
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, Const mxArray *prhs[])
{
	if (!nlhs && !nrhs) {
		helix_rebin_mex_help();
		return;
	}
	if (!helix_rebin_mex(nlhs, plhs, nrhs, prhs))
		mexErrMsgTxt("helix_rebin_mex");
}
#endif
