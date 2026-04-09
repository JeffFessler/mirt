// helix,fbp,mex.c
// Matlab mex gateway for helical CT FBP.
// Copyright 2010-07-26, Jeff Fessler, University of Michigan

#include "jf,mex,def.h"
#include "helix,fbp,def.h"

#define Usage	"usage error. see above"

static void helix_fbp_mex_help(void)
{
	printf("Usage for helical CT FBP mex routines.\n");

	printf("\n\
\n\
	image = function('helix,fbp', nx,ny,nz, dx,dy,dz, \n\
		offset_x, offset_y, offset_z, mask2, \n\
		dso, dsd, dfs, ds, dt, offset_s, offset_t, proj, beta, nthread)\n\
\n\
		image output is single [nx ny nz]\n\
		nx,ny,nz: (int32) image size\n\
		dx,dy,dz: (single) voxel size\n\
		offset_x,_y,_z: (single) center offset in pixels (usually 0)\n\
		mask2: (uint8) [nx,ny] 2D support mask\n\
		dso: (single) distance from source to isocenter\n\
		dsd: (single) distance from source to detector\n\
		ds: (single) horizontal ray spacing\n\
		dt: (single) vertical ray spacing\n\
		offset_s: (single) channel offset [pixels]\n\
		offset_t: (single) vertical offset on detector [pixels]\n\
		pitch: (single) todo: document\n\
		source_z0: (single) todo: document\n\
		orbit: (single) [degrees]\n\
		orbit_start: (single) [degrees]\n\
		proj: (single) [ns nt na] projection view at angle beta\n\
		nthread: (int32)\n\
\n");
}


static sof helix_fbp_mex_back(
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
Cmx mx_pitch,
Cmx mx_source_z0,
Cmx mx_orbit,
Cmx mx_orbit_start,
Cmx mx_proj, // [ns nt] or [nt ns] or [ns nt na] or [nt ns na]
Cmx mx_nthread)
{
	Call(mxIsInt32, (mx_nxyz))
	if (3 != mxGetM(mx_nxyz) * mxGetN(mx_nxyz))
		Fail("mx_nxyz must have 3 elements")
	cint nx = ((cint *) mxGetPr(mx_nxyz))[0];
	cint ny = ((cint *) mxGetPr(mx_nxyz))[1];
	cint nz = ((cint *) mxGetPr(mx_nxyz))[2];

	Call(mxIsRealSingle, (mx_dxyz))
	if (3 != mxGetM(mx_dxyz) * mxGetN(mx_dxyz))
		Fail("mx_dxyz must have 3 elements")

	Call(mxIsRealSingle, (mx_cxyz))
	if (3 != mxGetM(mx_cxyz) * mxGetN(mx_cxyz))
		Fail("mx_cxyz must have 3 elements")

	Call(mxIsUint8, (mx_mask))
	if (2 != mxGetNumberOfDimensions(mx_mask))
		Fail("mask must be 2d")
	if (nx != (int) mxGetM(mx_mask) || ny != (int) mxGetN(mx_mask))
		Fail("mask size not nx by ny?")

	Call(mxIsScalarSingle, (mx_dso))
	Call(mxIsScalarSingle, (mx_dsd))
	Call(mxIsScalarSingle, (mx_dfs))

	Call(mxIsRealSingle, (mx_dst))
	if (2 != mxGetM(mx_dst) * mxGetN(mx_dst))
		Fail("mx_dst must have 2 elements")

	Call(mxIsRealSingle, (mx_ost))
	if (2 != mxGetM(mx_ost) * mxGetN(mx_ost))
		Fail("mx_ost must have 2 elements")

	Call(mxIsScalarSingle, (mx_pitch))
	Call(mxIsScalarSingle, (mx_source_z0))

	Call(mxIsScalarSingle, (mx_orbit))
	Call(mxIsScalarSingle, (mx_orbit_start))

	Call(mxIsRealSingle, (mx_proj))
	if (3 != mxGetNumberOfDimensions(mx_proj))
		Fail("mx_proj must be 3d")

	Call(mxIsScalarInt32, (mx_nthread))
	cint nthread = mxGetInt(mx_nthread);

	mwSize odim[3] = {nx, ny, nz}; // usual order
	Call(plhs[0] = mxCreateNumericArray, (3, odim, mxSINGLE_CLASS, mxREAL))

//	cint chat = 0;

	byte *mask2 = NULL;
//	Mem0(ig->mask2, nx * ny)
//	Call(cbct_mask_init, (ig->mask2, (cbyte *) mxGetData(mx_mask),
//		nx, ny, nthread, chat))

	if (nthread != 1)
		Warn("only nthread=1 supported, ignoring")

	Call(helix_fbp_ssrb, (
			(float *) mxGetPr(plhs[0]), // output [nx ny nz]
			nx, ny, nz,
			mxGetPr_cfloat(mx_dxyz)[0], // dx
			mxGetPr_cfloat(mx_dxyz)[1], // dy
			mxGetPr_cfloat(mx_dxyz)[2], // dz
			mxGetPr_cfloat(mx_cxyz)[0], // offset_x
			mxGetPr_cfloat(mx_cxyz)[1], // offset_y
			mxGetPr_cfloat(mx_cxyz)[2], // offset_z
			mask2,
			mxGetSingle(mx_dso),
			mxGetSingle(mx_dsd),
			mxGetSingle(mx_dfs),
			mxGetDimensions(mx_proj)[0], // ns 
			mxGetDimensions(mx_proj)[1], // nt
			mxGetDimensions(mx_proj)[2], // na
			mxGetPr_cfloat(mx_dst)[0], // ds
			mxGetPr_cfloat(mx_dst)[1], // dt
			mxGetPr_cfloat(mx_ost)[0], // offset_s
			mxGetPr_cfloat(mx_ost)[1], // offset_t
			mxGetSingle(mx_pitch),
			mxGetSingle(mx_source_z0),
			mxGetSingle(mx_orbit),
			mxGetSingle(mx_orbit_start),
			mxGetPr_cfloat(mx_proj)))

//	Free0(mask2)
	Ok
}


// intermediate gateway routine
static sof helix_fbp_mex(
int nlhs, mxArray *plhs[],
int nrhs, Const mxArray *prhs[])
{
	char *arg;

	if (nrhs <= 1 || !mxIsChar(prhs[0])) // expect 'helix,fbp', ...
	{
		helix_fbp_mex_help();
		Call(mxu_arg, (nrhs, prhs))
		Fail(Usage)
	}

	if (nrhs == 1 && mxIsChar(prhs[0])) // 'help'
	{
		helix_fbp_mex_help();
		Ok
	}

	Call(arg = mxu_string, (prhs[0], "1st argument"))

	if (Streq(arg, "helix,fbp,ssrb"))
	{
		if (nrhs != 16 || nlhs != 1)
			Fail(Usage)
		Call(helix_fbp_mex_back, (plhs,
			prhs[1], prhs[2], prhs[3], prhs[4], prhs[5],
			prhs[6], prhs[7], prhs[8], prhs[9], prhs[10],
			prhs[11], prhs[12], prhs[13], prhs[14], prhs[15]))
	}

	else
		Fail1("command '%s' unknown", arg)

	Call(mxu_string_free, (arg))
	Ok
}


#ifdef Use_helix_fbp_mex
// gateway routine
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, Const mxArray *prhs[])
{
	if (!nlhs && !nrhs) {
		helix_fbp_mex_help();
		return;
	}
	if (!helix_fbp_mex(nlhs, plhs, nrhs, prhs))
		mexErrMsgTxt("helix_fbp_mex");
}
#endif
