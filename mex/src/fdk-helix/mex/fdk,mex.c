// fdk,mex.c
// Matlab mex gateway for FDK.
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


// fdk_mex_back()
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
Cmx mx_source_zs, // [1] or [na] by DK
Cmx mx_na1_half, // [1] na1_half by DK
Cmx mx_cone_par, // [1] cone_par by DK
Cmx mx_w3d, // [1] w3d by DK
Cmx mx_source_dz_per_rad, // [1] source_dz_per_rad by DK
Cmx mx_pitch, // [1] pitch by DK
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
	if (2 == mxGetNumberOfDimensions(mx_proj))
	{
		na = 1;
		Call(mxIsScalarDouble, (mx_beta))
		Call(mxIsScalarDouble, (mx_source_zs)) // DK
	}
	else if (3 == mxGetNumberOfDimensions(mx_proj))
	{
		na = mxGetDimensions(mx_proj)[2];
		Call(mxIsRealDouble, (mx_beta))
		Call(mxIsRealDouble, (mx_source_zs)) // DK
		int na_beta = (int) mxGetM(mx_beta) * mxGetN(mx_beta);
		if (na != na_beta)
			Fail2("beta length=%d vs na=size(proj,3)=%d", na_beta, na)
	}
	else
		Fail("mx_proj must be 2d or 3d")

	Call(mxIsScalarInt32, (mx_nthread))
	cint nthread = mxGetInt(mx_nthread);

	if (is_st)
	{
		mwSize odim[3] = {nx, ny, nz}; // usual order
		Call(plhs[0] = mxCreateNumericArray, (3, odim, mxSINGLE_CLASS, mxREAL))
		ns = mxGetDimensions(mx_proj)[0];
		nt = mxGetDimensions(mx_proj)[1];

	}
	else
	{
		mwSize odim[3] = {nz, nx, ny}; // trick: zxy order!
		Call(plhs[0] = mxCreateNumericArray, (3, odim, mxSINGLE_CLASS, mxREAL))
		nt = mxGetDimensions(mx_proj)[0];
		ns = mxGetDimensions(mx_proj)[1];
		// Note2("nt=%d ns=%d", nt, ns)
	}

	Call(mxIsScalarInt32, (mx_na1_half)) // DK
        cint na1_half = mxGetInt(mx_na1_half); // DK
	Call(mxIsScalarInt32, (mx_cone_par)) // DK
	cint cone_par = mxGetInt(mx_cone_par); // DK
	Call(mxIsScalarInt32, (mx_w3d)) // DK
        cint w3d = mxGetInt(mx_w3d); // DK
	Call(mxIsRealDouble, (mx_source_dz_per_rad)) // DK
        cfloat source_dz_per_rad = mxGetDouble(mx_source_dz_per_rad); // DK
	Call(mxIsRealSingle, (mx_pitch)) // DK
        cfloat pitch = mxGetSingle(mx_pitch); // DK

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

	if (is_st)
	{
		if (nthread != 1){
			Warn("only nthread=1 supported for fdk_st_, ignoring")
			Warn("source_zs not supported for fdk_st_, igonring") // DK
		}

		for (int ia=0; ia < na; ++ia)
			fdk_st_back1(
				(float *) mxGetData(plhs[0]),
				ig->nx, ig->ny, ig->nz,
				ig->dx, ig->dy, ig->dz,
				ig->offset_x, ig->offset_y, ig->offset_z,
				ig->mask2,
				cg->dso, cg->dsd, cg->dfs,
				cg->ns, cg->nt,
				cg->ds, cg->dt, cg->offset_s, cg->offset_t,
				mxGetPr_cfloat(mx_proj) + ia * nt * ns,
				mxGetPr_cdouble(mx_beta)[ia]);

	}

	else // new faster way
	{
		Call(fdk_ts_back_t, (
			(float *) mxGetData(plhs[0]),
			ig, cg, na,
			mxGetPr_cfloat(mx_proj),
			mxGetPr_cdouble(mx_beta),
			mxGetPr_cdouble(mx_source_zs), // DK
			na1_half, cone_par, w3d, source_dz_per_rad, pitch, // DK
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

	if (nrhs <= 1 || !mxIsChar(prhs[0])) // expect 'fdk,ts,back', ...
	{
		fdk_mex_help();
		Call(mxu_arg, (nrhs, prhs))
		Fail(Usage)
	}

	if (nrhs == 1 && mxIsChar(prhs[0])) // 'help'
	{
		fdk_mex_help();
		Ok
	}

	Call(arg = mxu_string, (prhs[0], "1st argument"))

	// fdk,{st|ts},back
	if (Streq(arg, "fdk,st,back") || Streq(arg, "fdk,ts,back"))
	{
		ctruf is_st = Streqn(arg, "fdk,st", 6);
		if (nrhs != 19 || nlhs != 1) // DK
			Fail(Usage)
		Call(fdk_mex_back, (plhs,
			prhs[1], prhs[2], prhs[3], prhs[4], prhs[5],
			prhs[6], prhs[7], prhs[8], prhs[9], prhs[10],
			prhs[11], prhs[12], 
			prhs[13], prhs[14], prhs[15], prhs[16], // DK
			prhs[17], prhs[18], // DK
			is_st))
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
	if (!nlhs && !nrhs)
	{
		fdk_mex_help();
		return;
	}
	if (!fdk_mex(nlhs, plhs, nrhs, prhs))
		mexErrMsgTxt("fdk_mex");
}
#endif
