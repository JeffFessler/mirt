/*
 * jmh,sse,mex.c mex interface to SSE projection / backprojection by
 * Jens Munk Hansen

 * $Id: jmh_sse_mex.cpp,v 1.25 2007/08/07 21:04:32 jmh Exp $
 */

/*
 * TODO:
 */
#ifdef _MSC_VER
 #pragma warning( disable : 4996 ) // suppress warning: 'sprintf' was declared deprecated
#endif

#define Need_jmh_sse_mex_gateway 1

#include "def,mexarg.h"

#include "tomo_mex.h"
#include "tomo_lut.h"

#include <math.h>

#define Usage   "usage error. see above"

static void jmh_sse_mex_help(void)
{
  printf("Usage for jmh,sse (sse-based projection) mex routines.\n\
\tproj = function('jmh,sse,proj', ...)\n        \
out:\n\
\tproj is [nb,na] single output sinogram\n\
in:\n\
\tnthread         int32, # threads\n\
\tnb              int32, # of radial samples\n\
\toffset_r        single radial sample offset [unitless]\n\
\tdr              single, ray spacing\n\
\tdx              single, pixel spacing in x, can be negative\n\
\tdy              single, pixel spacing in y, ''\n\
\toffset_x        single, x sample offset [unitless], e.g., 0 or 0.5\n\
\toffset_y        single, y ''\n\
\tmask    [nx,ny] uint8, support mask\n\
\tna              int32, # angles\n\
\torbit_low       single, projection angle [radians]\n  \
\torbit_high      single, projection angle [radians]\n\
\tnsubset         int32, # of subsets, the number of angles per subset must be divisible by 4\n\
\tisubset         int32, subset index, isubset = [0,....,nsubset-1]\n\
\timage   [nx,ny] single, input image (pixel coefficients)\n    \
\n\
\tback = function('jmh,sse,back', ...)\n        \
out:\n\
\tback is [nx,ny] single\n\
in:\n\
\teverything else as above, except:\n\
\tproj    [nb,na] single, input sinogram\n\
\n\
\tmatrix = function('jmh,sse,matrix',...)\n\
out:\n\
\tmatrix is [nx*ny*na*M,nx*ny*na*M] values and indices for system matrix\n\
in:\n\
\tnthread         int32, # threads\n\
\torbit_low       single, projection angle [radians]\n\
\torbit_high      single, projection angle [radians]\n\
\tnb              int32, # of radial samples\n\
\tdr              single, ray spacing\n\
\toffset_r        single, radial sample offset [unitless]\n\
\tsw              single, strip width units of dx\n\
\tnx              int32, # pixels in x\n\
\tdx              single, pixels spacing in x, must be positive\n\
\toffset_x        single, x sample offset [unitless], e.g., 0 or 0.5\n\
\tny              int32, # pixels in y\n\
\tdy              single, pixels spacing in y can be negative\n\
\toffset_y        single, y ''\n\
\tmask    [nx,ny] uint8, support mask\n\
\tna              int32, # angles\n\
\n");
}

char *mxu_string(const mxArray *mx, const char *arg) {
  char *cbuffer;
  int n = mxGetM(mx) * mxGetN(mx) + 1;

  if (!mxIsChar(mx))
    Fail1("%s must be char array", arg)

      Call(cbuffer = (char *) mxCalloc, (n, sizeof(char)))

      if (mxGetString(mx, cbuffer, n))
        Warn("bug with mxGetString")
          return cbuffer;
}

bool mxu_string_free(char *arg) {
  mxFree((void*)arg);
  Ok
}

/*
 * jmh_sse_proj_matrix() 
 */
static bool jmh_sse_matrix(mxArray *plhs[],
                           const mxArray* mx_nthread,
                           const mxArray* mx_orbit_low,
                           const mxArray* mx_orbit_high,
                           const mxArray* mx_nb,
                           const mxArray* mx_dr,
                           const mxArray* mx_offset_r,
                           const mxArray* mx_sw,
                           const mxArray* mx_nx,
                           const mxArray* mx_dx,
                           const mxArray* mx_offset_x,
                           const mxArray* mx_ny,
                           const mxArray* mx_dy,
                           const mxArray* mx_offset_y,
                           const mxArray* mx_mask,
                           const mxArray* mx_na) {

  // Change to unsigned and mxIsScalarInt32 to mxIsScalarUint32
  int na, nb, nx, ny, nthread;
  float dx, dy, dr;
  float offset_x, offset_y, offset_r;
  float sw;
  float orbit_low, orbit_high;
  int o_dims1[1];
  int o_dims2[1];
  
  /*
   * Check input types and sizes
   */
  
  int norb;
  const float* p_orbit;
  
  Call(mxIsScalarInt32, (mx_nthread))
  nthread = mxGetInt(mx_nthread);

  Call(mxIsScalarSingle, (mx_orbit_low))
  orbit_low = *((const float *) mxGetData(mx_orbit_low));

  Call(mxIsScalarSingle, (mx_orbit_high))
  orbit_high = *((const float *) mxGetData(mx_orbit_high));

  Call(mxIsScalarInt32, (mx_nb))
  nb = mxGetInt(mx_nb);
  
  Call(mxIsScalarSingle, (mx_dr))
  dr = *((const float *) mxGetData(mx_dr));

  Call(mxIsScalarSingle, (mx_offset_r))
  offset_r = *((const float *) mxGetData(mx_offset_r));

  Call(mxIsScalarSingle, (mx_sw))
  sw = *((const float *) mxGetData(mx_sw));

  Call(mxIsScalarInt32, (mx_nx))
  nx = mxGetInt(mx_nx);
  
  Call(mxIsScalarSingle, (mx_dx))
  dx = *((const float *) mxGetData(mx_dx));

  Call(mxIsScalarSingle, (mx_offset_x))
  offset_x = *((const float *) mxGetData(mx_offset_x));

  Call(mxIsScalarInt32, (mx_ny))
  ny = mxGetInt(mx_ny);

  Call(mxIsScalarSingle, (mx_dy))
  dy = *((const float *) mxGetData(mx_dy));

  Call(mxIsScalarSingle, (mx_offset_y))
  offset_y = *((const float *) mxGetData(mx_offset_y));

  Call(mxIsUint8, (mx_mask))
  if (nx != (int) mxGetM(mx_mask) || ny != (int) mxGetN(mx_mask))
    Fail("image / mask size mismatch")

  Call(mxIsScalarInt32, (mx_na))
  na = mxGetInt(mx_na);

#if 0
  if (1) {
    Note1("nthread=%d", nthread)
    Note1("orbit_low=%g", orbit_low)
    Note1("orbit_high=%g", orbit_high)
    Note1("na=%d", na)
    Note1("nb=%d", nb)
    Note1("dr=%g", dr)      
    Note1("sw=%g", sw)      
    Note1("offset_r=%g", offset_r)
    Note1("offset_x=%g", offset_x)
    Note1("offset_y=%g", offset_y)
    Note4("nx,ny=%d,%d dx,dy=%g,%g", nx, ny, dx, dy)
  }
#endif

  unsigned int M = (unsigned int) ceil((dx * sqrt(2.0f) + sw) / dr);

  int n_odim = 1;
  
  o_dims1[0] = nx*ny*na*M;
  o_dims2[0] = nx*ny*na*M;

  Call(plhs[0] = mxCreateNumericArray,
       (n_odim, (const mwSize*)o_dims1, mxSINGLE_CLASS, mxREAL))

#ifndef _LARGE_INDICES
  Call(plhs[1] = mxCreateNumericArray,
       (n_odim, (const mwSize*)o_dims2, mxUINT16_CLASS, mxREAL))
#else
  Call(plhs[1] = mxCreateNumericArray,
       (n_odim, (const mwSize*)o_dims2, mxUINT32_CLASS, mxREAL))
#endif
    
  // Initialize common variables - necessary
  LUT::initialize_LUTs(na, 1, orbit_low,orbit_high, dx, dr, sw);
  
#ifndef _LARGE_INDICES
  Call(tomo_strip_parallel_2d_mex, (na,
                                    nb, dr, offset_r, sw,
                                    nx, dx, offset_x,
                                    ny, dy, offset_y,
                                    (cbyte *) mxGetData(mx_mask), ((float *) mxGetData(plhs[0])), ((unsigned short*) mxGetData(plhs[1]))))
#else
  Call(tomo_strip_parallel_2d_mex, (na,
                                    nb, dr, offset_r, sw,
                                    nx, dx, offset_x,
                                    ny, dy, offset_y,
                                    (cbyte *) mxGetData(mx_mask), ((float *) mxGetData(plhs[0])), ((unsigned long*) mxGetData(plhs[1]))))
#endif
    
  LUT::free_LUTs();
  
  Ok
}


/*
 * jmh_sse_proj_mex()
 */
static bool jmh_sse_proj_mex(mxArray *plhs[],
                             Cmx mx_nthread,
                             Cmx mx_nb,
                             Cmx mx_offset_r,
                             Cmx mx_dr,
                             Cmx mx_dx,
                             Cmx mx_dy,
                             Cmx mx_offset_x,
                             Cmx mx_offset_y,
                             Cmx mx_mask,
                             Cmx mx_na,
                             Cmx mx_orbit_low,
                             Cmx mx_orbit_high,
                             Cmx mx_nsubset,
                             Cmx mx_isubset,
                             Cmx mx_image) {
  int ll, LL = 1;
  int nx, ny, nb, na, nthread;
  float dx, dy;
  float offset_x, offset_y, offset_r;
  float dr;
  float orbit_low, orbit_high;
  int nsubset, isubset;



  /*
   * Check input types and sizes
   */
  Call(mxIsScalarInt32, (mx_nthread))
  nthread = mxGetInt(mx_nthread);

  Call(mxIsScalarInt32, (mx_nb))
  nb = mxGetInt(mx_nb);

  Call(mxIsScalarSingle, (mx_offset_r))
  offset_r = *((const float *) mxGetData(mx_offset_r));

  Call(mxIsScalarSingle, (mx_dr))
  dr = *((const float *) mxGetData(mx_dr));
  
  Call(mxIsScalarSingle, (mx_dx))
  dx = *((const float *) mxGetData(mx_dx));

  Call(mxIsScalarSingle, (mx_dy))
  dy = *((const float *) mxGetData(mx_dy));

  Call(mxIsScalarSingle, (mx_offset_x))
  offset_x = *((const float *) mxGetData(mx_offset_x));

  Call(mxIsScalarSingle, (mx_offset_y))
  offset_y = *((const float *) mxGetData(mx_offset_y));


  Call(mxIsRealSingle, (mx_image))
  nx = mxGetM(mx_image);
  ny = mxGetDimensions(mx_image)[1];
      
      /* trick: can't use mxGetN since it returns product */

  Call(mxIsUint8, (mx_mask))
  if (nx != (int) mxGetM(mx_mask) || ny != (int) mxGetN(mx_mask))
    Fail("image / mask size mismatch")

  Call(mxIsScalarInt32, (mx_na))
  na = mxGetInt(mx_na);

  Call(mxIsScalarSingle, (mx_orbit_low))
  orbit_low = *((const float *) mxGetData(mx_orbit_low));

  Call(mxIsScalarSingle, (mx_orbit_high))
  orbit_high = *((const float *) mxGetData(mx_orbit_high));

  Call(mxIsScalarInt32, (mx_nsubset))
  nsubset = mxGetInt(mx_nsubset);

  Call(mxIsScalarInt32, (mx_isubset))
  isubset = mxGetInt(mx_isubset);

  if ((isubset < 0) || !(isubset < nsubset)) Fail("invalid subset")

  if ((na / nsubset) % 4 != 0)
    Fail("(na/nsubset) % 4 != 0")                                                   
                                               
#if 0
  if (1) {
    Note1("nthread=%d", nthread)
    Note2("isubset,nsubset=%d,%d", isubset,nsubset)
    Note1("na=%d", na)
    Note1("nb=%d", nb)
    Note1("offset_r=%g", offset_r)
    Note1("offset_x=%g", offset_x)
    Note1("offset_y=%g", offset_y)
    Note4("nx,ny=%d,%d dx,dy=%g,%g", nx, ny, dx, dy)
  }
#endif
   
  /*
   * output dimensions; allowing for multiple input images
   */
  {
    int n_idim, *o_dims;
    cint *i_dims;
    n_idim = mxGetNumberOfDimensions(mx_image);
    Call(i_dims = (cint*) mxGetDimensions, (mx_image))

    //  Mem0(o_dims, n_idim)
    Alloc0(o_dims, int, n_idim, sizeof(*(o_dims)))
  
    Bcopy(i_dims, o_dims, n_idim)
    o_dims[0] = nb;
    o_dims[1] = na/nsubset;
    for (ll=2; ll < n_idim; ++ll)
      LL *= i_dims[ll];

    /* allocate output */
    Call(plhs[0] = mxCreateNumericArray,
         (n_idim, (const mwSize*)o_dims, mxSINGLE_CLASS, mxREAL))
    Free0(o_dims)
  }

  // Initialize common variables - necessary
  LUT::initialize_LUTs(na, nsubset, orbit_low, orbit_high, dx, dr, dr); // dx,dr,sw
  //    mexPrintf("LL: %d\n",LL);
  for (ll=0; ll < LL; ++ll) {
    Call(jmh_sse_proj_threaded,
         (((float *) mxGetData(plhs[0])) + ll * nb * na,
          nb, offset_r,((cfloat *) mxGetData(mx_image)) + ll * nx * ny,
          (cbyte *) mxGetData(mx_mask),
          nx, ny, dx, dy, offset_x, offset_y,
          nsubset, isubset, na, nthread))
  }
  LUT::free_LUTs();
  Ok
}


/*
 * jmh_sse_back_mex()
 */
static bool jmh_sse_back_mex(mxArray *plhs[],
                             Cmx mx_nthread,
                             Cmx mx_nb,
                             Cmx mx_offset_r,
                             Cmx mx_dr,
                             Cmx mx_dx,
                             Cmx mx_dy,
                             Cmx mx_offset_x,
                             Cmx mx_offset_y,
                             Cmx mx_mask,
                             Cmx mx_na,
                             Cmx mx_orbit_low,
                             Cmx mx_orbit_high,
                             Cmx mx_nsubset,
                             Cmx mx_isubset,
                             Cmx mx_proj)
{
  int ll, LL = 1;
  int nx, ny, nb, na, nthread;
  float dx, dy, dr;
  float offset_x, offset_y, offset_r;
  int nsubset, isubset;
  float orbit_low, orbit_high;

  /*
   * Check input types and sizes
   */
  Call(mxIsScalarInt32, (mx_nthread))
  nthread = mxGetInt(mx_nthread);
  
  Call(mxIsScalarInt32, (mx_nb))
  nb = mxGetInt(mx_nb);

  Call(mxIsScalarSingle, (mx_offset_r))
  offset_r = *((const float *)mxGetData(mx_offset_r));

  Call(mxIsScalarSingle, (mx_dr))
  dr = *((const float *) mxGetData(mx_dr));

  Call(mxIsScalarSingle, (mx_dx))
  dx = *((const float *)mxGetData(mx_dx));

  Call(mxIsScalarSingle, (mx_dy))
  dy = *((const float *)mxGetData(mx_dy));

  Call(mxIsScalarSingle, (mx_offset_x))
  offset_x = *((const float *)mxGetData(mx_offset_x));

  Call(mxIsScalarSingle, (mx_offset_y))
  offset_y = *((const float *)mxGetData(mx_offset_y));

  Call(mxIsScalarInt32, (mx_nsubset))
  nsubset = mxGetInt(mx_nsubset);

  Call(mxIsScalarInt32, (mx_isubset))
  isubset = mxGetInt(mx_isubset);

  if ((isubset < 0) || !(isubset < nsubset))
    Fail("subset out of range [0:nsubset-1]")

  Call(mxIsRealSingle, (mx_proj))

  Call(mxIsScalarInt32, (mx_na))
  na = mxGetInt(mx_na);

  if (nb != (int) mxGetM(mx_proj))
    Fail("nb mismatch")
      
  // trick: can't use mxGetN since it returns product
  if (na != mxGetDimensions(mx_proj)[1]*nsubset)
    Fail("na mismatch: na != nsubset*dim(proj,2)")

  if ((na / nsubset) % 4 != 0)
    Fail("(na/nsubset) % 4 != 0")    

  Call(mxIsUint8, (mx_mask))
  nx = mxGetM(mx_mask);
  ny = mxGetN(mx_mask);
    
  Call(mxIsScalarSingle, (mx_orbit_low))
  orbit_low = *((const float *) mxGetData(mx_orbit_low));

  Call(mxIsScalarSingle, (mx_orbit_high))
  orbit_high = *((const float *) mxGetData(mx_orbit_high));

#if 0
  if (1) {
    Note1("nthread=%d", nthread)
    Note2("isubset,nsubset=%d,%d", isubset,nsubset)
    Note1("na=%d", na)
    Note1("nb=%d", nb)
    Note1("offset_r=%g", offset_r)
    Note1("offset_x=%g", offset_x)
    Note1("offset_y=%g", offset_y)
    Note4("nx,ny=%d,%d dx,dy=%g,%g", nx, ny, dx, dy)
  }
#endif
      
  /*
   * output dimensions; allowing for multiple input sinograms
   */
  {
    int n_idim, *o_dims;
    cint *i_dims;
    n_idim = mxGetNumberOfDimensions(mx_proj);
    Call(i_dims = (cint*) mxGetDimensions, (mx_proj))
        
        //            Mem0(o_dims, n_idim)
    Alloc0(o_dims, int, n_idim, sizeof(*(o_dims)))

    Bcopy(i_dims, o_dims, n_idim)
    o_dims[0] = nx;
    o_dims[1] = ny;
    for (ll=2; ll < n_idim; ++ll)
      LL *= i_dims[ll];

    //    mexPrintf("%d\n",n_idim);
    /* allocate output */
    Call(plhs[0] = mxCreateNumericArray,
         (n_idim, (const mwSize*)o_dims, mxSINGLE_CLASS, mxREAL))
          Free0(o_dims)
  }

  LUT::initialize_LUTs(na, nsubset, orbit_low, orbit_high, dx, dr, dr); // dx, dr, sw
  
  for (ll=0; ll < LL; ++ll) {
    Call(jmh_sse_back_threaded,
         (((float *) mxGetData(plhs[0])) + ll * nx * ny,
          (cbyte *) mxGetData(mx_mask),
          nx, ny, dx, dy, offset_x, offset_y,
          ((cfloat *) mxGetData(mx_proj)) + ll * nb * na,
          nb, offset_r,
          nsubset, isubset, na,
          nthread))
  }
  LUT::free_LUTs();
  Ok
}


/*
 * jmh_sse_mex()
 */
bool jmh_sse_mex(cint nlhs, mxArray *plhs[], cint nrhs, const mxArray *prhs[]) {

  char *type;

  /* Check number of arguments */
        
  if (nrhs < 1 || nlhs < 1 || !mxIsChar(prhs[0])) {
                jmh_sse_mex_help();
                Fail(Usage)
  }

  Call(type = mxu_string, (prhs[0], "1st argument"))
    
  if Streq(type, "jmh,sse,proj") {
      if ((nrhs != 16) || (nlhs != 1)) {
          jmh_sse_mex_help();
          Fail(Usage)
      }
      Call(jmh_sse_proj_mex, (plhs, prhs[1], prhs[2], prhs[3],
                              prhs[4], prhs[5], prhs[6], prhs[7], prhs[8], prhs[9],
                              prhs[10],prhs[11],prhs[12],prhs[13],prhs[14],prhs[15]))
  }  
  else if Streq(type, "jmh,sse,back") {
      if ((nrhs != 16) || (nlhs != 1)){
        jmh_sse_mex_help();
        Fail(Usage)
      }
      Call(jmh_sse_back_mex, (plhs, prhs[1], prhs[2], prhs[3],
                              prhs[4], prhs[5], prhs[6], prhs[7], prhs[8], prhs[9],
                              prhs[10],prhs[11],prhs[12],prhs[13],prhs[14],prhs[15]))
  }
  else
    if Streq(type, "jmh,sse,matrix") {
        if ((nrhs != 16) || (nlhs != 2)) {
          jmh_sse_mex_help();
          Fail(Usage)
        }
        Call(jmh_sse_matrix, (plhs,prhs[1],prhs[2], prhs[3], prhs[4], prhs[5],
                              prhs[6], prhs[7], prhs[8], prhs[9], prhs[10],
                              prhs[11],prhs[12],prhs[13],prhs[14],prhs[15]))
    }
    else {
      jmh_sse_mex_help();
      Fail(Usage)
    }
  
  Call(mxu_string_free, (type))
    
  Ok
}

#if defined(Need_jmh_sse_mex_gateway)
/* gateway routine */
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	if (!nlhs && !nrhs) {
		jmh_sse_mex_help();
		return;
	}
	if (!jmh_sse_mex(nlhs, plhs, nrhs, prhs))
		mexErrMsgTxt("jmh_sse_mex()");
}
#endif

#undef Usage
