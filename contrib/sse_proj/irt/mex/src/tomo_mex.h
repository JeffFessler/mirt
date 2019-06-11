// $Id: tomo_mex.h,v 1.23 2007/08/07 21:04:32 jmh Exp $
#ifndef TOMO_MEX_H
#define TOMO_MEX_H

//#define _LARGE_INDICES 1

#ifndef _NO_SSE
 #include "extintrin.h"
#else
 // Standard C++ headers
 #include <float.h>
#endif

template <class R> class sino_info {
 public:
  sino_info() : dx(R(1.0)), dy(R(1.0)),
                offset_x(R(0.0)), offset_y(R(0.0)), offset_r(R(0.0)) {};
  size_t nx;
  R      dx;
  R      offset_x;
  size_t ny;
  R      dy;
  R      offset_y;
  size_t nr;
  R      dr;
  R      sw;
  R      offset_r;
  size_t na;
  R      orbit_low;
  R      orbit_high;
};

template <class R> class alg_info {
 public:
  alg_info() : nthread(1), nsave(0), nsub(1) {};
  size_t niter;
  size_t nsub;
  size_t nthread;
  size_t nsave;
  R pixmax;
  R relax0;
  R relax1;
};

// Temporary struct for storing data
template <class R> class sino_data {
 public:
  sino_data() : data(NULL), randoms(NULL), eff(NULL), mask(NULL), init(NULL), ypi(NULL){};
  ~sino_data() {if(data) _mm_free(data); if(randoms) _mm_free(randoms); if(eff) _mm_free(eff);
  if(ypi) _mm_free(ypi); if(mask) _mm_free(mask); if(init) _mm_free(init);};
  R* data; // Precorrected noisy data 
  R* randoms;
  R* eff; // Measured using transmission
  R* ypi; // Uncorrected noisy data
  unsigned char* mask;
  R* init;
};

#ifndef _LARGE_INDICES
size_t tomo_strip_parallel_2d_mex_sub(const size_t nangles, const size_t nsub, const size_t isub,
                                  const size_t nb, const float dr, const float offset_r, const float sw,
                                  const size_t nx, const float dx, const float offset_x,
                                  const size_t ny, const float dy, const float offset_y,
                                  const unsigned char* mask, float* values, unsigned short* indices);
#else
size_t tomo_strip_parallel_2d_mex_sub(const size_t nangles, const size_t nsub, const size_t isub,                          
                              const size_t nb, const float dr, const float offset_r, const float sw,
                              const size_t nx, const float dx, const float offset_x,
                              const size_t ny, const float dy, const float offset_y,
                              const unsigned char* mask, float* values, unsigned long* indices);
#endif                              
                              
#ifndef _LARGE_INDICES
size_t tomo_strip_parallel_2d_mex(const size_t nangles,
                                  const size_t nb, const float dr, const float offset_r, const float sw,
                                  const size_t nx, const float dx, const float offset_x,
                                  const size_t ny, const float dy, const float offset_y,
                                  const unsigned char* mask, float* values, unsigned short* indices);
#else
size_t tomo_strip_parallel_2d_mex(const size_t nangles,
                                  const size_t nb, const float dr, const float offset_r, const float sw,
                                  const size_t nx, const float dx, const float offset_x,
                                  const size_t ny, const float dy, const float offset_y,
                                  const unsigned char* mask, float* values, unsigned long* indices);
#endif

size_t jmh_sse_proj_threaded_old(float* proj, const size_t nr, const float offset_r, const float *image,
                           const unsigned char* mask, const size_t nx, const size_t ny, const float dx,
                           const float dy, const float offset_x, const float offset_y, const size_t na, const size_t nthread);

size_t jmh_sse_proj_threaded(float* proj, const size_t nr, const float offset_r, const float *image,
                           const unsigned char* mask, const size_t nx, const size_t ny, const float dx,
                           const float dy, const float offset_x, const float offset_y, const size_t nsub, const size_t isub,
                           const size_t na, const size_t nthread);

size_t jmh_sse_proj_threaded_unrolled(float* proj, const size_t nr, const float offset_r,	const float	*image,
                           const unsigned char* mask, const size_t	nx, const size_t	ny, const float dx,
                           const float dy, const float offset_x, const float offset_y, const size_t na, const int nthread, int ithread = -1);


bool jmh_sse_back_threaded(float *image, const unsigned char* mask, const size_t nx, const size_t ny,
                           const float dx, const float dy, const float offset_x, const float offset_y,
                           const float *proj, const size_t nr, const float offset_r, const size_t nsub, const size_t isub,
                           const size_t na, const size_t nthread);

bool tomo_back_full_rank(float *image, const unsigned char* mask, const size_t nx, const size_t ny,
                         const float dx, const float dy, const float offset_x, const float offset_y,
                         const float *proj, const size_t nr, const float offset_r, const size_t nsub, const size_t isub,
                         const size_t na, const size_t nthread);

/*
#if defined(_WIN32) && defined(_MSC_VER)
 #define _THREAD_RETURN_TYPE unsigned __stdcall
#elif defined(__linux__)
 #define _THREAD_RETURN_TYPE void*
#endif

_THREAD_RETURN_TYPE proj_thread_function(void *arg);
_THREAD_RETURN_TYPE back_thread_function(void *arg);
*/

#endif
