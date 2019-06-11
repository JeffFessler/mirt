// $Id: tomo_lut.cpp,v 1.13 2007/08/03 13:55:21 jmh Exp $
#include "tomo_lut.h"

#ifdef __ICC
 #pragma warning(disable:981) // Operands are evaluated in unspecified order
#endif

#define __MINMAX_DEFINED        // use STL's generic min and max templates
#define __USE_STL               

#define _NO_MULTI_THREADING 1

// STL include files - include STL files first!
#include <algorithm> 

#ifdef _WIN32
 #define NOMINMAX
#endif

// C++ standard include files
#include <math.h> //#include <sse_tomo_strip_2d/math.h>

#ifndef M_PI
 #define M_PI        3.14159265358979323846F
#endif

inline float deg2rad(float deg) {
  return deg*M_PI/180.0f;
}

#include <string.h>
//#include "threads.h"

//#include "mm_malloc.h"
#include "malloc.h"

#ifndef M_PI
 #define M_PI        3.14159265358979323846F
#endif



namespace LUT {
  // Doesn't initialize variables - I know
  volatile float* fcang = NULL;
  volatile float* fsang = NULL;
  volatile float* fdmax = NULL;
  volatile float* fdbrk = NULL;
  volatile float* ftaumax = NULL;
  volatile float* flang = NULL;
    
  // Default parameters - initialized using initialize_LUTs
  float dr = 0.0f;
  float sw = 0.0f;
  float dx = 0.0f;
  float norm = 0.0f;

  __m128 dx128 = _mm_setzero_ps();
  __m128 dr128 = _mm_setzero_ps();
  
  // Used by tomo_strip
  __m128 norm128 = _mm_setzero_ps();
  __m128 sw128 = _mm_setzero_ps();

#ifndef _NO_SSE

// Stripwidth: should be equal to dr or less
// TODO: introduce dy to allow different resolution in x and y   
bool initialize_LUTs(size_t nangles, size_t nsub, float orbit_low, float orbit_high,
                     float idx, float idr, float isw) {
 
  unsigned int i;
  
  dr = idr;
  sw = isw;
  dx = idx;
  
  norm = fabs(dx*dx)/sw; // fabs(dx*dy)/sw

  if ( (nangles % 4 != 0) || ((nangles / nsub) % 4 != 0) || (nangles % nsub != 0))
    return false;

  float* angles;
  
  // Angles  
  angles = (float*) _mm_malloc( (size_t) ((nangles) * sizeof(float)),16);

  size_t iangle;
  
  //std::cout << nangles << std::endl;
  
  // Arrange angles for ordered subsets {s1,s2,....}: [s1(0), s1(1),...., s1(#a/#subsets),s2(0),...] 
  for (i=0;i<nsub;i++) {
    for (unsigned int j=0 ; j<nangles/nsub ; j++) {
      iangle = i + j*nsub;
      #ifndef Mmex
        angles[nangles/nsub*i+j] = deg2rad(orbit_low + (orbit_high - orbit_low)*iangle/nangles);
      #else
        angles[nangles/nsub*i+j] = orbit_low + (orbit_high - orbit_low)*iangle/nangles;
      #endif
    }
  }

  dx128   = _mm_set_ps1(dx);
  dr128   = _mm_set_ps1(dr);
  sw128   = _mm_set_ps1(sw);
  norm128 = _mm_set_ps1(norm);
  
  // Allocate LUTs
  fcang   = (float*)_mm_malloc((size_t) ((nangles/4) * sizeof(__m128)),16);
  fsang   = (float*)_mm_malloc((size_t) ((nangles/4) * sizeof(__m128)),16);
  fdmax   = (float*)_mm_malloc((size_t) ((nangles/4) * sizeof(__m128)),16);
  fdbrk   = (float*)_mm_malloc((size_t) ((nangles/4) * sizeof(__m128)),16);
  ftaumax = (float*)_mm_malloc((size_t) ((nangles/4) * sizeof(__m128)),16);
  flang   = (float*)_mm_malloc((size_t) ((nangles/4) * sizeof(__m128)),16);
  
  __m128* cang = (__m128*)fcang;
  __m128* sang = (__m128*)fsang;
  __m128* dmax = (__m128*)fdmax;
  __m128* dbrk = (__m128*)fdbrk;
  __m128* taumax = (__m128*)ftaumax;
  __m128* lang = (__m128*)flang;

  // Initiate LUTs
  for (i=0; i< nangles/4 ; i++) {      

    // Used by tomo_strip_parallel_2d - make sure intrinsic sin/cos are not disabled -mfpmath=387 
    cang[i] = _mm_set_ps(cos(*(angles+3)),cos(*(angles+2)),cos(*(angles+1)),cos(*(angles)));
    sang[i] = _mm_set_ps(sin(*(angles+3)),sin(*(angles+2)),sin(*(angles+1)),sin(*(angles)));

    __m128 abscang = _mm_fabs_ps(cang[i]);
    __m128 abssang = _mm_fabs_ps(sang[i]);
    
    // Used by square_strip_int
    dmax[i] = _mm_mul_ps(_mm_add_ps(abscang,abssang),half128);
    dbrk[i] = _mm_mul_ps(_mm_fabs_ps(_mm_sub_ps(abscang,abssang)),half128); 
     
    // Used by tomo_strip_parallel_2d
    taumax[i] = _mm_div_ps(
                  _mm_add_ps(
                    _mm_mul_ps(
                      dmax[i],
                      dx128),
                    _mm_mul_ps(LUT::sw128,half128)),
                  dr128);

    __m128 denominator = _mm_max_ps(abscang,abssang);
    __m128 mask        = _mm_cmpeq_ps( denominator, _mm_setzero_ps() );
    denominator = _mm_or_ps(denominator, _mm_and_ps( mask, eps128 )); 

    // Used by square_strip_int
    lang[i] = _mm_rcp_ps(denominator);

    // Fix if dmax(dbrk==dmax) add 10*eps, not necessary
    dmax[i] = _mm_add_ps(_mm_and_ps(_mm_cmpeq_ps(dbrk[i],dmax[i]),_mm_mul_ps(eps128,_mm_set_ps1(10.0f))),dmax[i]);
    angles+=4;
  }
  angles -= nangles;

  // Free angles
  _mm_free((char*)angles);
  _mm_empty();
  _mm_sfence();
  return true;
}

#else
bool initialize_LUTs(size_t inangles, float iorbit_low, float iorbit_high,
		     float idx, float idr, float isw) {

  // TODO: arrange in struct
  orbit_low  = iorbit_low;
  orbit_high = iorbit_high;  
  dx = idx;
  dr = idr;
  sw = isw;
  nangles = inangles;
  
  norm = dx*dx/sw;

  float* angles;

  // Angles  
  angles = (float*) _mm_malloc( (size_t) ((nangles) * sizeof(float)));

  size_t iangle;
  for (unsigned int i=0;i<nangles;i++) {
    iangle = i/nsub + i%nsub;
    #ifndef Mmex
      angles[i] = deg2rad(orbit_low + (orbit_high - orbit_low)*iangle/nangles);
    #else
      angles[i] = orbit_low + (orbit_high - orbit_low)*iangle/nangles;
    #endif
  }

  // Allocate LUTs
  fcang   = (float*)_mm_malloc((size_t) ((nangles) * sizeof(__m128)));
  fsang   = (float*)_mm_malloc((size_t) ((nangles) * sizeof(__m128)));
  fdmax   = (float*)_mm_malloc((size_t) ((nangles) * sizeof(__m128)));
  fdbrk   = (float*)_mm_malloc((size_t) ((nangles) * sizeof(__m128)));
  ftaumax = (float*)_mm_malloc((size_t) ((nangles) * sizeof(__m128)));
  flang   = (float*)_mm_malloc((size_t) ((nangles) * sizeof(__m128)));
  
  // Initiate LUTs
  for (unsigned int i=0;i<nangles/4;i++) {      
    
    // Used by tomo_strip_parallel_2d
    cang[i] = cos(*angles);
    sang[i] = sin(*angles);

    float abscang = fabs(cang[i]);
    float abssang = fabs(sang[i]);
    
    // Used by square_strip_int
    dmax[i] = 0.5*(abscang+abssang);
    dbrk[i] = 0.5*(abscang-abssang);
    
    // Used by tomo_strip_parallel_2d
    taumax[i] = (dmax[i]*dx+0.5*sw)/dr;

    float denominator = max(abscang,abssang);
    // Remove branch
    denominator = (denominator == 0.0) ? eps : denominator; 

    // Used by square_strip_int
    lang[i] = 1.0/denominator;

    // Fix if dmax(dbrk==dmax) add 10*eps, not necessary
    // Remove branch
    dmax[i] = (dbrk == dmax) ? dmax[i] + 10*eps : dmax[i];
    angles++;
  }

  angles -= nangles;

  // Free angles
  _mm_free((char*)angles);
  _mm_empty();
  return true;

}
#endif
  
  
void free_LUTs() {
  if (fcang)
    _mm_free((char*)fcang);
  if (fsang)
    _mm_free((char*)fsang);  
  if (flang)
    _mm_free((char*)flang);  
  if (ftaumax)
    _mm_free((char*)ftaumax);  
  if (fdmax)
    _mm_free((char*)fdmax);  
  if (fdbrk)
    _mm_free((char*)fdbrk);
}

} // end namespace LUT

