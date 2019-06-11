/**********************************************************/
/* $Id: tomo_mex.cpp,v 1.55 2007/08/07 21:43:29 jmh Exp $ */
/* Copyright 2007-7-02, Jens Munk Hansen, NRU             */
/**********************************************************/

#include "tomo_lut.h"
#include "tomo_strip.h"
#include "tomo_mex.h"

#ifdef __ICC
 #pragma warning(disable:981)  // Operands are evaluated in unspecified order
 #pragma warning(disable:1599) // Hides variables
#endif

#define _SPARSE 1

#ifdef __RESTRICT
 #define RESTRICT restrict
#else
 #define RESTRICT __declspec(restrict)
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
#include <math.h>

#ifdef __linux__
 #include "mm_malloc.h"
#endif

#include <malloc.h>

// Original function tested againt Matlab source of J. A. Fessler

// Remove nthread, orbit_low, orbit_high

#ifndef _LARGE_INDICES
size_t tomo_strip_parallel_2d_mex(const size_t na, 
                              const size_t nb, const float dr, const float offset_r, const float sw,
                              const size_t nx, const float dx, const float offset_x,
                              const size_t ny, const float dy, const float offset_y,
                              const unsigned char* mask, float* values, unsigned short* indices) {                         
#else
size_t tomo_strip_parallel_2d_mex(const size_t na, 
                              const size_t nb, const float dr, const float offset_r, const float sw,
                              const size_t nx, const float dx, const float offset_x,
                              const size_t ny, const float dy, const float offset_y,
                              const unsigned char* mask, float* values, unsigned long* indices) {                         
#endif

  size_t M;
  float wx, wy, wb;
  float *xs, *ys;

  __m128 wb128, tau;
  __m128i ib_min;
  
  __m64  nb64;
  __m128i nb128;

  __m128 *cang, *sang, *taumax;
  __m128* m_Array;
  
  __m128* values128 = (__m128*) values;

#ifndef _LARGE_INDICES
  __m64* m_iArray;
  __m64* indices64 = (__m64*) indices;
#else
  __m128i* m_iArray;
  __m128i* indices64 = (__m128i*) indices;
#endif
  
  if ((na % 4 != 0))
    return 0;
  
  // Number of strips that intersects a pixel for a fixed angle
  M = (size_t) ceil((dx * sqrt(2.0f) + sw) / dr);
  
  // Pixel centers
  wx = dx*(nx-1.0f)/2.0f + offset_x;
  wy = dy*(ny-1.0f)/2.0f + offset_y;

  // R center
  wb = (nb-1.0f)/2.0f + offset_r;
  
  xs = (float*) _mm_malloc( (size_t) ((nx) * sizeof(float)),16);
  ys = (float*) _mm_malloc( (size_t) ((ny) * sizeof(float)),16);
  
  for (unsigned int i=0;i<nx;i++) xs[i] = dx * ((float) i)-wx;
  for (unsigned int i=0;i<ny;i++) ys[i] = dy * ((float) i)-wy;
  
#ifndef _LARGE_INDICES
  m_iArray = (__m64*) _mm_malloc(M*na/4 * sizeof(__m64),16); 
  __m64 max_indices = _mm_set1_pi16(((short)nb)*((short)na)-1);
#else
  m_iArray = (__m128i*) _mm_malloc(M*na/4 * sizeof(__m128i),16); 
  __m128i max_indices = _mm_set1_epi32(((int)nb)*((int)na)-1);
#endif

  nb64  = _mm_set1_pi16((short)nb);
  nb128  = _mm_set1_epi32((int)nb);
  wb128 = _mm_set_ps1(wb);

  // LUT references
  cang   = (__m128*) LUT::fcang;
  sang   = (__m128*) LUT::fsang;
  taumax = (__m128*) LUT::ftaumax;

  _mm_empty(); // Before _mm_malloc
  m_Array = (__m128*) _mm_malloc(M*na/4 * sizeof(__m128),16); 

  // TODO: restrict arguments!!! or register them explicitly  
  for (unsigned int j=0;j<ny;j++) {
    for (unsigned int i=0;i<nx;i++) {
      for (unsigned int k=0;k<na/4;k++) {        
        // Project pixel centers, (x * cos(a) + y * sin(a)) / dr + wb
        tau = _mm_add_ps(
                _mm_div_ps(
                  _mm_add_ps(                                 
                    _mm_mul_ps(cang[k],_mm_set_ps1((float) xs[i] )),
                    _mm_mul_ps(sang[k],_mm_set_ps1((float) ys[j] ))),
                  LUT::dr128),                                       // Divide by dr
                wb128);                                              // Add 'r-offset'

                        
        // Get first r-index (minus one)
        ib_min = _mm_cvttps_epi32(_mm_floor_ps(                                      
                   _mm_sub_ps(                                       // Subtract max r-diff
                     tau,
                     taumax[k])));

        // TODO: Use either small or large indices
        
        // Data index (#angle * nb + #r)
#ifndef _LARGE_INDICES    
        register __m64 dataoffset = 
                   _mm_mullo_pi16(
                     _mm_add_pi16(
                       zero2three64i,
                       _mm_set1_pi16(4*(short)k)),
                     nb64);
#else
        register __m128i dataoffset32 =
                        _mm_mul32_epi(
                          _mm_add_epi32(
                            zero2three128i,
                            _mm_set1_epi32(4*(int)k)),
                          nb128);        
#endif
        for (unsigned int l=0;l<M;l++)  { // Loop over r-values that intersects pixel
          // R-indices
          __m128i ib = _mm_add_epi32(ib_min,_mm_set1_epi32(l+1)); 

          // Calculate data indices(r,angle)
#ifndef _LARGE_INDICES
          __m64 indices16 = 
            _mm_add_pi16(
              dataoffset,
              _mm_movepi64_pi64(
                _mm_packs_epi32(
                  ib,
                  _mm_setzero_si128()))); // _mm_cvtps_pi16(fib)); // no truncation for __m64

          // Mask out indices that are too large to avoid access violation
          __m64 imask = _m_pand(_m_pcmpgtw(max_indices,indices16),_m_pcmpgtw(indices16,mones64));
          
          indices16 = _mm_sel_pi16(max_indices,indices16,imask);            
                   
#else
          register __m128i indices32 =
                      _mm_add_epi32(
                        dataoffset32,
                        ib);
                        
          // Mask out indices that are too large to avoid access violation
          __m128i imask = _mm_and_si128(_mm_cmpgt_epi32(max_indices,indices32),_mm_cmpgt_epi32(indices32,mones128));
          
          indices32 = _mm_sel_pi32(max_indices,indices32,imask);            
                        
#endif
          // Mask out some indices (0 <= r-index < nb) - TODO: use _mm_neg_pi
          __m128i rmask = _mm_and_si128(_mm_cmpgt_epi32(nb128,ib),_mm_cmpgt_epi32(ib,mones128));

          // R-values 
          __m128 res = _mm_mul_ps(
                         _mm_sub_ps(
                           _mm_cvtepi32_ps(ib),
                           tau),
                         LUT::dr128); 
          
          // Compute integrals
          __m128 myres = sse_square_strip_int(res, k);

          // Filter according to r-mask
#ifdef __linux__
          myres = _mm_sel_ps(_mm_setzero_ps(), myres, (__m128)rmask );
#else
          myres = _mm_sel_ps(_mm_setzero_ps(), myres, *(__m128*)&rmask );
#endif

          // Multiply by support mask
          myres = _mm_mul_ps(myres,_mm_cvtpu8_ps(_mm_set1_pi8(mask[i+j*nx])));
          // Dangerous, might need _mm_empty() after _mm_cvtpu8_ps

#ifndef NDEBUG          
          _mm_empty();
//          std::cout << "Indices: " << indices16 << std::endl;
          std::cout << "Values: " << myres << std::endl;
#endif   
            
          // Store results
          m_Array[M*k+l] = myres; // (a1,r1), (a1,r2), etc 
#ifndef _LARGE_INDICES
          m_iArray[M*k+l] = indices16;
#else
          m_iArray[M*k+l] = indices32;
#endif
        }
      }
      _mm_empty(); // Think this is necessary for memcpy
      
      // TODO: non-temporal cache writes using prefetching and SSE streaming
      memcpy((char*)(values128),(char*)m_Array,M*(na/4)*sizeof(__m128));
#ifndef _LARGE_INDICES
      memcpy((char*)(indices64),(char*)m_iArray,M*(na/4)*sizeof(__m64));
#else
      memcpy((char*)(indices64),(char*)m_iArray,M*(na/4)*sizeof(__m128i));
#endif
      // Expensive access to heap
      values128+=M*na/4;
      indices64+=M*na/4;
    }
  }

  // Move back output pointer
  values128 -= nx*ny*M*na/4;
  indices64 -= nx*ny*M*na/4;

  _mm_empty();

  // Flush temporary vectors
  _mm_free(xs);
  _mm_free(ys);
  
  // Flush temporary output
  _mm_free(m_Array);
  _mm_free(m_iArray);
  
  return nx*ny*M*na/4;
}

#ifndef _LARGE_INDICES
size_t tomo_strip_parallel_2d_mex_sub(const size_t na, const size_t nsub, const size_t isub,
                              const size_t nb, const float dr, const float offset_r, const float sw,
                              const size_t nx, const float dx, const float offset_x,
                              const size_t ny, const float dy, const float offset_y,
                              const unsigned char* mask, float* values, unsigned short* indices) {                         
#else
size_t tomo_strip_parallel_2d_mex_sub(const size_t na, const size_t nsub, const size_t isub,
                              const size_t nb, const float dr, const float offset_r, const float sw,
                              const size_t nx, const float dx, const float offset_x,
                              const size_t ny, const float dy, const float offset_y,
                              const unsigned char* mask, float* values, unsigned long* indices) {                         
#endif

  unsigned int M;
  float wx, wy, wb;
  float *xs, *ys;

  __m128 wb128, tau;
  __m128i ib_min;
  
  __m64  nb64;
  __m128i nb128;

  __m128 *cang, *sang, *taumax;
  __m128* m_Array;
  
  __m128* values128 = (__m128*) values;
  
#ifndef _LARGE_INDICES
  __m64* m_iArray;
  __m64* indices64 = (__m64*) indices;
#else
  __m128i* m_iArray;
  __m128i* indices64 = (__m128i*) indices;
#endif
    
  if ( (na % 4 != 0) || ((na / nsub) % 4 != 0) || (na % nsub != 0))
    return 0;

  // Number of strips that intersects a pixel for a fixed angle
  M = (unsigned int) ceil((dx * sqrt(2.0f) + sw) / dr);
  
  // Pixel centers
  wx = dx*(nx-1.0f)/2.0f + offset_x;
  wy = dy*(ny-1.0f)/2.0f + offset_y;

  // R center
  wb = (nb-1.0f)/2.0f + offset_r;
  // x- and y-values - if derived on the fly --> unresolved memory problems???
  
  xs = (float*) _mm_malloc( (size_t) ((nx) * sizeof(float)),16);
  ys = (float*) _mm_malloc( (size_t) ((ny) * sizeof(float)),16);
  
  for (unsigned int i=0;i<nx;i++) xs[i] = dx * ((float) i)-wx;
  for (unsigned int i=0;i<ny;i++) ys[i] = dy * ((float) i)-wy;

#ifndef _LARGE_INDICES
  m_iArray = (__m64*) _mm_malloc(M*na/4 * sizeof(__m64),16); 
  __m64 max_indices = _mm_set1_pi16(((short)nb)*((short)na)-1);
#else
  m_iArray = (__m128i*) _mm_malloc(M*na/4 * sizeof(__m128i),16); 
  __m128i max_indices = _mm_set1_epi32(((int)nb)*((int)na)-1);
#endif

  nb64  = _mm_set1_pi16((short)nb);
  nb128  = _mm_set1_epi32((int)nb);
  wb128 = _mm_set_ps1(wb);
  
  // LUT references
  cang   = (__m128*) LUT::fcang + isub*na/nsub/4;
  sang   = (__m128*) LUT::fsang + isub*na/nsub/4;
  taumax = (__m128*) LUT::ftaumax + isub*na/nsub/4;

  _mm_empty(); // Before _mm_malloc
  m_Array = (__m128*) _mm_malloc(M*na/4 * sizeof(__m128),16); 

  // TODO: restrict arguments!!! or register them explicitly  
  for (unsigned int j=0;j<ny;j++) {
    for (unsigned int i=0;i<nx;i++) {
      for (unsigned int k=0;k<na/nsub/4;k++) {        
        // Project pixel centers, (x * cos(a) + y * sin(a)) / dr + wb
        tau = _mm_add_ps(
                _mm_div_ps(
                  _mm_add_ps(                                 
                    _mm_mul_ps(cang[k],_mm_set_ps1((float) xs[i] )),
                    _mm_mul_ps(sang[k],_mm_set_ps1((float) ys[j] ))),
                  LUT::dr128),                                       // Divide by dr
                wb128);                                              // Add 'r-offset'

                        
        // Get first r-index (minus one)
        ib_min = _mm_cvttps_epi32(_mm_floor_ps(                                      
                   _mm_sub_ps(                                       // Subtract max r-diff
                     tau,
                     taumax[k])));

        // TODO: Use either small or large indices
        
        // Data index (#angle * nb + #r)
#ifndef _LARGE_INDICES
        register __m64 dataoffset = 
                  _mm_mullo_pi16(
                    _mm_add_pi16(
                      _mm_mullo_pi16(zero2three64i,_mm_set1_pi16((short)nsub)),
                      _mm_set1_pi16(((short)isub)+4*((short)nsub)*((short)k))), // LUT offsets[k]
                    nb64);        
        
#else
        register __m128i dataoffset32 =
                        _mm_mul32_epi(
                          _mm_add_epi32(
                            zero2three128i,
                            _mm_set1_epi32(4*(int)k)),
                          nb128);
                          
                                  
#endif
        for (unsigned int l=0;l<M;l++)  { // Loop over r-values that intersects pixel
          // R-indices
          __m128i ib = _mm_add_epi32(ib_min,_mm_set1_epi32(l+1)); 

          // Calculate data indices(r,angle)
#ifndef _LARGE_INDICES
          __m64 indices16 = 
            _mm_add_pi16(
              dataoffset,
              _mm_movepi64_pi64(
                _mm_packs_epi32(
                  ib,
                  _mm_setzero_si128()))); // _mm_cvtps_pi16(fib)); // no truncation for __m64

          // Mask out indices that are too large to avoid access violation
          __m64 imask = _m_pand(_m_pcmpgtw(max_indices,indices16),_m_pcmpgtw(indices16,mones64));
          
          indices16 = _mm_sel_pi16(max_indices,indices16,imask);            
                   
#else
          register __m128i indices32 =
                      _mm_add_epi32(
                        dataoffset32,
                        ib);
                           
          // Mask out indices that are too large to avoid access violation
          __m128i imask = _mm_and_si128(_mm_cmpgt_epi32(max_indices,indices32),_mm_cmpgt_epi32(indices32,mones128));
          
          indices32 = _mm_sel_pi32(max_indices,indices32,imask);  
#endif
          // Mask out some indices (0 <= r-index < nb) - TODO: use _mm_neg_pi
          __m128i rmask = _mm_and_si128(_mm_cmpgt_epi32(nb128,ib),_mm_cmpgt_epi32(ib,mones128));

          // R-values 
          __m128 res = _mm_mul_ps(
                         _mm_sub_ps(
                           _mm_cvtepi32_ps(ib),
                           tau),
                         LUT::dr128); 
          
          // Compute integrals
          __m128 myres = sse_square_strip_int(res, k);

          // Filter according to r-mask
#ifdef __linux__
          myres = _mm_sel_ps(_mm_setzero_ps(), myres, (__m128)rmask );
#else
          myres = _mm_sel_ps(_mm_setzero_ps(), myres, *(__m128*)&rmask );
#endif

          // Multiply by support mask
          myres = _mm_mul_ps(myres,_mm_cvtpu8_ps(_mm_set1_pi8(mask[i+j*nx])));
          // Dangerous, might need _mm_empty() after _mm_cvtpu8_ps

#ifndef NDEBUG          
#ifndef _LARGE_INDICES
          std::cout << "Indices: " << indices16 << std::endl;
#else
          std::cout << "Indices: " << indices32 << std::endl;
#endif

          std::cout << "Values: " << myres << std::endl;
#endif          
          
          // Store results
          m_Array[M*k+l] = myres; // (a1,r1), (a1,r2), etc 
#ifndef _LARGE_INDICES
          m_iArray[M*k+l] = indices16;
#else
          m_iArray[M*k+l] = indices32;
#endif
        }
      }

      _mm_empty(); // Think this is necessary for memcpy (ICC warning)
      
      // TODO: non-temporal cache writes using prefetching and SSE streaming
      memcpy((char*)(values128),(char*)m_Array,M*(na/4)*sizeof(__m128));
#ifndef _LARGE_INDICES
      memcpy((char*)(indices64),(char*)m_iArray,M*(na/4)*sizeof(__m64));
#else
      memcpy((char*)(indices64),(char*)m_iArray,M*(na/4)*sizeof(__m128i));
#endif

      // Expensive access to heap
      values128+=M*na/4;
      indices64+=M*na/4;
    }
  }

  // Move back output pointer
  values128 -= nx*ny*M*na/4;
  indices64 -= nx*ny*M*na/4;

  _mm_empty();
  _mm_sfence();

  // Flush temporary vectors
  _mm_free(xs);
  _mm_free(ys);
  
  // Flush temporary output
  _mm_free(m_Array);
  _mm_free(m_iArray);
  
  return nx*ny*M*na/4;
}

// Remove nthread, orbit_low, orbit_high

size_t jmh_sse_proj_threaded(float* proj, const size_t nb, const float offset_r, const float *image,
                           const unsigned char* mask, const size_t nx, const size_t ny, const float dx,
                           const float dy, const float offset_x, const float offset_y, const size_t na, const size_t nthread) {

  unsigned int M;
  float wx, wy, wb;

  __m128 *cang, *sang, *taumax;
  
  __m128 wb128, tau;
  __m128i ib_min;
  
  __m64  nb64;
  __m128i nb128;

  if ((na % 4 != 0))
    return 0;

  // Number of strips that intersects a pixel for a fixed angle
  M = (unsigned int) ceil((dx * sqrt(2.0f) + LUT::sw) / LUT::dr);
  
  // Pixel centers
  wx = dx*(nx-1.0f)/2.0f + offset_x;
  wy = dy*(ny-1.0f)/2.0f + offset_y;

  // R center
  wb = (nb-1.0f)/2.0f + offset_r;
  wb128 = _mm_set_ps1(wb);
  
  float *xs = (float*) _mm_malloc( (size_t) ((nx) * sizeof(float)),16);
  float *ys = (float*) _mm_malloc( (size_t) ((ny) * sizeof(float)),16);
  
  for (unsigned int i=0;i<nx;i++) xs[i] = dx * ((float) i)-wx;
  for (unsigned int i=0;i<ny;i++) ys[i] = dy * ((float) i)-wy;

  nb64  = _mm_set1_pi16((short)nb);
  nb128  = _mm_set1_epi32((int)nb);

  // LUT references
  cang   = (__m128*) LUT::fcang;
  sang   = (__m128*) LUT::fsang;
  taumax = (__m128*) LUT::ftaumax;

#ifndef _LARGE_INDICES
  register __m64 max_indices = _mm_set1_pi16(((short)nb)*((short)na)-1);
#else
  register __m128i max_indices = _mm_set1_epi32(((int)nb)*((int)na)-1);
#endif

  for (unsigned int k=0;k<na/4;k++) {
    
    // Data index (angle * #nb + r)
#ifndef _LARGE_INDICES
    // angle * #nb
    register __m64 dataoffset = 
      _mm_mullo_pi16(
        _mm_add_pi16(
          zero2three64i,
          _mm_set1_pi16(4*(short)k)),
        nb64);
#else
    register __m128i dataoffset32 =
      _mm_mul32_epi(
        _mm_add_epi32(
          zero2three128i,
          _mm_set1_epi32(4*(int)k)),
        nb128);        
#endif
    for (unsigned int j=0;j<ny;j++) {
      for (unsigned int i=0;i<nx;i++) {
        
        // Project pixel centers, (x * cos(a) + y * sin(a)) / dr + wb
        tau = _mm_add_ps(
                _mm_div_ps(
                  _mm_add_ps(                                 
                    _mm_mul_ps(cang[k],_mm_set_ps1((float) xs[i] )),
                    _mm_mul_ps(sang[k],_mm_set_ps1((float) ys[j] ))),
                  LUT::dr128),                                       // Divide by dr
                wb128);                                              // Add 'r-offset'
                        
        // Get first r-index (minus one)
        ib_min = _mm_cvttps_epi32(_mm_floor_ps(                                      
                   _mm_sub_ps(                                       // Subtract max r-diff
                     tau,
                     taumax[k])));

        
        for (unsigned int l=0;l<M;l++)  { // Loop over r-values that intersects pixel
          // R-indices
          __m128i ib = _mm_add_epi32(ib_min,_mm_set1_epi32(l+1));    // Add the 'minus' one

          // Calculate data indices(r,angle)
#ifndef _LARGE_INDICES
          __m64 indices16 = 
            _mm_add_pi16(
              dataoffset,
              _mm_movepi64_pi64(
                _mm_packs_epi32(
                  ib,
                  _mm_setzero_si128()))); // _mm_cvtps_pi16(fib)); // no truncation for __m64 

          // Mask out indices that are too large to avoid access violation
          __m64 imask = _m_pand(_m_pcmpgtw(max_indices,indices16),_m_pcmpgtw(indices16,mones64));
          
          indices16 = _mm_sel_pi16(max_indices,indices16,imask);
#else
          register __m128i indices32 =
                      _mm_add_epi32(
                        dataoffset32,
                        ib);

          // Mask out indices that are too large to avoid access violation
          __m128i imask = _mm_and_si128(_mm_cmpgt_epi32(max_indices,indices32),_mm_cmpgt_epi32(indices32,mones128));
          
          indices32 = _mm_sel_pi32(max_indices,indices32,imask);
#endif

          // Mask out some indices (0 <= r-index < nb) - TODO: use _mm_neg_pi
          __m128i rmask = _mm_and_si128(_mm_cmpgt_epi32(nb128,ib),_mm_cmpgt_epi32(ib,mones128));

          // R-values 
          __m128 res = _mm_mul_ps(
                         _mm_sub_ps(
                           _mm_cvtepi32_ps(ib),
                           tau),
                         LUT::dr128); 
          
          
          // Compute integrals
          __m128 sqint = sse_square_strip_int(res, k);

          // Filter according to r-mask
#ifdef __linux__
          sqint = _mm_sel_ps(_mm_setzero_ps(), sqint, (__m128)rmask );
#elif (defined(_MSC_VER) && defined(_WIN32))
          sqint = _mm_sel_ps(_mm_setzero_ps(), sqint, *(__m128*)&rmask );
#endif

          // Multiply by support mask -- no EMMS of _mm_cvtpu8_ps
          sqint = _mm_mul_ps(sqint,_mm_cvtpu8_ps(_mm_set1_pi8(mask[i+j*nx])));

          // Multiply by value
          sqint = _mm_mul_ps(sqint,_mm_set1_ps(image[i+j*nx]));
#ifdef __linux__
          float cproj[4] __attribute__ ((aligned (16)));
#elif (defined(_MSC_VER) && defined(_WIN32))         
          __declspec(align(16)) float cproj[4];
#endif

#ifndef _LARGE_INDICES
          // Load old result
          cproj[0] = proj[_mm_extract_pi16(indices16,0)];
          cproj[1] = proj[_mm_extract_pi16(indices16,1)];
          cproj[2] = proj[_mm_extract_pi16(indices16,2)];
          cproj[3] = proj[_mm_extract_pi16(indices16,3)];

#else
          // Load old result
          cproj[0] = proj[_mm_cvtsi128_si32(indices32)];
          indices32 = _mm_shuffle_epi32(indices32, _MM_SHUFFLE (0, 3, 2, 1));
          cproj[1] = proj[_mm_cvtsi128_si32(indices32)];
          indices32 = _mm_shuffle_epi32(indices32, _MM_SHUFFLE (0, 3, 2, 1));
          cproj[2] = proj[_mm_cvtsi128_si32(indices32)];
          indices32 = _mm_shuffle_epi32(indices32, _MM_SHUFFLE (0, 3, 2, 1));
          cproj[3] = proj[_mm_cvtsi128_si32(indices32)];
          indices32 = _mm_shuffle_epi32(indices32, _MM_SHUFFLE (0, 3, 2, 1));

#endif
          // Add new contrib to old result
          sqint = _mm_add_ps(_mm_load_ps((float*)cproj),sqint);

          // Store results - expensive
          _mm_store_ps((float*)cproj,sqint);

#ifndef _LARGE_INDICES
          proj[_mm_extract_pi16(indices16,0)] = cproj[0];
          proj[_mm_extract_pi16(indices16,1)] = cproj[1];
          proj[_mm_extract_pi16(indices16,2)] = cproj[2];
          proj[_mm_extract_pi16(indices16,3)] = cproj[3];
#else
          proj[_mm_cvtsi128_si32(indices32)] = cproj[0];
          indices32 = _mm_shuffle_epi32(indices32, _MM_SHUFFLE (0, 3, 2, 1));
          proj[_mm_cvtsi128_si32(indices32)] = cproj[1];
          indices32 = _mm_shuffle_epi32(indices32, _MM_SHUFFLE (0, 3, 2, 1));
          proj[_mm_cvtsi128_si32(indices32)] = cproj[2];
          indices32 = _mm_shuffle_epi32(indices32, _MM_SHUFFLE (0, 3, 2, 1));
          proj[_mm_cvtsi128_si32(indices32)] = cproj[3];
#endif
        } // end for (unsigned int l=0;l<M;l++)
      } // end for (unsigned int i=0;i<nx;i++)
    } // end for (unsigned int j=0;j<ny;j++)
  }
  
  _mm_empty();
  _mm_sfence();
  
  // Flush temporary vectors
  _mm_free(xs);
  _mm_free(ys);

  return true;
}

// Remove nthread, orbit_low, orbit_high

size_t jmh_sse_proj_threaded(float* proj, const size_t nb, const float offset_r, const float *image,
                           const unsigned char* mask, const size_t nx, const size_t ny, const float dx,
                           const float dy, const float offset_x, const float offset_y, const size_t nsub, const size_t isub,
                           const size_t na, const size_t nthread) {

  unsigned int M;
  float wx, wy, wb;

  __m128 *cang, *sang, *taumax;
  
  __m128 wb128, tau;
  __m128i ib_min;
  
  __m64  nb64;
  __m128i nb128;
  
  if ( (na % 4 != 0) || ((na / nsub) % 4 != 0) || (na % nsub != 0))
    return false;

  _mm_empty();
  
  // Number of strips that intersects a pixel for a fixed angle
  M = (unsigned int) ceil((dx * sqrt(2.0f) + LUT::sw) / LUT::dr);
  
  // Pixel centers
  wx = (nx-1.0f)/2.0f + offset_x;
  wy = (ny-1.0f)/2.0f + offset_y;

  // R center
  wb = (nb-1.0f)/2.0f + offset_r;
  
  float *xs = (float*) _mm_malloc( (size_t) ((nx) * sizeof(float)),16);
  float *ys = (float*) _mm_malloc( (size_t) ((ny) * sizeof(float)),16);
  
  for (unsigned int i=0;i<nx;i++) xs[i] = dx * (((float) i)-wx);
  for (unsigned int i=0;i<ny;i++) ys[i] = dy * (((float) i)-wy);

  wb128 = _mm_set_ps1(wb);
  nb64  = _mm_set1_pi16((short)nb);
  nb128  = _mm_set1_epi32((int)nb);

  // LUT references
  cang   = (__m128*) LUT::fcang + isub*na/nsub/4;
  sang   = (__m128*) LUT::fsang + isub*na/nsub/4;
  taumax = (__m128*) LUT::ftaumax + isub*na/nsub/4;

#ifndef _LARGE_INDICES
 #ifdef _SPARSE
  __m64 max_indices = _mm_set1_pi16(((short)nb)*((short)na/(short)nsub)-1);
 #else
  __m64 max_indices = _mm_set1_pi16(((short)nb)*((short)na)-1);
 #endif
#else
 #ifdef _SPARSE
  __m128i max_indices = _mm_set1_epi32(((int)nb)*((int)na/nsub)-1);
 #else
  __m128i max_indices = _mm_set1_epi32(((int)nb)*((int)na)-1);
 #endif
#endif

  for (unsigned int k = 0 ; k < na/4/nsub ; k++) {
    // Data index (#angle * nb + #r)
#ifndef _LARGE_INDICES

 #ifdef _SPARSE
    register __m64 dataoffset = 
                   _mm_mullo_pi16(
                     _mm_add_pi16(
                       zero2three64i,
                       _mm_set1_pi16(4*(short)k)),
                     nb64);
 #else
    register __m64 dataoffset = 
      _mm_mullo_pi16(
        _mm_add_pi16(
          _mm_mullo_pi16(zero2three64i,_mm_set1_pi16((short)nsub)),
          _mm_set1_pi16((short)isub+4*(short)nsub*(short)k)), // LUT offsets[k]
        nb64);
 #endif

#else

 #ifdef _SPARSE
    register __m128i dataoffset32 =
      _mm_mul32_epi(
        _mm_add_epi32(
          zero2three128i,
          _mm_set1_epi32(4*(int)k)),
        nb128);        
 #else
    register __m128i dataoffset32 =
      _mm_mul32_epi(
        _mm_add_epi32(
          _mm_mullo_epi32(zero2three128i,_mm_set1_epi32((int)nsub)),
          _mm_set1_epi32((int)isub+4*(int)nsub*(int)k)),
        nb128);        
 #endif

#endif

    // Reset projection - find another way (check if it is easier to re-allocate)
   
    for (unsigned int j=0;j<ny;j++) {
      for (unsigned int i=0;i<nx;i++) {
        
        // Project pixel centers, (x * cos(a) + y * sin(a)) / dr + wb
        tau = _mm_add_ps(
                _mm_div_ps(
                  _mm_add_ps(                                 
                    _mm_mul_ps(cang[k],_mm_set_ps1((float) xs[i] )),
                    _mm_mul_ps(sang[k],_mm_set_ps1((float) ys[j] ))),
                  LUT::dr128),                                       // Divide by dr
                wb128);                                              // Add 'r-offset'
                        
        // Get first r-index (minus one)
        ib_min = _mm_cvttps_epi32(_mm_floor_ps(                                      
                   _mm_sub_ps(                                       // Subtract max r-diff
                     tau,
                     taumax[k])));

        
        for (unsigned int l=0;l<M;l++)  { // Loop over r-values that intersects pixel
          // R-indices
          __m128i ib = _mm_add_epi32(ib_min,_mm_set1_epi32(l+1)); 

          // Calculate data indices(r,angle)
#ifndef _LARGE_INDICES
          __m64 indices16 = 
            _mm_add_pi16(
              dataoffset,
              _mm_movepi64_pi64(
                _mm_packs_epi32(
                  ib,
                  _mm_setzero_si128()))); // _mm_cvtps_pi16(fib)); // no truncation for __m64 

          // Mask out indices that are too large to avoid access violation
          __m64 imask = _m_pand(_m_pcmpgtw(max_indices,indices16),_m_pcmpgtw(indices16,mones64));
          indices16 = _mm_sel_pi16(max_indices,indices16,imask);
#else
          register __m128i indices32 =
                      _mm_add_epi32(
                        dataoffset32,
                        ib);

          // Mask out indices that are too large to avoid access violation
          __m128i imask = _mm_and_si128(_mm_cmpgt_epi32(max_indices,indices32),_mm_cmpgt_epi32(indices32,mones128));
          
          indices32 = _mm_sel_pi32(max_indices,indices32,imask);            
          
#endif

          // Mask out some indices (0 <= r-index < nb) - TODO: use _mm_neg_pi
          __m128i rmask = _mm_and_si128(_mm_cmpgt_epi32(nb128,ib),_mm_cmpgt_epi32(ib,mones128));

          // R-values 
          __m128 res = _mm_mul_ps(
                         _mm_sub_ps(
                           _mm_cvtepi32_ps(ib),
                           tau),
                         LUT::dr128); 
          
          
          // Compute integrals
          __m128 sqint = sse_square_strip_int(res, isub*na/nsub/4+k);

          // Filter according to r-mask
#ifdef __linux__
          sqint = _mm_sel_ps(_mm_setzero_ps(), sqint, (__m128)rmask );
#elif (defined(_MSC_VER) && defined(_WIN32))
          sqint = _mm_sel_ps(_mm_setzero_ps(), sqint, *(__m128*)&rmask );
#endif

          // Multiply by support mask -- no EMMS of _mm_cvtpu8_ps
          sqint = _mm_mul_ps(sqint,_mm_cvtpu8_ps(_mm_set1_pi8(mask[i+j*nx])));

          // Multiply by value
          sqint = _mm_mul_ps(sqint,_mm_set1_ps(image[i+j*nx]));
#ifdef __linux__
          float cproj[4] __attribute__ ((aligned (16)));
#elif (defined(_MSC_VER) && defined(_WIN32))         
          __declspec(align(16)) float cproj[4];
#endif

#ifndef _LARGE_INDICES
          // Load old result
          cproj[0] = proj[_mm_extract_pi16(indices16,0)];
          cproj[1] = proj[_mm_extract_pi16(indices16,1)];
          cproj[2] = proj[_mm_extract_pi16(indices16,2)];
          cproj[3] = proj[_mm_extract_pi16(indices16,3)];
#else
          // Load old result
          cproj[0] = proj[_mm_cvtsi128_si32(indices32)];
          indices32 = _mm_shuffle_epi32(indices32, _MM_SHUFFLE (0, 3, 2, 1));
          cproj[1] = proj[_mm_cvtsi128_si32(indices32)];
          indices32 = _mm_shuffle_epi32(indices32, _MM_SHUFFLE (0, 3, 2, 1));
          cproj[2] = proj[_mm_cvtsi128_si32(indices32)];
          indices32 = _mm_shuffle_epi32(indices32, _MM_SHUFFLE (0, 3, 2, 1));
          cproj[3] = proj[_mm_cvtsi128_si32(indices32)];
          indices32 = _mm_shuffle_epi32(indices32, _MM_SHUFFLE (0, 3, 2, 1));
#endif
          // Add new contrib to old result
          sqint = _mm_add_ps(_mm_load_ps((float*)cproj),sqint);

          // Store results - expensive
          _mm_store_ps((float*)cproj,sqint);

#ifndef _LARGE_INDICES
          proj[_mm_extract_pi16(indices16,0)] = cproj[0];
          proj[_mm_extract_pi16(indices16,1)] = cproj[1];
          proj[_mm_extract_pi16(indices16,2)] = cproj[2];
          proj[_mm_extract_pi16(indices16,3)] = cproj[3];
#else
          proj[_mm_cvtsi128_si32(indices32)] = cproj[0];
          indices32 = _mm_shuffle_epi32(indices32, _MM_SHUFFLE (0, 3, 2, 1));
          proj[_mm_cvtsi128_si32(indices32)] = cproj[1];
          indices32 = _mm_shuffle_epi32(indices32, _MM_SHUFFLE (0, 3, 2, 1));
          proj[_mm_cvtsi128_si32(indices32)] = cproj[2];
          indices32 = _mm_shuffle_epi32(indices32, _MM_SHUFFLE (0, 3, 2, 1));
          proj[_mm_cvtsi128_si32(indices32)] = cproj[3];

#endif
        }
      }
    }
  }
  _mm_empty();
  _mm_sfence();

  // Flush temporary vectors
  _mm_free(xs);
  _mm_free(ys);
  return true;
}

// Remove nthread, orbit_low, orbit_high

// Experiment for un-rolling an inner loop in hope that the compiler will keep FPU busy at all times
size_t jmh_sse_proj_threaded_unrolled(float* proj, const size_t nb, const float offset_r, const float *image,
                           const unsigned char* mask, const size_t nx, const size_t ny, const float dx,
                           const float dy, const float offset_x, const float offset_y, const size_t na, const int nthread, int ithread) {

  unsigned int M;
  float wx, wy, wb;

  __m128 *cang, *sang, *taumax;
  
  __m128 wb128, tau[4];
  __m128i ib_min[4];
  
  __m64  nb64;
  __m128i nb128;

#ifndef _NO_MULTI_THREADING 
  
  if ((ithread < 0) && (nthread == 1)) {
    ithread = 0;
  }

  if ((nthread != 1) && (ithread < 0)){

    unsigned int threadID;

    #ifdef _WIN32
      HANDLE threads[MAXTHREADS];
    #else
      pthread_t threads[MAXTHREADS];
    #endif

    jmh_sse_proj_params TP[2];
    TP[0].proj = proj;
    TP[0].nb = nb;
    TP[0].offset_r = offset_r;
    TP[0].image = *const_cast<float**>(&image);
    TP[0].mask = *const_cast<unsigned char**>(&mask);
    TP[0].nx = nx;
    TP[0].ny = ny;
    TP[0].dx = dx;
    TP[0].dy = dy;
    TP[0].offset_x = offset_x;
    TP[0].offset_y = offset_y;
    TP[0].angles = *const_cast<float**>(&angles);
    TP[0].na = na;
    TP[0].nthread = nthread;
    TP[0].ithread = 0;

    TP[1] = TP[0];
    TP[1].ithread = 1;

    // Start threads
    for (int thr=0; thr < nthread; thr++) { 
      #ifdef _WIN32
        threads[thr] = (HANDLE)_beginthreadex( NULL, 0, &proj_thread_function, &TP[thr], 0, &threadID );
      #else
        pthread_create( &(mythread[thr]), NULL, proj_thread_function, (void*)&(thread_arg[thr])))
      #endif
    }

    // Wait for threads to complete
    for (int thr = 0; thr < nthread; thr++) {
      Wait_thread(threads[thr]);
    } 
    return true;
  }
  else {

#endif
    // No multi-threading

  // Number of strips that intersects a pixel for a fixed angle
  M = (unsigned int) ceil((dx * sqrt(2.0f) + LUT::sw) / LUT::dr);
  
  // Pixel centers
  wx = dx*(nx-1.0f)/2.0f + offset_x;
  wy = dy*(ny-1.0f)/2.0f + offset_y;

  // R center
  wb = (nb-1.0f)/2.0f + offset_r;
  wb128 = _mm_set_ps1(wb);
  
  // x- and y-values - if derived on the fly --> unresolved memory problems (EMMS)
  
#ifdef __linux__
  float *xs = (float*) _mm_malloc( (size_t) ((nx) * sizeof(float)),16);
  float *ys = (float*) _mm_malloc( (size_t) ((ny) * sizeof(float)),16);
#elif (defined(_MSC_VER) && defined(_WIN32))         
  __declspec(align(16)) float *xs = (float*) _mm_malloc( (size_t) ((nx) * sizeof(float)),16);
  __declspec(align(16)) float *ys = (float*) _mm_malloc( (size_t) ((ny) * sizeof(float)),16);
#endif

  
  for (unsigned int i=0;i<nx;i++) xs[i] = dx * ((float) i)-wx;
  for (unsigned int i=0;i<ny;i++) ys[i] = dy * ((float) i)-wy;
  //  for (unsigned int i=ny;i>0;i--) ys[i] = dy * ((float) (i-1))-wy;

  nb64  = _mm_set1_pi16((short)nb);
  nb128  = _mm_set1_epi32((int)nb);

  // LUT references
  cang   = (__m128*) LUT::fcang;
  sang   = (__m128*) LUT::fsang;
  taumax = (__m128*) LUT::ftaumax;

  __m128 sqint[4];

  // Create temporary image and apply mask is way too expensive

#ifndef _NO_MULTI_THREADING  
  for (int k = ithread * ((int)na/4/nthread) ; k <  ithread * ((int)na/4/nthread) + (int)na/4/nthread ; k++) {
#else
  for (unsigned int k = 0 ; k <  na/4 ; k++) {
#endif
    // Data index (#angle * nb + #r)
#ifndef _LARGE_INDICES
    register __m64 dataoffset = 
      _mm_mullo_pi16(
        _mm_add_pi16(
          zero2three64i,
          _mm_set1_pi16(4*(short)k)),
        nb64);
    __m64 max_indices = _mm_set1_pi16((short)nb*(short)na-1);
#else
    register __m128i dataoffset32 =
      _mm_mul32_epi(
        _mm_add_epi32(
          zero2three128i,
          _mm_set1_epi32(4*(int)k)),
        nb128);        
#endif
    // TODO: Prefetch data (smaller portions)

    _mm_prefetch((char*)(&proj[nx*ny*k]+nb*4*sizeof(float)),_MM_HINT_NTA);
    for (unsigned int j=0;j<ny;j++) {
      for (unsigned int i=0;i<nx;i+=4) { // Un-roll'ed by four (experiment) 
        _mm_prefetch((char*)(&image[i+j*nx]+4*sizeof(float)),_MM_HINT_NTA);

        // Temporary image values
#ifdef __linux__
        float timage[4] __attribute__ ((aligned (16)));
#elif (defined(_MSC_VER) && defined(_WIN32))
        __declspec(align(16)) float timage[4];
#endif

        // Load mask - prefetching is carried out by compiler
        __m128 tmask = _mm_cvtpu8_ps(
                         _mm_set_pi8(
                           0,0,0,0,
                           mask[i+j*nx],mask[1+i+j*nx],mask[2+i+j*nx],mask[3+i+j*nx]));

        // Apply mask
        __m128 ttimage = _mm_mul_ps(tmask,_mm_load_ps(&image[i+j*nx]));

        // Store masked image values
        _mm_store_ps((float*)timage,ttimage);

        // Project pixel centers, (x * cos(a) + y * sin(a)) / dr + wb
        tau[0] = _mm_add_ps(
                   _mm_div_ps(
                     _mm_add_ps(                                 
                       _mm_mul_ps(cang[k],_mm_set_ps1((float) xs[i] )),
                       _mm_mul_ps(sang[k],_mm_set_ps1((float) ys[j] ))),
                     LUT::dr128),                                       // Divide by dr
                   wb128);                                              // Add 'r-offset'

        tau[1] = _mm_add_ps(_mm_div_ps(_mm_add_ps(_mm_mul_ps(cang[k],_mm_set_ps1((float) xs[i+1] )),
                                                  _mm_mul_ps(sang[k],_mm_set_ps1((float) ys[j] ))),
                                       LUT::dr128), wb128);
        tau[2] = _mm_add_ps(_mm_div_ps(_mm_add_ps(_mm_mul_ps(cang[k],_mm_set_ps1((float) xs[i+2] )),
                                                  _mm_mul_ps(sang[k],_mm_set_ps1((float) ys[j] ))),
                                       LUT::dr128), wb128);
        tau[3] = _mm_add_ps(_mm_div_ps(_mm_add_ps(_mm_mul_ps(cang[k],_mm_set_ps1((float) xs[i+3] )),
                                                  _mm_mul_ps(sang[k],_mm_set_ps1((float) ys[j] ))),
                                       LUT::dr128), wb128);

        // Get first r-index (minus one)
        ib_min[0] = _mm_cvttps_epi32(_mm_floor_ps(                                      
                                       _mm_sub_ps(              // Subtract max r-diff
                                         tau[0],
                                         taumax[k])));

        ib_min[1] = _mm_cvttps_epi32(_mm_floor_ps(_mm_sub_ps(tau[1],taumax[k])));
        ib_min[2] = _mm_cvttps_epi32(_mm_floor_ps(_mm_sub_ps(tau[2],taumax[k])));
        ib_min[3] = _mm_cvttps_epi32(_mm_floor_ps(_mm_sub_ps(tau[3],taumax[k])));

        for (unsigned int l=0;l<M;l++)  { // Loop over r-values that intersects pixel
          // R-indices
          __m128i ib[4];
          // R-values
          __m128 res[4];
          // R-mask (0 =< r < nb)
          __m128i rmask[4];
          // Data indices(r,angle)
          __m64 indices16[4];

          // Current values
#ifdef __linux__
          float cproj[4] __attribute__ ((aligned(16)));
#elif (defined(_MSC_VER) && defined(_WIN32))
          __declspec(align(16)) float cproj[4];
#endif

          // R-indices
          ib[0] = _mm_add_epi32(ib_min[0],_mm_set1_epi32(l+1)); 
          ib[1] = _mm_add_epi32(ib_min[1],_mm_set1_epi32(l+1)); 
          ib[2] = _mm_add_epi32(ib_min[2],_mm_set1_epi32(l+1)); 
          ib[3] = _mm_add_epi32(ib_min[3],_mm_set1_epi32(l+1)); 

          // Calculate data indices(r,angle)
#ifndef _LARGE_INDICES
          indices16[0] = 
            _mm_add_pi16(
              dataoffset,
              _mm_movepi64_pi64(
                _mm_packs_epi32(
                  ib[0],
                  _mm_setzero_si128()))); // _mm_cvtps_pi16(fib)); // no truncation for __m64 

          // Mask out indices that are too large to avoid AV
          __m64 imask = _m_pcmpgtw(max_indices,indices16[0]);
          indices16[0] = _mm_sel_pi16(max_indices,indices16[0],imask);
#else
          register __m128i indices32 =
                      _mm_add_epi32(
                        dataoffset32,
                        ib[0]);
#endif
          // Calculate data indices(r,angle)
#ifndef _LARGE_INDICES
          indices16[1] = _mm_add_pi16(dataoffset,_mm_movepi64_pi64(_mm_packs_epi32(ib[1],_mm_setzero_si128())));
          imask = _m_pcmpgtw(max_indices,indices16[1]);
          indices16[1] = _mm_sel_pi16(max_indices,indices16[1],imask);
          indices16[2] = _mm_add_pi16(dataoffset,_mm_movepi64_pi64(_mm_packs_epi32(ib[2],_mm_setzero_si128())));
          imask = _m_pcmpgtw(max_indices,indices16[2]);
          indices16[2] = _mm_sel_pi16(max_indices,indices16[2],imask);
          indices16[3] = _mm_add_pi16(dataoffset,_mm_movepi64_pi64(_mm_packs_epi32(ib[3],_mm_setzero_si128())));
          imask = _m_pcmpgtw(max_indices,indices16[3]);
          indices16[3] = _mm_sel_pi16(max_indices,indices16[3],imask);
#endif

          // R-mask (0 =< r < nb)
          rmask[0] = _mm_and_si128(_mm_cmpgt_epi32(nb128,ib[0]),_mm_cmpgt_epi32(ib[0],mones128));

          // R-values 
          res[0] = _mm_mul_ps(_mm_sub_ps(_mm_cvtepi32_ps(ib[0]), tau[0]), LUT::dr128); 
          rmask[1] = _mm_and_si128(_mm_cmpgt_epi32(nb128,ib[1]), _mm_cmpgt_epi32(ib[1],mones128));
          res[1] = _mm_mul_ps(_mm_sub_ps(_mm_cvtepi32_ps(ib[1]), tau[1]), LUT::dr128); 
          rmask[2] = _mm_and_si128(_mm_cmpgt_epi32(nb128,ib[2]), _mm_cmpgt_epi32(ib[2],mones128));
          res[2] = _mm_mul_ps(_mm_sub_ps(_mm_cvtepi32_ps(ib[2]), tau[2]), LUT::dr128); 
          rmask[3] = _mm_and_si128(_mm_cmpgt_epi32(nb128,ib[3]), _mm_cmpgt_epi32(ib[3],mones128));
          res[3] = _mm_mul_ps(_mm_sub_ps(_mm_cvtepi32_ps(ib[3]), tau[3]), LUT::dr128); 
          
          // Compute integrals
          sqint[0] = sse_square_strip_int(res[0], k);

          // Filter according to r-mask
#ifdef __linux__
          sqint[0] = _mm_sel_ps(_mm_setzero_ps(), sqint[0], (__m128)rmask[0] );
          sqint[1] = sse_square_strip_int(res[1], k);
          sqint[1] = _mm_sel_ps(_mm_setzero_ps(), sqint[0], (__m128)rmask[1] );
          sqint[2] = sse_square_strip_int(res[2], k);
          sqint[2] = _mm_sel_ps(_mm_setzero_ps(), sqint[0], (__m128)rmask[2] );
          sqint[3] = sse_square_strip_int(res[3], k);
          sqint[3] = _mm_sel_ps(_mm_setzero_ps(), sqint[0], (__m128)rmask[3] );
#else
          // Mask out some indices (0 <= r-index < nb)
          sqint[0] = _mm_sel_ps(_mm_setzero_ps(), sqint[0], *(__m128*)&rmask[0]);
          sqint[1] = sse_square_strip_int(res[1], k);
          sqint[1] = _mm_sel_ps(_mm_setzero_ps(), sqint[1], *(__m128*)&rmask[1]);
          sqint[2] = sse_square_strip_int(res[2], k);
          sqint[2] = _mm_sel_ps(_mm_setzero_ps(), sqint[2], *(__m128*)&rmask[2]);
          sqint[3] = sse_square_strip_int(res[3], k);
          sqint[3] = _mm_sel_ps(_mm_setzero_ps(), sqint[3], *(__m128*)&rmask[3]);
#endif

          // Multiply by value... use two-vectors to one matrix
          sqint[0] = _mm_mul_ps(sqint[0],_mm_set1_ps(timage[0]));
          sqint[1] = _mm_mul_ps(sqint[1],_mm_set1_ps(timage[1]));
          sqint[2] = _mm_mul_ps(sqint[2],_mm_set1_ps(timage[2]));
          sqint[3] = _mm_mul_ps(sqint[3],_mm_set1_ps(timage[3]));

          cproj[0] = proj[_mm_extract_pi16(indices16[0], 0)];
          cproj[1] = proj[_mm_extract_pi16(indices16[0], 1)];
          cproj[2] = proj[_mm_extract_pi16(indices16[0], 2)];
          cproj[3] = proj[_mm_extract_pi16(indices16[0], 3)];

          // Add to old result
          sqint[0] = _mm_add_ps(_mm_load_ps((float*)cproj),sqint[0]);

          // Store results - expensive 
          _mm_store_ps((float*)cproj,sqint[0]);

          proj[_mm_extract_pi16(indices16[0],0)] = cproj[0];
          cproj[0] = proj[_mm_extract_pi16(indices16[1], 0)];
          proj[_mm_extract_pi16(indices16[0],1)] = cproj[1];
          cproj[1] = proj[_mm_extract_pi16(indices16[1], 1)];
          proj[_mm_extract_pi16(indices16[0],2)] = cproj[2];
          cproj[2] = proj[_mm_extract_pi16(indices16[1], 2)];
          proj[_mm_extract_pi16(indices16[0],3)] = cproj[3];
          cproj[3] = proj[_mm_extract_pi16(indices16[1], 3)];

          // Add to old result
          sqint[1] = _mm_add_ps(_mm_load_ps((float*)cproj),sqint[1]);

          // Store results - expensive 
         _mm_store_ps((float*)cproj,sqint[1]);

         proj[_mm_extract_pi16(indices16[1],0)] = cproj[0];
         cproj[0] = proj[_mm_extract_pi16(indices16[2], 0)];
         proj[_mm_extract_pi16(indices16[1],1)] = cproj[1];
         cproj[1] = proj[_mm_extract_pi16(indices16[2], 1)];
         proj[_mm_extract_pi16(indices16[1],2)] = cproj[2];
         cproj[2] = proj[_mm_extract_pi16(indices16[2], 2)];
         proj[_mm_extract_pi16(indices16[1],3)] = cproj[3];
         cproj[3] = proj[_mm_extract_pi16(indices16[2], 3)];

          // add to old result
          sqint[2] = _mm_add_ps(_mm_load_ps((float*)cproj),sqint[2]);

          // Store results - expensive 
         _mm_store_ps((float*)cproj,sqint[2]);

         proj[_mm_extract_pi16(indices16[2],0)] = cproj[0];
         cproj[0] = proj[_mm_extract_pi16(indices16[3], 0)];
         proj[_mm_extract_pi16(indices16[2],1)] = cproj[1];
         cproj[1] = proj[_mm_extract_pi16(indices16[3], 1)];
         proj[_mm_extract_pi16(indices16[2],2)] = cproj[2];
         cproj[2] = proj[_mm_extract_pi16(indices16[3], 2)];
         proj[_mm_extract_pi16(indices16[2],3)] = cproj[3];
         cproj[3] = proj[_mm_extract_pi16(indices16[3], 3)];

          // add to old result
          sqint[3] = _mm_add_ps(_mm_load_ps((float*)cproj),sqint[3]);

          // Store results - expensive 
         _mm_store_ps((float*)cproj,sqint[3]);

         proj[_mm_extract_pi16(indices16[3],0)] = cproj[0];
         proj[_mm_extract_pi16(indices16[3],1)] = cproj[1];
         proj[_mm_extract_pi16(indices16[3],2)] = cproj[2];
         proj[_mm_extract_pi16(indices16[3],3)] = cproj[3];

        }
      }
    }
  }

  // Flush temporary vectors
  _mm_free(xs);
  _mm_free(ys);

  _mm_empty();
  _mm_sfence();            // Guarantees that every preceding store is globally visible
                           // before any subsequent store. 
  return true;
#ifndef _NO_MULTI_THREADING
  }
#endif
}

// Remove nthread, orbit_low, orbit_high

bool jmh_sse_back_threaded(float *image, const unsigned char* mask, const size_t nx, const size_t ny,
                           const float dx, const float dy, const float offset_x, const float offset_y,
                           const float *proj, const size_t nr, const float offset_r, const size_t nsub, const size_t isub,
                           const size_t na, const size_t nthread) {
  
  unsigned int M;
  float wx, wy, wb;

  __m128 *cang, *sang, *taumax;
  
  __m128 wb128, tau;
  __m128i ib_min;
  
  __m64  nb64;
  __m128i nb128;
  
  if ((na % 4 != 0) || ((na / nsub) % 4 != 0))
    return false;

  // Number of strips that intersects a pixel for a fixed angle
  M = (unsigned int) ceil((dx * sqrt(2.0f) + LUT::sw) / LUT::dr);

  // Pixel centers
  wx = (nx-1.0f)/2.0f + offset_x;
  wy = (ny-1.0f)/2.0f + offset_y;

  // R center
  wb = (nr-1.0f)/2.0f + offset_r;

  _mm_empty(); // Before floating point operations
  
  float *xs = (float*) _mm_malloc( (size_t) ((nx) * sizeof(float)),16);
  float *ys = (float*) _mm_malloc( (size_t) ((ny) * sizeof(float)),16);
  
  for (unsigned int i=0;i<nx;i++) xs[i] = dx * ((float) i-wx);
  for (unsigned int i=0;i<ny;i++) ys[i] = dy * ((float) i-wy);
  
#ifndef _LARGE_INDICES
#ifdef _SPARSE
  __m64 max_indices = _mm_set1_pi16(((short)nr)*((short)na/(short)nsub)-1);
#else
  __m64 max_indices = _mm_set1_pi16(((short)nr)*((short)na)-1);
#endif
#else
#ifdef _SPARSE
  __m128i max_indices = _mm_set1_epi32(((int)nr)*((int)na/nsub)-1);
#else
  __m128i max_indices = _mm_set1_epi32(((int)nr)*((int)na)-1);
#endif
#endif

  wb128 = _mm_set_ps1(wb);
  nb64  = _mm_set1_pi16((short)nr);
  nb128  = _mm_set1_epi32((int)nr);

  // LUT references
  cang   = (__m128*) LUT::fcang;
  sang   = (__m128*) LUT::fsang;
  taumax = (__m128*) LUT::ftaumax;

  // If not values are NaNs
#ifdef Mmex
  for (unsigned int i=0 ; i < nx*ny ; i++) {
    image[i] = 0.0f;
  }
#endif
  
  
  for (unsigned int k = 0 ; k < na/4/nsub ; k++) {
    // Data index (#angle * nb + #r)
#ifndef _LARGE_INDICES
#ifdef _SPARSE
    register __m64 dataoffset = 
                   _mm_mullo_pi16(
                     _mm_add_pi16(
                       zero2three64i,
                       _mm_set1_pi16(4*(short)k)),
                     nb64);
#else
    register __m64 dataoffset = 
      _mm_mullo_pi16(
        _mm_add_pi16(
          _mm_mullo_pi16(zero2three64i,_mm_set1_pi16((short)nsub)),
          _mm_set1_pi16((short)isub+4*(short)nsub*(short)k)), // LUT offsets[k]
        nb64);
#endif
#else
#ifdef _SPARSE
    register __m128i dataoffset32 =
      _mm_mul32_epi(
        _mm_add_epi32(
          zero2three128i,
          _mm_set1_epi32(4*(int)k)),
        nb128);        
#else
    register __m128i dataoffset32 = 
      _mm_mul32_epi(
        _mm_add_epi32(
          _mm_mullo_epi32(zero2three128i,_mm_set1_epi32((int)nsub)),
          _mm_set1_epi32((int)isub+4*(int)nsub*(int)k)),
        nb128);        
#endif
#endif
    for (unsigned int j=0;j<ny;j++) {
      for (unsigned int i=0;i<nx;i++) {
        //image[i+j*nx] = 0.0f;
        // Project pixel centers, (x * cos(a) + y * sin(a)) / dr + wb
        tau = _mm_add_ps(
                _mm_div_ps(
                  _mm_add_ps(                                 
                    _mm_mul_ps(cang[k+isub*na/nsub/4],_mm_set_ps1((float) xs[i] )),
                    _mm_mul_ps(sang[k+isub*na/nsub/4],_mm_set_ps1((float) ys[j] ))),
                  LUT::dr128),                                       // Divide by dr
                wb128);                                              // Add 'r-offset'
                        
        // Get first r-index (minus one)
        ib_min = _mm_cvttps_epi32(_mm_floor_ps(                                      
                   _mm_sub_ps(                                       // Subtract max r-diff
                     tau,
                     taumax[k+isub*na/nsub/4])));

        for (unsigned int l=0;l<M;l++)  { // Loop over r-values that intersects pixel
          // R-indices
          __m128i ib = _mm_add_epi32(ib_min,_mm_set1_epi32(l+1)); 

          // Calculate data indices(r,angle)
#ifndef _LARGE_INDICES
          __m64 indices16 = 
            _mm_add_pi16(
              dataoffset,
              _mm_movepi64_pi64(
                _mm_packs_epi32(
                  ib,
                  _mm_setzero_si128()))); // _mm_cvtps_pi16(fib)); // no truncation for __m64 

          // Mask out indices that are too large to avoid access violation
          __m64 imask = _m_pand(_m_pcmpgtw(max_indices,indices16),_m_pcmpgtw(indices16,mones64));
          
          indices16 = _mm_sel_pi16(max_indices,indices16,imask);
#else
          register __m128i indices32 =
                      _mm_add_epi32(
                        dataoffset32,
                        ib);
#endif
          // Mask out some indices (0 <= r-index < nb) - TODO: use _mm_neg_pi
          __m128i rmask = _mm_and_si128(_mm_cmpgt_epi32(nb128,ib),_mm_cmpgt_epi32(ib,mones128));

          // R-values 
          __m128 res = _mm_mul_ps(
                         _mm_sub_ps(
                           _mm_cvtepi32_ps(ib),
                           tau),
                         LUT::dr128); 
          
          // Compute integrals
          __m128 sqint = sse_square_strip_int(res, k+isub*na/nsub/4);

          // Filter according to r-mask
#ifdef __linux__
          sqint = _mm_sel_ps(_mm_setzero_ps(), sqint, (__m128)rmask );
#else
          sqint = _mm_sel_ps(_mm_setzero_ps(), sqint, *(__m128*)&rmask );
#endif

          // Multiply by support mask - do this once afterwards, test this
                sqint = _mm_mul_ps(sqint,_mm_cvtpu8_ps(_mm_set1_pi8(mask[i+j*nx])));

#ifdef __linux__
          float cproj[4] __attribute__ ((aligned(16)));
#elif (defined(_MSC_VER) && defined(_WIN32))
          __declspec(align(16)) float cproj[4];
#endif

#ifndef _LARGE_INDICES
          cproj[0] = proj[_mm_extract_pi16(indices16,0)];
          cproj[1] = proj[_mm_extract_pi16(indices16,1)];
          cproj[2] = proj[_mm_extract_pi16(indices16,2)];
          cproj[3] = proj[_mm_extract_pi16(indices16,3)];
#else
          cproj[0] = proj[_mm_cvtsi128_si32(indices32)];
          indices32 = _mm_shuffle_epi32(indices32, _MM_SHUFFLE (0, 3, 2, 1));
          cproj[1] = proj[_mm_cvtsi128_si32(indices32)];
          indices32 = _mm_shuffle_epi32(indices32, _MM_SHUFFLE (0, 3, 2, 1));
          cproj[2] = proj[_mm_cvtsi128_si32(indices32)];
          indices32 = _mm_shuffle_epi32(indices32, _MM_SHUFFLE (0, 3, 2, 1));
          cproj[3] = proj[_mm_cvtsi128_si32(indices32)];
          indices32 = _mm_shuffle_epi32(indices32, _MM_SHUFFLE (0, 3, 2, 1));
#endif
          sqint = _mm_mul_ps(sqint,_mm_load_ps((float*)cproj));
          
                // TODO: Examine shift vs. shuffle
          /*          
          __m128 shifted = *(__m128*)&_mm_srli_si128(*(__m128i*)&sqint,8);
          sqint = _mm_add_ps(shifted,sqint);
          shifted = *(__m128*)&_mm_srli_si128(*(__m128i*)&sqint,4); 
          sqint = _mm_add_ps(shifted,sqint);
          shifted = _mm_load_ps1((float*)&image[i+j*nx]);
          sqint = _mm_add_ps(shifted,sqint);        
          _mm_store_ss((float*)&image[i+j*nx],sqint);
          */

          __m128 vsum = _mm_shuffle_ps(sqint, sqint, _MM_SHUFFLE (0, 1, 2, 3));
          vsum = _mm_add_ps(vsum,sqint); /* { 0+3, 1+2, 2+1, 3+0, } */
          sqint = _mm_movehl_ps(sqint, vsum);
          vsum = _mm_add_ps(vsum, sqint); /* { 0+3+2+1, 1+2+3+0, 2+1+2+1. 3+0+3+0 } */
          sqint = _mm_load_ps1((float*)&image[i+j*nx]);
          vsum = _mm_add_ps(vsum,sqint);        
          _mm_store_ss((float*)&image[i+j*nx],vsum);
        }
      }
    }
  }

  // Flush temporary vectors
  _mm_free(xs);
  _mm_free(ys);
  
  _mm_empty();
  _mm_sfence();
  return true;
}


bool tomo_back_full_rank(float *image, const unsigned char* mask, const size_t nx, const size_t ny,
                         const float dx, const float dy, const float offset_x, const float offset_y,
                         const float *proj, const size_t nr, const float offset_r, const size_t nsub, const size_t isub,
                         const size_t na, const size_t nthread) {
  
  unsigned int M;
  float wx, wy, wb;

  __m128 *cang, *sang, *taumax;
  
  __m128 wb128, tau;
  __m128i ib_min;
  
  __m64  nb64;
  __m128i nb128;
  
  if ((na % 4 != 0) || ((na / nsub) % 4 != 0))
    return false;

  _mm_empty(); // Before floating point operations

  // Number of strips that intersects a pixel for a fixed angle
  M = (unsigned int) ceil((dx * sqrt(2.0f) + LUT::sw) / LUT::dr);

  // Pixel centers
  wx = (nx-1.0f)/2.0f + offset_x;
  wy = (ny-1.0f)/2.0f + offset_y;

  // R center
  wb = (nr-1.0f)/2.0f + offset_r;

  _mm_empty(); // Before floating point operations
  
  float *xs = (float*) _mm_malloc( (size_t) ((nx) * sizeof(float)),16);
  float *ys = (float*) _mm_malloc( (size_t) ((ny) * sizeof(float)),16);
  
  for (unsigned int i=0;i<nx;i++) xs[i] = dx * ((float) i-wx);
  for (unsigned int i=0;i<ny;i++) ys[i] = dy * ((float) i-wy);

  wb128 = _mm_set_ps1(wb);
  nb64  = _mm_set1_pi16((short)nr);
  nb128  = _mm_set1_epi32((int)nr);

  // LUT references
  cang   = (__m128*) LUT::fcang;
  sang   = (__m128*) LUT::fsang;
  taumax = (__m128*) LUT::ftaumax;

#ifndef _LARGE_INDICES
  __m64 max_indices = _mm_set1_pi16(((short)nr)*((short)na)-1);
#else
  __m128i max_indices = _mm_set1_epi32(((int)nr)*((int)na)-1);
#endif

for (unsigned int k = 0 ; k < na/4/nsub ; k++) {
    // Data index (#angle * nb + #r)
#ifndef _LARGE_INDICES
    register __m64 dataoffset = 
      _mm_mullo_pi16(
        _mm_add_pi16(
          _mm_mullo_pi16(zero2three64i,_mm_set1_pi16((short)nsub)),
          _mm_set1_pi16((short)isub+4*(short)nsub*(short)k)), // LUT offsets[k]
        nb64);
#else
    register __m128i dataoffset32 =
      _mm_mul32_epi(
        _mm_add_epi32(
          _mm_mul32_epi(zero2three128i,_mm_set1_epi32((int)nsub)),
          _mm_set1_epi32((int)isub+4*(int)nsub*(int)k)),
        nb128);
#endif
    for (unsigned int j=0;j<ny;j++) {
      for (unsigned int i=0;i<nx;i++) {

        // Project pixel centers, (x * cos(a) + y * sin(a)) / dr + wb
        tau = _mm_add_ps(
                _mm_div_ps(
                  _mm_add_ps(                                 
                    _mm_mul_ps(cang[k+isub*na/nsub/4],_mm_set_ps1((float) xs[i] )),
                    _mm_mul_ps(sang[k+isub*na/nsub/4],_mm_set_ps1((float) ys[j] ))),
                  LUT::dr128),                                       // Divide by dr
                wb128);                                              // Add 'r-offset'
                        
        // Get first r-index (minus one)
        ib_min = _mm_cvttps_epi32(_mm_floor_ps(                                      
                   _mm_sub_ps(                                       // Subtract max r-diff
                     tau,
                     taumax[k+isub*na/nsub/4])));

        for (unsigned int l=0;l<M;l++)  { // Loop over r-values that intersects pixel
          // R-indices
          __m128i ib = _mm_add_epi32(ib_min,_mm_set1_epi32(l+1)); 

          // Calculate data indices(r,angle)
#ifndef _LARGE_INDICES
          __m64 indices16 = 
            _mm_add_pi16(
              dataoffset,
              _mm_movepi64_pi64(
                _mm_packs_epi32(
                  ib,
                  _mm_setzero_si128()))); // _mm_cvtps_pi16(fib)); // no truncation for __m64 

          // Mask out indices that are too large to avoid access violation
          __m64 imask = _m_pand(_m_pcmpgtw(max_indices,indices16),_m_pcmpgtw(indices16,mones64));
          
          indices16 = _mm_sel_pi16(max_indices,indices16,imask);
#else
          register __m128i indices32 =
                      _mm_add_epi32(
                        dataoffset32,
                        ib);
#endif
          // Mask out some indices (0 <= r-index < nb) - TODO: use _mm_neg_pi
          __m128i rmask = _mm_and_si128(_mm_cmpgt_epi32(nb128,ib),_mm_cmpgt_epi32(ib,mones128));

          // R-values 
          __m128 res = _mm_mul_ps(
                         _mm_sub_ps(
                           _mm_cvtepi32_ps(ib),
                           tau),
                         LUT::dr128); 
          
          // Compute integrals
          __m128 sqint = sse_square_strip_int(res, k+isub*na/nsub/4);

          // Filter according to r-mask
#ifdef __linux__
          sqint = _mm_sel_ps(_mm_setzero_ps(), sqint, (__m128)rmask );
#else
          sqint = _mm_sel_ps(_mm_setzero_ps(), sqint, *(__m128*)&rmask );
#endif

          // Multiply by support mask - do this once afterwards, test this
                sqint = _mm_mul_ps(sqint,_mm_cvtpu8_ps(_mm_set1_pi8(mask[i+j*nx])));

#ifdef __linux__
          float cproj[4] __attribute__ ((aligned(16)));
#elif (defined(_MSC_VER) && defined(_WIN32))
          __declspec(align(16)) float cproj[4];
#endif

#ifndef _LARGE_INDICES
          cproj[0] = proj[_mm_extract_pi16(indices16,0)];
          cproj[1] = proj[_mm_extract_pi16(indices16,1)];
          cproj[2] = proj[_mm_extract_pi16(indices16,2)];
          cproj[3] = proj[_mm_extract_pi16(indices16,3)];
#else
          cproj[0] = proj[_mm_cvtsi128_si32(indices32)];
          indices32 = _mm_shuffle_epi32(indices32, _MM_SHUFFLE (0, 3, 2, 1));
          cproj[1] = proj[_mm_cvtsi128_si32(indices32)];
          indices32 = _mm_shuffle_epi32(indices32, _MM_SHUFFLE (0, 3, 2, 1));
          cproj[2] = proj[_mm_cvtsi128_si32(indices32)];
          indices32 = _mm_shuffle_epi32(indices32, _MM_SHUFFLE (0, 3, 2, 1));
          cproj[3] = proj[_mm_cvtsi128_si32(indices32)];
          indices32 = _mm_shuffle_epi32(indices32, _MM_SHUFFLE (0, 3, 2, 1));
#endif

          sqint = _mm_mul_ps(sqint,_mm_load_ps((float*)cproj));
          
          __m128 vsum = _mm_shuffle_ps(sqint, sqint, _MM_SHUFFLE (0, 1, 2, 3));
          vsum = _mm_add_ps(vsum,sqint); /* { 0+3, 1+2, 2+1, 3+0, } */
          sqint = _mm_movehl_ps(sqint, vsum);
          vsum = _mm_add_ps(vsum, sqint); /* { 0+3+2+1, 1+2+3+0, 2+1+2+1. 3+0+3+0 } */
          sqint = _mm_load_ps1((float*)&image[i+j*nx]);
          vsum = _mm_add_ps(vsum,sqint);        
          _mm_store_ss((float*)&image[i+j*nx],vsum);
        }
      }
    }
  }

  // Flush temporary vectors
  _mm_free(xs);
  _mm_free(ys);
  
  _mm_empty();
  _mm_sfence();
  return true;
}

