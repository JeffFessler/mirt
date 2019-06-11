// $Id: tomo_strip.cpp,v 1.11 2007/08/07 21:43:29 jmh Exp $
#include "tomo_strip.h"
#include "tomo_lut.h"

#ifdef __ICC
 #pragma warning(disable:981) // Operands are evaluated in unspecified order
 #pragma warning(disable:1599) // Hides variables
#endif

#ifndef _NO_SSE

// Must not be called directly

// extern inline __m128 sse_square_strip_int(__m128 cradius, int k) { 
// linux
extern __m128 sse_square_strip_int(__m128 cradius, size_t k) {

  __m128 a,b,mask,temp,a1,b1;

  __m128 output;

  register __m128* dmax;
  register __m128* dbrk;
  register __m128* lang;

  // LUT references
  dmax = (__m128*) LUT::fdmax;
  dbrk = (__m128*) LUT::fdbrk;
  lang = (__m128*) LUT::flang;

  a = _mm_div_ps(_mm_sub_ps(cradius,_mm_mul_ps(LUT::sw128,half128)),LUT::dx128);
  b = _mm_div_ps(_mm_add_ps(cradius,_mm_mul_ps(LUT::sw128,half128)),LUT::dx128);

  // First 'triangle or rhomb', g1
  a1 = _mm_max_ps(a,_mm_neg_ps(dmax[k]));
  b1 = _mm_min_ps(b,_mm_neg_ps(dbrk[k]));
  mask = _mm_cmpgt_ps(b1,a1);

  temp = _mm_mul_ps(
           _mm_div_ps(_mm_mul_ps(lang[k],half128),_mm_sub_ps(dmax[k],dbrk[k])),
             _mm_sub_ps(_mm_square_ps(_mm_add_ps(b1,dmax[k])),
                        _mm_square_ps(_mm_add_ps(a1,dmax[k]))));
     
  output = _mm_sel_ps(_mm_setzero_ps(), temp, mask );

  // Rectangular region, g2
  a1 = _mm_max_ps(a,_mm_neg_ps(dbrk[k]));
  b1 = _mm_min_ps(b,dbrk[k]);

  mask = _mm_cmpgt_ps(b1,a1);

  temp = _mm_add_ps(
           _mm_mul_ps(lang[k],_mm_sub_ps(b1,a1)),
           output);
  output = _mm_sel_ps(output, temp, mask);
  
  // Last 'triangle or rhomb', g3
  a1 = _mm_max_ps(a,dbrk[k]);
  b1 = _mm_min_ps(b,dmax[k]);
  mask = _mm_cmpgt_ps(b1,a1);
  temp = _mm_add_ps(
           _mm_mul_ps(
             _mm_div_ps(_mm_mul_ps(lang[k],half128),_mm_sub_ps(dmax[k],dbrk[k])),
             _mm_sub_ps(_mm_square_ps(_mm_sub_ps(a1,dmax[k])),
                        _mm_square_ps(_mm_sub_ps(b1,dmax[k])))),
           output); 

  output = _mm_sel_ps(output, temp, mask);

  // Normalize by strip-width
  output = _mm_mul_ps(output,LUT::norm128); 

  // Flush SSE registers
  // _mm_empty();

  return output;
}

#else
inline static float square_strip_int(float cradius, size_t k) {

  float a,b,temp,a1,b1;
  bool mask;

  float output;

  register float* dmax;
  register float* dbrk;
  register float* lang;

  // LUT references
  dmax = (float*) LUT::fdmax;
  dbrk = (float*) LUT::fdbrk;
  lang = (float*) LUT::flang;

  a = (cradius-0.5f*LUT::sw)/LUT::dx;
  b = (cradius+0.5f*LUT::sw)/LUT::dx;

  // First 'triangle or rhomb', g1
  a1 = std::max(a,-dmax[k]);
  b1 = std::max(b,-dbrk[k]);

  mask = (b1 >= a1);

  temp = (0.5f*lang[k]/(dmax[k]-dbrk[k])) *((b1+dmax[k])*(b1+dmax[k])- (a1+dmax[k])*(a1+dmax[k]));

  // TODO: Avoid branches
  if (mask)
    output = temp;
  else
    output = 0.0f;

  // Rectangular region, g2
  a1 = std::max(a,-dbrk[k]);
  b1 = std::min(b,dbrk[k]);

  mask = (b1 >= a1);

  temp = (lang[k]*(b1-a1))+output;
  if (mask)
    output = temp;
  
  // Last 'triangle or rhomb', g3
  a1 = std::max(a,dbrk[k]);
  b1 = std::min(b,dmax[k]);
  mask = (b1 >= a1);

  temp = ((0.5f*lang[k] / (dmax[k]-dbrk[k]))*
   (( ((a1-dmax[k])*(a1-dmax[k]) )-
      ((b1-dmax[k])*(b1-dmax[k])) ))+output); 

  if (mask)
    output = temp;

  // Normalize by strip-width
  output = output*LUT::norm;

  return output;
}

#endif
