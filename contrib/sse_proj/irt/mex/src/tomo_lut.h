// $Id: tomo_lut.h,v 1.6 2007/08/03 13:55:21 jmh Exp $

#ifndef TOMO_LUT_H
#define TOMO_LUT_H

#ifndef _NO_SSE
 #include "extintrin.h"
#else
 #include <float.h>
#endif

namespace LUT {

// Look-up tables
extern volatile float* fcang;
extern volatile float* fsang;
extern volatile float* fdmax;
extern volatile float* fdbrk;
extern volatile float* ftaumax;
extern volatile float* flang;

// Used by square_strip_int (SSE not available)
extern float dx;
extern float sw;
extern float norm;

// Used by sse_square_strip_int
extern __m128 dx128;
extern __m128 sw128;
extern __m128 norm128;

// Used by projectors
extern float dr;

// Used by projectors (SSE not available)
extern __m128 dr128;

 
bool initialize_LUTs(size_t nangles, size_t nsub, float orbit_low, float orbit_high,
                     float idx, float idr, float isw);

/*
extern inline bool initialize_LUTs(size_t inangles, float iorbit_low, float iorbit_high,
                     float idx, float idr, float isw);

*/
void free_LUTs();

} // end namespace LUT

#endif
