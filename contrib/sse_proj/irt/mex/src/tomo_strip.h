
// $Id: tomo_strip.h,v 1.9 2007/08/07 21:43:29 jmh Exp $
#ifndef TOMO_STRIP_H
#define TOMO_STRIP_H

#ifndef _NO_SSE
 #include "extintrin.h"
#endif

#if (defined(_MSC_VER) && defined(_WIN32))

#ifndef _NO_SSE
 extern __forceinline __m128 sse_square_strip_int(__m128 cradius, size_t k);
#else
 static __forceinline float square_strip_int(float cradius, size_t k);
#endif

#elif defined(__linux__)

#ifndef _NO_SSE
 extern __m128 sse_square_strip_int(__m128 cradius, size_t k);// __attribute__ ((always_inline));
#else
 static inline float square_strip_int(float cradius, int k) __attribute__ ((always_inline));
#endif

#endif

#endif
