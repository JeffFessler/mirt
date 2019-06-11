// $Id: extintrin.h,v 1.27 2007/08/07 21:43:29 jmh Exp $
#ifndef EXTINTRIN_H
#define EXTINTRIN_H

#ifdef _MSC_VER
 #pragma warning( disable : 4238 ) // suppress warning C4238: nonstandard extension used : class rvalue used as lvalue
 #pragma warning( disable : 4127 ) // suppress warning C4127: conditional expression is constant
#endif

// MMX
//#include <mmintrin.h>

// SSE
#include <xmmintrin.h>
// SSE2
#include <emmintrin.h>
// SSE3
//#include <pmmintrin.h>

//#define _mm_extract_pi32(A, N) __builtin_arm_textrmsw ((__v2si)(A), (N))

#ifndef _MM_SHUFFLE
 #define _MM_SHUFFLE(z,y,x,w)    ((z << 6) | (y << 4) | (x << 2) | w)
#endif

// STL headers first
#ifndef Mmex
 #include <iostream>
#endif

// Standard C++ headers
#include <stdio.h>
#include <stdlib.h>
#include <float.h>

const __m128 min128  = _mm_set_ps1(FLT_MAX);  
const __m128 max128  = _mm_set_ps1(FLT_MIN);  
const __m128 eps128  = _mm_set_ps1(FLT_EPSILON);
const __m128 half128 = _mm_set_ps1(0.5f);
const __m64 ones64 = _mm_set1_pi16(1);
const __m64 mones64 = _mm_set1_pi16(-1);
const __m128 ones128 = _mm_set_ps1(1.0f);
const __m128i ones128i = _mm_set1_epi32(1);
const __m128i zero2three128i = _mm_set_epi32(3,2,1,0);
const __m64 zero2three64i = _mm_set_pi16(3,2,1,0);
const __m128i mones128 = _mm_set1_epi32(-1);

#if (defined(_MSC_VER) && defined(_WIN32))
  static __forceinline __m128 _mm_fabs_ps(__m128 x);
  static __forceinline __m128 _mm_neg_ps(__m128 x);
  static __forceinline __m128i _mm_neg_si(__m128i x);
  static __forceinline __m64 _mm_neg_si(__m64 x);
  static __forceinline __m128 _mm_trunc_ps(__m128 v );
  static __forceinline __m128 _mm_floor_ps(__m128 v );
  static __forceinline __m128 _mm_ceil_ps(__m128 v );
  static __forceinline __m128 _mm_sel_ps(__m128 a, __m128 b, __m128 mask );
  static __forceinline __m128 _mm_square_ps(__m128 a);
  static __forceinline int _mm_any_eq( __m128 a, __m128 b );
  static __forceinline float _mm_vsum_ps(__m128 v);
  static __m128 _mm_pow_fixed_point_ps(__m128 x, int exponent);
  static __forceinline __m128 _mm_pow_fixed_point_ps(__m128 x, float exponent);
#else
  static inline __m128 _mm_fabs_ps(__m128 x) __attribute__ ((always_inline));
  static inline __m128 _mm_neg_ps(__m128 x) __attribute__ ((always_inline));
  static inline __m128i _mm_neg_si(__m128i x) __attribute__ ((always_inline));
  static inline __m64 _mm_neg_si(__m64 x) __attribute__ ((always_inline));
  static inline __m128 _mm_trunc_ps(__m128 v ) __attribute__ ((always_inline));
  static inline __m128 _mm_floor_ps(__m128 v ) __attribute__ ((always_inline));
  static inline __m128 _mm_ceil_ps(__m128 v )  __attribute__ ((always_inline));
  static inline __m128 _mm_sel_ps(__m128 a, __m128 b, __m128 mask ) __attribute__ ((always_inline));
  static inline __m128 _mm_square_ps(__m128 a) __attribute__ ((always_inline));
  static inline    int _mm_any_eq( __m128 a, __m128 b ) __attribute__ ((always_inline));
  static inline float _mm_vsum_ps(__m128 v) __attribute__ ((always_inline));
  static __m128 _mm_pow_fixed_point_ps(__m128 x, int exponent);
  static inline __m128 _mm_pow_fixed_point_ps(__m128 x, float exponent) __attribute__ ((always_inline));
#endif

static inline __m128 _mm_fabs_ps(__m128 x)
{
#if (defined(_MSC_VER) && defined(_WIN32))
  static __declspec(align(16)) int clear_signmask[4] = {0x7FFFFFFFL,0x7FFFFFFFL,0x7FFFFFFFL,0x7FFFFFFFL};
#elif defined(__linux__)
 	static int clear_signmask[4] __attribute__ ((aligned (16)))=
    {0x7FFFFFFFL,0x7FFFFFFFL,0x7FFFFFFFL,0x7FFFFFFFL};
#endif
 	return _mm_and_ps(x,_mm_load_ps((float *) clear_signmask));
}

static inline __m128 _mm_fabs_shift_ps(__m128 x)
{
#if (defined(_MSC_VER) && defined(_WIN32) && !defined(__ICC))
  __m128 b = *(__m128*) &_mm_srli_epi32( _mm_slli_epi32( *(__m128i*) &x, 1 ), 1 );
#elif (defined(__linux__) || defined(__ICC))   
  __m128 b = (__m128) _mm_srli_epi32( _mm_slli_epi32( (__m128i) x, 1 ), 1 );
#endif
  return b;
}

static inline __m128i _mm_neg_si(__m128i x) {
  register unsigned int uiMask = 0xFFFFFFFFL;
  return _mm_add_epi32(_mm_xor_si128(x, _mm_set1_epi32(uiMask)),ones128i);
}
static inline __m64 _mm_neg_si(__m64 x) {
  register unsigned short uiMask = 0xFFFF;
  return _mm_adds_pi16(_mm_xor_si64(x,_mm_set1_pi16(uiMask)),ones64);  
}

// Flip sign of a packed single
static inline __m128 _mm_neg_ps(__m128 x)
{
#if (defined(_MSC_VER) && defined(_WIN32))
  static __declspec(align(16)) int clear_signmask[4] = {0x80000000L,0x80000000L,0x80000000L,0x80000000L};
#elif defined(__linux__)
  static int clear_signmask[4] __attribute__ ((aligned (16)))=
    {0x80000000L,0x80000000L,0x80000000L,0x80000000L};
#endif
//  _mm_empty();
  return _mm_xor_ps(x,_mm_load_ps((float *) clear_signmask));
   
/*
  register unsigned int uiMask = 0x80000000L;
  return _mm_xor_ps( _mm_set_ps1( *(float*)(&uiMask) ), x);
  return _mm_xor_ps( _mm_set_ps1((int)uiMask, x);
*/
}

// Similar to _mm_cvttepi32_ps
static inline __m128 _mm_trunc_ps(__m128 v ) {

#if (defined(_MSC_VER) && defined(_WIN32))
     // register unsigned int twoTo23s = 0x4b000000L;
     // __m128 twoTo23 = _mm_load_ps1((float*)(&twoTo23s));
     static const __m128 twoTo23 = _mm_set_ps1(1 << 23);
#elif defined(__linux__)
    static const __m128 twoTo23 = (__m128){ 0x1.0p23f, 0x1.0p23f, 0x1.0p23f, 0x1.0p23f };
#endif

    // b = _mm_fabs_ps(v)
    __m128 b = _mm_fabs_shift_ps(v);
    // The essence of the floor routine
    __m128 d = _mm_sub_ps( _mm_add_ps( b, twoTo23 ), twoTo23 );
    // 1 if v >= 2**23
    __m128 largeMaskE = _mm_cmpgt_ps( b, twoTo23 );
    // Check for possible off by one error
    __m128 g = _mm_cmplt_ps( b, d );
    // Convert positive check result to -1.0, negative to 0.0
#if (defined(_MSC_VER) && defined(_WIN32))
    __m128 h = _mm_cvtepi32_ps( *(__m128i*) &g );
#elif defined(__linux__)
    __m128 h = _mm_cvtepi32_ps( (__m128i) g );
#endif
    // Add in the error if there is one
    __m128 t = _mm_add_ps( d, h );

    // Put the sign bit back
#if (defined(_MSC_VER) && defined(_WIN32) && !defined(__ICC))
    __m128 sign = *((__m128*) &_mm_slli_epi32( _mm_srli_epi32( *(__m128i*) &v, 31), 31 ));
#elif (defined(__linux__) || defined(__ICC))
    __m128 sign = (__m128) _mm_slli_epi32( _mm_srli_epi32( (__m128i) v, 31), 31 );
#endif
    
    t = _mm_or_ps( t, sign );

    // Select between output result and input value based on _mm_fabs_ps(v) >= 2**23
    v = _mm_and_ps( v, largeMaskE );
    t = _mm_andnot_ps( largeMaskE, t );

    return _mm_or_ps( t, v );
}

inline __m128 _mm_floor_ps(__m128 v ) {

#if (defined(_MSC_VER) && defined(_WIN32))
    static const __m128 twoTo23 = _mm_set_ps1(1 << 23);
#elif (defined(__linux__) || defined(__ICC))
    static const __m128 twoTo23 = (__m128){ 0x1.0p23f, 0x1.0p23f, 0x1.0p23f, 0x1.0p23f };
#endif

    // b = _mm_fabs_ps(v)
    __m128 b = _mm_fabs_shift_ps(v);
    // The essence of the floor routine
    __m128 d = _mm_sub_ps( _mm_add_ps( _mm_add_ps( _mm_sub_ps( v, twoTo23 ), twoTo23 ), twoTo23 ), twoTo23 );
    // 1 if v >= 2**23
    __m128 largeMaskE = _mm_cmpgt_ps( b, twoTo23 );
    // Check for possible off by one error
    __m128 g = _mm_cmplt_ps( v, d );
    // Convert positive check result to -1.0, negative to 0.0
#if (defined(_MSC_VER) && defined(_WIN32) && !defined(__ICC))
    __m128 h = _mm_cvtepi32_ps( *(__m128i*) &g );
#elif (defined(__linux__)|| defined(__ICC))
    __m128 h = _mm_cvtepi32_ps( (__m128i) g );
#endif
    // Add in the error if there is one
    __m128 t = _mm_add_ps( d, h );
    // Select between output result and input value based on v >= 2**23
    v = _mm_and_ps( v, largeMaskE );
    t = _mm_andnot_ps( largeMaskE, t );

    return _mm_or_ps( t, v ); 
}
   
static inline __m128 _mm_ceil_ps(__m128 v ) {

#if (defined(_MSC_VER) && defined(_WIN32))
    static const __m128 twoTo23 = _mm_set_ps1(1 << 23);
#elif defined(__linux__)
    static const __m128 twoTo23 = (__m128){ 0x1.0p23f, 0x1.0p23f, 0x1.0p23f, 0x1.0p23f };
#endif

  __m128 d = _mm_sub_ps( _mm_add_ps( _mm_add_ps( _mm_sub_ps( v, twoTo23 ), twoTo23 ), twoTo23 ), twoTo23 );
  __m128 largeMaskE = _mm_cmpgt_ps( v, twoTo23 );
  __m128 g = _mm_cmpgt_ps( v, d );
#if (defined(_MSC_VER) && defined(_WIN32) && !defined(__ICC))
  __m128 h = _mm_cvtepi32_ps( *(__m128i*) &g );
#elif (defined(__linux__)|| defined(__ICC))
  __m128 h = _mm_cvtepi32_ps((__m128i)g);
#endif
  __m128 t = _mm_sub_ps( d, h );

  v = _mm_and_ps( v, largeMaskE );
  t = _mm_andnot_ps( largeMaskE, t );
  return _mm_or_ps( t, v );
}

static inline __m128 _mm_sel_ps(__m128 a, __m128 b, __m128 mask ) {
    b = _mm_and_ps( b, mask );
    a = _mm_andnot_ps( mask, a );
    return _mm_or_ps( a, b );
}

static inline __m64 _mm_sel_pi16(__m64 a, __m64 b, __m64 mask ) {
    b = _m_pand( b, mask );
    a = _m_pandn( mask, a );
    return _m_por( a, b );
}

static inline __m128i _mm_sel_pi32(__m128i a, __m128i b, __m128i mask ) {
    b = _mm_and_si128( b, mask );
    a = _mm_andnot_si128( mask, a );
    return _mm_or_si128( a, b );
}

//
static inline __m128 _mm_square_ps(__m128 a) {
  return _mm_mul_ps(a,a);
}

// Test a==b for each float in a & b
static inline int _mm_any_eq( __m128 a, __m128 b ) {
    register __m128 mask = _mm_cmpeq_ps( a, b );
    //copy top bit of each result to maskbits
    return _mm_movemask_ps( mask ) != 0;
}

static __m128 _mm_pow_fixed_point_ps(__m128 x, int exponent)
{
	__m128 rslt=ones128;									// x^0=1.0
	int xp=abs(exponent);
	if (xp & 3)												// fraction present?
	{
		__m128 sq_rt=_mm_sqrt_ps(x);
		if (xp & 1)											// .25?
			rslt=_mm_sqrt_ps(sq_rt);						// x^.25
		if (xp & 2)
			rslt=_mm_mul_ps(rslt,sq_rt);
	}
	xp>>=2;													// strip fraction
	__m128 curpower=x;										// curpower iterates through  x,x^2,x^4,x^8,x^16...

	while(1)
	{
		if (xp & 1)
			rslt=_mm_mul_ps(rslt,curpower);
		xp>>=1;
		if (xp)
			curpower=_mm_mul_ps(curpower,curpower);
		else
			break;
	}
	if (exponent<0)
		return _mm_rcp_ps(rslt);							// pow(x,-b)=1/pow(x,b)
	else
		return rslt;
}

// _mm_pow_fixed_point_ps - raise an sse register to a power.  This is analogous to the C pow()
// function, with some restictions: fractional exponents are only handled with 2 bits of precision.
// Basically, fractions of 0,.25,.5, and .75 are handled. _mm_pow_fixed_point_ps(x,.30) will be
// the same as _mm_pow_fixed_point_ps(x,.25). negative and fractional powers are handled by the
// SSE reciprocal and square root approximation instructions and so are not especially accurate --
// Note that this routine does not raise numeric exceptions because it uses SSE--- This routine is
// O(log2(exponent)).
inline __m128 _mm_pow_fixed_point_ps(__m128 x, float exponent)
{
	return _mm_pow_fixed_point_ps(x,(int) (4.0*exponent));
}



// _mm_mul32_epi: multiply 2 signed or unsigned 32 bit integer vector
static inline __m128i _mm_mul32_epi(const __m128i a, const __m128i b) {
  __m128i a13, b13, prod02, prod13, prod01, prod23, prod0123;
  a13 = _mm_shuffle_epi32(a, 0xF5);              // (-,a3,-,a1)
  b13 = _mm_shuffle_epi32(b, 0xF5);              // (-,b3,-,b1)
  prod02 = _mm_mul_epu32(a, b);                  // (-,a2*b2,-,a0*b0)
  prod13 = _mm_mul_epu32(a13, b13);              // (-,a3*b3,-,a1*b1)
  prod01 = _mm_unpacklo_epi32(prod02,prod13);    // (-,-,a1*b1,a0*b0)
  prod23 = _mm_unpackhi_epi32(prod02,prod13);    // (-,-,a3*b3,a2*b2)
  prod0123 = _mm_unpacklo_epi64(prod01,prod23);  // (ab3,ab2,ab1,ab0)
  return prod0123;
}

//SSE2: multiply a * b and return full result in high and low result
static void _mm_mul_full( __m128i &highResult, __m128i &lowResult, __m128i a, __m128i b)
{
    __m128i hi = _mm_mulhi_epi16( a, b );        // (a7*b7[16:31],a6*b6[16:31],a5*b5[16:31],a4*b4[16:31],a3*b3[16:31],a2*b2[16:31],a1*b1[16:31],a0*b0[16:31])
    __m128i low = _mm_mullo_epi16( a, b );       // (a7*b7[0:15] ,a6*b6[0:15] ,a5*b5[0:15] ,a4*b4[0:15],a3*b3[0:15] ,a2*b2[0:15] ,a1*b1[0:15] ,a0*b0[0:15])
    highResult = _mm_unpacklo_epi16( hi, low ); // (a3*b3[0:15],a3*b3[16:31],a2*b2[0:15],a2*b2[16:31],a1*b1[0:15],a1*b1[16:31],a0*b0[0:15],a0*b0[16:31])
    lowResult = _mm_unpackhi_epi16( hi, low );  // (a7*b7[0:15],a7*b7[16:31],a6*b6[0:15],a6*b6[16:31],a5*b5[0:15],a5*b5[16:31],a4*b4[0:15],a4*b4[16:31])
}

#ifndef __linux__
static inline bool operator==(__m128 a, __m128 b) {
  register __m128 mask = _mm_cmpeq_ps( a, b );
  return _mm_movemask_ps( mask ) == 0x0F;
}
#endif

#ifndef Mmex

inline std::ostream& operator<<(std::ostream& out, __m128 vf) {
#if (defined(_MSC_VER) && defined(_WIN32) && !defined(__ICC))
  out.flush();
  out << "__m128[3]: " << (float) vf.m128_f32[3] << " __m128[2]: " << (float) vf.m128_f32[2] << " __m128[1]: " << (float) vf.m128_f32[1] << " __m128[0]: " << (float) vf.m128_f32[0] << " " << std::endl;
#elif (defined(_WIN32) && defined(__ICC))
  __declspec(align(16)) float output[4];
  _mm_stream_ps((float*)output,vf);
  _mm_empty();
  out << "__m128[3]: " << (float) output[3] << " __m128[2]: " << (float) output[2] \
      << " __m128[1]: " << (float) output[1] << " __m128[0]: " << (float) output[0] << " " << std::endl;  
#elif defined(__linux__)
  float output[4] __attribute__ ((aligned (16)));
  _mm_stream_ps((float*)output,vf);
  _mm_empty();
  out << "__m128[3]: " << (float) output[3] << " __m128[2]: " << (float) output[2] \
      << " __m128[1]: " << (float) output[1] << " __m128[0]: " << (float) output[0] << " " << std::endl;
#endif
  return out;
}

inline std::ostream& operator<<(std::ostream& out, __m128i vi) {
#if (defined(_MSC_VER) && defined(_WIN32) && !defined(__ICC))
  _mm_empty();
  out << "__m128i[3]: " << vi.m128i_i32[3] << " __m128i[2]: " << vi.m128i_i32[2] \
      << " __m128i[1]: " << vi.m128i_i32[1] << " __m128i[0]: " << vi.m128i_i32[0] << std::endl;
#elif (defined(_WIN32) && defined(__ICC))
  __declspec(align(16)) int output[4];
  _mm_stream_si128((__m128i*)output,vi);
  _mm_empty();
  out << "__m128[3]: " << (int) output[3] << " __m128[2]: " << (int) output[2] \
      << " __m128[1]: " << (int) output[1] << " __m128[0]: " << (int) output[0] << " " << std::endl;  
#elif defined(__linux__)
  int output[4] __attribute__ ((aligned (16)));
  _mm_stream_si128((__m128i*)output,vi);
  _mm_empty();
  out << "__m128i[3]: " << (int) output[3] << " __m128i[2]: " << (int) output[2] \
      << " __m128i[1]: " << (int) output[1] << " __m128i[0]: " << (int) output[0] << " " << std::endl;
#endif
  return out;
}

inline std::ostream& operator<<(std::ostream& out, __m64 vi) {
#if (defined(_MSC_VER) && defined(_WIN32) && !defined(__ICC))
  _mm_empty();
  out << "__m64[3]: " << vi.m64_i16[3] << " __m64[2]: " << vi.m64_i16[2] \
      << " __m64[1]: " << vi.m64_i16[1] << " __m64[0]: " << vi.m64_i16[0] << std::endl;
#elif (defined(_WIN32) && defined(__ICC))
  __declspec(align(16)) short int output[4];
  _mm_stream_pi((__m64*)output,vi);
  _mm_empty(); // before stream operations
  out << "__m64[3]: " << output[3] << " __m64[2]: " << output[2] \
      << " __m64[1]: " << output[1] << " __m64[0]: " << output[0] << std::endl;
#elif defined(__linux__)
  short int output[4] __attribute__ ((aligned (16)));
  _mm_stream_pi((__m64*)output,vi);
  _mm_empty();
  out << "__m64[3]: " << output[3] << " __m64[2]: " << output[2] \
      << " __m64[1]: " << output[1] << " __m64[0]: " << output[0] << std::endl;
#endif
  return out;
}

#endif

#ifdef _WIN32

template<unsigned N, unsigned I>
class CosineSSE {
public:
  static inline __m128 cos() {

  return _mm_set_ps(1-(I*2*M_PI/N)*(I*2*M_PI/N)/2*(1-(I*2*M_PI/N)*(I*2*M_PI/N)/3/4*
        (1-(I*2*M_PI/N)*(I*2*M_PI/N)/5/6*(1-(I*2*M_PI/N)*(I*2*M_PI/N)/7/8*
        (1-(I*2*M_PI/N)*(I*2*M_PI/N)/9/10*(1-(I*2*M_PI/N)*(I*2*M_PI/N)/11/12*
        (1-(I*2*M_PI/N)*(I*2*M_PI/N)/13/14*(1-(I*2*M_PI/N)*(I*2*M_PI/N)/15/16*
        (1-(I*2*M_PI/N)*(I*2*M_PI/N)/17/18*(1-(I*2*M_PI/N)*(I*2*M_PI/N)/19/20*
        (1-(I*2*M_PI/N)*(I*2*M_PI/N)/21/22*(1-(I*2*M_PI/N)*(I*2*M_PI/N)/23/24
        ))))))))))),1-((I+1)*2*M_PI/N)*((I+1)*2*M_PI/N)/2*(1-((I+1)*2*M_PI/N)*((I+1)*2*M_PI/N)/3/4*
        (1-((I+1)*2*M_PI/N)*((I+1)*2*M_PI/N)/5/6*(1-((I+1)*2*M_PI/N)*((I+1)*2*M_PI/N)/7/8*
        (1-((I+1)*2*M_PI/N)*((I+1)*2*M_PI/N)/9/10*(1-((I+1)*2*M_PI/N)*((I+1)*2*M_PI/N)/11/12*
        (1-((I+1)*2*M_PI/N)*((I+1)*2*M_PI/N)/13/14*(1-((I+1)*2*M_PI/N)*((I+1)*2*M_PI/N)/15/16*
        (1-((I+1)*2*M_PI/N)*((I+1)*2*M_PI/N)/17/18*(1-((I+1)*2*M_PI/N)*((I+1)*2*M_PI/N)/19/20*
        (1-((I+1)*2*M_PI/N)*((I+1)*2*M_PI/N)/21/22*(1-((I+1)*2*M_PI/N)*((I+1)*2*M_PI/N)/23/24
        ))))))))))),1-((I+2)*2*M_PI/N)*((I+2)*2*M_PI/N)/2*(1-((I+2)*2*M_PI/N)*((I+2)*2*M_PI/N)/3/4*
        (1-((I+2)*2*M_PI/N)*((I+2)*2*M_PI/N)/5/6*(1-((I+2)*2*M_PI/N)*((I+2)*2*M_PI/N)/7/8*
        (1-((I+2)*2*M_PI/N)*((I+2)*2*M_PI/N)/9/10*(1-((I+2)*2*M_PI/N)*((I+2)*2*M_PI/N)/11/12*
        (1-((I+2)*2*M_PI/N)*((I+2)*2*M_PI/N)/13/14*(1-((I+2)*2*M_PI/N)*((I+2)*2*M_PI/N)/15/16*
        (1-((I+2)*2*M_PI/N)*((I+2)*2*M_PI/N)/17/18*(1-((I+2)*2*M_PI/N)*((I+2)*2*M_PI/N)/19/20*
        (1-((I+2)*2*M_PI/N)*((I+2)*2*M_PI/N)/21/22*(1-((I+2)*2*M_PI/N)*((I+2)*2*M_PI/N)/23/24
        ))))))))))),1-((I+3)*2*M_PI/N)*((I+3)*2*M_PI/N)/2*(1-((I+3)*2*M_PI/N)*((I+3)*2*M_PI/N)/3/4*
        (1-((I+3)*2*M_PI/N)*((I+3)*2*M_PI/N)/5/6*(1-((I+3)*2*M_PI/N)*((I+3)*2*M_PI/N)/7/8*
        (1-((I+3)*2*M_PI/N)*((I+3)*2*M_PI/N)/9/10*(1-((I+3)*2*M_PI/N)*((I+3)*2*M_PI/N)/11/12*
        (1-((I+3)*2*M_PI/N)*((I+3)*2*M_PI/N)/13/14*(1-((I+3)*2*M_PI/N)*((I+3)*2*M_PI/N)/15/16*
        (1-((I+3)*2*M_PI/N)*((I+3)*2*M_PI/N)/17/18*(1-((I+3)*2*M_PI/N)*((I+3)*2*M_PI/N)/19/20*
        (1-((I+3)*2*M_PI/N)*((I+3)*2*M_PI/N)/21/22*(1-((I+3)*2*M_PI/N)*((I+3)*2*M_PI/N)/23/24
        ))))))))))));
  }
};

template<unsigned N, unsigned I>
class AbsCosineSSE {
public:
  static inline __m128 abscos() {
    return _mm_fabs_ps(CosineSSE<N,I>::cos());
  }
};
// For absolute values check if 0 < N/I < PI and ad a minus if not
template<unsigned N, unsigned I>
class SineSSE {
public:
  static inline __m128 sin() {

    return _mm_set_ps((I*2*M_PI/N)*(1-(I*2*M_PI/N)*(I*2*M_PI/N)/2/3*(1-(I*2*M_PI/N)*
           (I*2*M_PI/N)/4/5*(1-(I*2*M_PI/N)*(I*2*M_PI/N)/6/7*(1-(I*2*M_PI/N)*
           (I*2*M_PI/N)/8/9*(1-(I*2*M_PI/N)*(I*2*M_PI/N)/10/11*(1-(I*2*M_PI/N)*
           (I*2*M_PI/N)/12/13*(1-(I*2*M_PI/N)*(I*2*M_PI/N)/14/15*
           (1-(I*2*M_PI/N)*(I*2*M_PI/N)/16/17*
           (1-(I*2*M_PI/N)*(I*2*M_PI/N)/18/19*(1-(I*2*M_PI/N)*
           (I*2*M_PI/N)/20/21)))))))))),
           ((I+1)*2*M_PI/N)*(1-((I+1)*2*M_PI/N)*((I+1)*2*M_PI/N)/2/3*(1-((I+1)*2*M_PI/N)*
           ((I+1)*2*M_PI/N)/4/5*(1-((I+1)*2*M_PI/N)*((I+1)*2*M_PI/N)/6/7*(1-((I+1)*2*M_PI/N)*
           ((I+1)*2*M_PI/N)/8/9*(1-((I+1)*2*M_PI/N)*((I+1)*2*M_PI/N)/10/11*(1-((I+1)*2*M_PI/N)*
           ((I+1)*2*M_PI/N)/12/13*(1-((I+1)*2*M_PI/N)*((I+1)*2*M_PI/N)/14/15*
           (1-((I+1)*2*M_PI/N)*((I+1)*2*M_PI/N)/16/17*
           (1-((I+1)*2*M_PI/N)*((I+1)*2*M_PI/N)/18/19*(1-((I+1)*2*M_PI/N)*
           ((I+1)*2*M_PI/N)/20/21)))))))))),
           ((I+2)*2*M_PI/N)*(1-((I+2)*2*M_PI/N)*((I+2)*2*M_PI/N)/2/3*(1-((I+2)*2*M_PI/N)*
           ((I+2)*2*M_PI/N)/4/5*(1-((I+2)*2*M_PI/N)*((I+2)*2*M_PI/N)/6/7*(1-((I+2)*2*M_PI/N)*
           ((I+2)*2*M_PI/N)/8/9*(1-((I+2)*2*M_PI/N)*((I+2)*2*M_PI/N)/10/11*(1-((I+2)*2*M_PI/N)*
           ((I+2)*2*M_PI/N)/12/13*(1-((I+2)*2*M_PI/N)*((I+2)*2*M_PI/N)/14/15*
           (1-((I+2)*2*M_PI/N)*((I+2)*2*M_PI/N)/16/17*
           (1-((I+2)*2*M_PI/N)*((I+2)*2*M_PI/N)/18/19*(1-((I+2)*2*M_PI/N)*
           ((I+2)*2*M_PI/N)/20/21)))))))))),
           ((I+3)*2*M_PI/N)*(1-((I+3)*2*M_PI/N)*((I+3)*2*M_PI/N)/2/3*(1-((I+3)*2*M_PI/N)*
           ((I+3)*2*M_PI/N)/4/5*(1-((I+3)*2*M_PI/N)*((I+3)*2*M_PI/N)/6/7*(1-((I+3)*2*M_PI/N)*
           ((I+3)*2*M_PI/N)/8/9*(1-((I+3)*2*M_PI/N)*((I+3)*2*M_PI/N)/10/11*(1-((I+3)*2*M_PI/N)*
           ((I+3)*2*M_PI/N)/12/13*(1-((I+3)*2*M_PI/N)*((I+3)*2*M_PI/N)/14/15*
           (1-((I+3)*2*M_PI/N)*((I+3)*2*M_PI/N)/16/17*
           (1-((I+3)*2*M_PI/N)*((I+3)*2*M_PI/N)/18/19*(1-((I+3)*2*M_PI/N)*
           ((I+3)*2*M_PI/N)/20/21)))))))))));
  }
};



template<unsigned N, unsigned I>
class AbsSineSSE {
public:
  static inline __m128 abssin() {
    return _mm_fabs_ps(SineSSE<N,I>::sin());
  }
};
#endif

#endif



