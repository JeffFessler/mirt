// defs-env.h 2015-07-24,16:56:18

// defs-top.h
// top-level portability definitions
//
// Copyright Dec. 1994, Jeff Fessler, University of Michigan

#ifndef DefsTop
#define DefsTop

// set the following flag for compiling on PC/NT/x86
#if defined(Is_pc)
#define	Userand
#define	Need_uint
#define Need_proto_lgamma
#define Need_proto_erf
#define Need_time_h
#endif

// Must include pthread.h *before* my other definitions, because pthread.h
// sets a __need_timespec flag which affects declarations in time.h
#ifdef Use_thread
// include "def,thread.h" used to go here
#include <pthread.h>
#endif

#include <stdio.h>
#include <string.h>
#if defined(Need_assert)
#include <assert.h>
#endif
#include <sys/types.h>
#include <ctype.h>
#if defined(__osf__) | defined(Need_time_h)
#include <time.h>
#endif
#if defined(Need_ieeefp)
#include <ieeefp.h> // needed for solaris for finite()
#endif
#if defined(Need_dynlib)
#include <dlfcn.h> // needed for dynamic libraries
#endif

#ifdef Use_mpi
#include <mpi.h>
#endif

#if defined(hphp)
#define _INCLUDE_XOPEN_SOURCE
#define Userand48
	extern double drand48(void);
#if !defined(Mmex)
#	define Need_uint
#endif
#endif

#if defined(Need_uint)
	typedef unsigned int	uint;
	typedef unsigned short	ushort;
#endif
#if defined(Need_ulong)
	typedef unsigned long	ulong;
#endif

#if defined(titan)
#define Provide_strstr // no strstr() on titan
#define Need_proto_strstr
#include <vmath.h>
#else
#include <math.h>
#endif

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif
#ifndef M_SQRT1_2
#define M_SQRT1_2 0.70710678118654752440
#endif
#ifndef M_SQRT2
#define M_SQRT2 1.4142135623730950488
#endif

#ifndef CLOCKS_PER_SEC // stupid SunOS doesn't have it
#define CLOCKS_PER_SEC 1000000
#endif

#ifndef Clock
#ifdef Use_global_clock0
#define	GlobalClock0declare	clock_t GlobalClock0;
#define	GlobalClock0extern	extern clock_t GlobalClock0;
#define	Clock0	GlobalClock0 = clock();
#define	Clock	( (clock() - GlobalClock0) / (double) CLOCKS_PER_SEC )
#else
#define	Clock0	;
#define	Clock	-1.
#define	GlobalClock0declare
#define	GlobalClock0extern
#endif
#endif

// No stdlib on titan!
#if defined(titan)
	extern void exit();
	extern char *calloc(), *malloc();
	extern char *memcpy(), *memset();
	extern int free();
#else
#	include <stdlib.h>
#endif

// Other externs
#ifdef Need_protos_sol2 // needed for sol2
#	define Need_proto_erf
#	define Need_proto_lgamma
#	define Need_proto_random
#endif

#if defined(titan) | defined(sun)
	extern int bzero(), bcopy();
#endif


#if !defined(Const)
#	define Const const
#endif

#if 0 // todo
#if !defined(noConst)
#	define Const const
#else
#	define Const
#endif
#endif


#if !defined(Size_t)
#	define Size_t size_t
#endif
typedef Const Size_t cSize_t;

#if 0 // todo
#if defined(_SIZE_T_)
#	define Size_t size_t
#else
#	define Size_t unsigned
#endif
#endif


#if !defined(MAXint)
#	define MAXchar		127
#	define MAXbyte		255
#	define MAXshort		32767
#	define MAXint		2147483647
#	define MAXint32		2147483647
#	define MAXuint32	4294967295
#	define MAXfloat		1e30
#	define MAXdouble	1e300
#endif

#ifndef SHRT_MAX
#define SHRT_MAX	32767
#endif
#ifndef SHRT_MIN
#define SHRT_MIN	(-32768)
#endif

// define all the globals that are needed by the macros
#ifndef GlobalsDefinedHere
#	define GlobalsDefinedHere \
		FILE *Global_fp_exist; \
		GlobalClock0declare \
		int Mpi_rank=0, Mpi_numprocs=1;
#endif

// my other definitions
#if 0
#include "def,macro.h"
#include "def,type.h"
#include "def,inline.h"
#include "def,alloc.h"
#include "def,proto.h"
#endif

#endif // DefsTop
// def,type.h
// type definitions
// Copyright 1995-2, Jeff Fessler, University of Michigan

#if !defined(Defs_type) && !defined(No_defs_type)
#define Defs_type

// uint will come from sys/types.h

// return type for functions that can fail or succeed:
#if 0
	typedef enum {sof_failure, sof_success} sof; // Failure | Success
#else // the enum version does not work with cuda 2.2 compiler, so just use int:
	typedef int sof; // Failure=0 | Success=1
#	define sof_failure 0
#	define sof_success 1
#endif

#define Success sof_success
#define Failure sof_failure
#define Ok return Success;

#if 0 // someday maybe use this version
typedef enum {truf_false, truf_true} truf; // False | True
#else
typedef int truf; // False=0 | True=1
#define truf_false 0
#define truf_true 1
#endif

#ifndef True
#define True 1
#endif
#ifndef False
#define False 0
#endif

typedef Const truf ctruf;


typedef short		word;	// two byte integers
typedef unsigned short	uword;	// two byte integers, unsigned
typedef unsigned char	byte;

typedef double TypeCalc; // fastest floating point type

typedef Const char	cchar;
typedef Const byte	cbyte;
typedef Const short	cshort;
typedef Const int	cint;
typedef Const long	clong;
typedef Const unsigned int	cuint;
typedef Const float	cfloat;
typedef Const double	cdouble;
typedef Const void	cvoid;

typedef cchar		*pcchar;
typedef const pcchar	*pcpcchar; // for argv
typedef cfloat		*pcfloat;

#endif // Defs_type
// def,macro.h
// favorite macro definitions
// Copyright 1994-12, Jeff Fessler, University of Michigan

#ifndef Defs_macro
#define Defs_macro

#define	Chat		((chat > 0) ? (chat-1) : 0)

#define	Calloc(n,s)	calloc((unsigned) (n), (unsigned) (s))

#define	Swab(s,d,n)	{ \
	if ((s) == (d)) Exit("swab problem") \
	swab((char *)(s), (char *)(d), n); }

#define	Nonull(str)	((str) ? (str) : "null")

#define Argflt(n,def)	(argc > n ? atof(argv[n]) : def)
#define Argint(n,def)	(argc > n ? atoi(argv[n]) : def)
#define Argstr(n,def)	(argc > n ? argv[n] : def)
#define Argstn(n)	((argc > n && strcmp(argv[n], "-")) ? argv[n] : NULL)

#define Isdigit(c)	(('0' <= (c)) && ((c) <= '9'))

#if !defined(Isfinite)
#define Isfinite(x)	isfinite(x) // from finite to isfinite 2008-10-30
#endif

#if !defined(Isinf)
#define Isinf(x)	isinf(x)
#endif


// messages: using variadic macro feature of c99

#ifdef Use_thread

// thread safe version (if posix threads) so thread messages do not mingle
extern void jf_thread_print( // in jf,thread1.c
FILE *stream, 
cchar *how, // usually "Note" or "WARN" or "FAIL"
cchar *file, // usually __FILE_
cint line, // usually __LINE__
cchar *format, ...);

#else // alternate general macro if no posix threads available

#define jf_thread_print(stream, how, file, line, ...) \
	{ \
	(void) fprintf(stream, "%s %s %d: ", how, file, line); \
	(void) fprintf(stream, __VA_ARGS__); \
	(void) fprintf(stream, "\n"); \
	(void) fflush(stream); \
	}

#endif

#ifndef Note
#	define Note(...) { \
	jf_thread_print(stdout, "Note", __FILE__, __LINE__, __VA_ARGS__); }

#define Note1(msg, a1)				Note(msg, a1)
#define Note2(msg, a1, a2)			Note(msg, a1, a2)
#define Note3(msg, a1, a2, a3)			Note(msg, a1, a2, a3)
#define Note4(msg, a1, a2, a3, a4)		Note(msg, a1, a2, a3, a4)
#define Note5(msg, a1, a2, a3, a4, a5)		Note(msg, a1, a2, a3, a4, a5)

#endif


#ifndef Warn
#	define Warn(...) { \
	jf_thread_print(stderr, "WARN", __FILE__, __LINE__, __VA_ARGS__); }

#define Warn1(msg, a1)				Warn(msg, a1)
#define Warn2(msg, a1, a2)			Warn(msg, a1, a2)
#define Warn3(msg, a1, a2, a3)			Warn(msg, a1, a2, a3)
#define Warn4(msg, a1, a2, a3, a4)		Warn(msg, a1, a2, a3, a4)
#endif


#ifndef Fail_nolock
#	define Fail_nolock(...) { \
	(void) fprintf(stderr, "FAIL %s %d: ", __FILE__, __LINE__); \
	(void) fprintf(stderr, __VA_ARGS__); \
	(void) fprintf(stderr, "\n"); \
	(void) fflush(stderr); \
	}
#endif

#ifndef Fail
#	define Fail(...) { \
	jf_thread_print(stderr, "FAIL", __FILE__, __LINE__, __VA_ARGS__); \
	return Failure; }

#define Fail1(msg, a1)				Fail(msg, a1)
#define Fail2(msg, a1, a2)			Fail(msg, a1, a2)
#define Fail3(msg, a1, a2, a3)			Fail(msg, a1, a2, a3)
#define Fail4(msg, a1, a2, a3, a4)		Fail(msg, a1, a2, a3, a4)
#define Fail5(msg, a1, a2, a3, a4, a5)		Fail(msg, a1, a2, a3, a4, a5)
#define Fail6(msg, a1, a2, a3, a4, a5, a6) Fail(msg, a1, a2, a3, a4, a5, a6)
#endif


#ifndef Exit
#	define Exit(msg)	{ \
	(void)fprintf(stderr, "EXIT %s %d: %s\n", __FILE__, __LINE__, msg); \
	exit(-1); }
#endif


#ifndef Call
#define Call(fun, arg)		{ if (!(fun arg)) Fail(#fun) }
#define Call1(fun, arg, str)	{ if (!(fun arg)) Fail1(#fun" %s", str) }
#endif

// And those pesky IO functions....
#if defined(titan)
#	define Fread(p,s,n,f)	fread ((char *)(p), (Size_t)(s), (int)(n), f)
#	define Fwrite(p,s,n,f)	fwrite((char *)(p), (Size_t)(s), (int)(n), f)
#elif defined(sun)
#	define Fread(p,s,n,f)	fread ((char *)(p), (int)(s), (int)(n), f)
#	define Fwrite(p,s,n,f)	fwrite((char *)(p), (int)(s), (int)(n), f)
#else
#	define Fread(p,s,n,f)	fread ((void *)(p), (Size_t)(s), (Size_t)(n), f)
#	define Fwrite(p,s,n,f)	fwrite((cvoid *)(p), (Size_t)(s), (Size_t)(n), f)
#	define Strncat(s1,s2,n)	strncat(s1, s2, (int) (n))
#endif

#if defined(titan) | defined(sun)
#define Memset(p,c,n)		memset((char *)(p), (int)(c), (int)(n));
#define Memcpy(t,f,n)		memcpy((char *)(t), (char *)(f), (int)(n));
#define Strncat(s1,s2,n)	strncat(s1, s2, (Size_t) (n))
#define Bzero(p,n)	bzero((char *) (p), (int) ((n) * sizeof(*(p))));
#define Bcopy(s,d,n)	bcopy((char *) (s), (char *) (d), (int) ((n) * sizeof(*(s))));
#else
#define Memset(p,v,n)	(void)memset((void *) (p), v, (Size_t) (n));
#define Memcpy(d,s,n)	{ if ((d) == (s)) Exit("Identical memcpy pointers") \
		else (void)memcpy((void *) (d), (cvoid *) (s), (Size_t) (n)); }
#define Bzero(p,n)	Memset(p, 0, (n) * sizeof(*(p)))
#define Bcopy(s,d,n)	Memcpy(d, s, (n) * sizeof(*(s)))
#define Bcopyvoid(src,dst,n,siz)	Memcpy(dst, src, (n) * (siz))
#endif

#if 0 // superceded by versions in io library
#ifndef FileExist
	// (!fclose(fopen(name, "rb"))) did not work on linux
	extern FILE *Global_fp_exist;
#define FileExist(name)		((Global_fp_exist = fopen(name, "rb"), \
				Global_fp_exist && !fclose(Global_fp_exist)))
#endif

// macro to check to see if we can write a file to this name.
// use this with caution since if the file exists, it will be truncated.
#ifndef FileWriteable
#define FileWriteable(s)	((Global_fp_exist = fopen(s, "wb"), \
				Global_fp_exist && fclose(Global_fp_exist)))
#endif
#endif

// The suffix "0" means return 0 on error.
#define	Fopen0(fp,n,t)		{if ( !(fp = fopen(n, t)) )	\
					Fail1("fopen '%s'", n)}
#define	Fopen0r(fp,n)		Fopen0(fp, n, "rb")
#define	Fopen0w(fp,n)		Fopen0(fp, n, "wb")
#define	Fopen0e(fp,n)		Fopen0(fp, n, "rb+") // for editing
#define	Fflush0(fp)		{if (fflush(fp))	Fail("fflush")}
// seek from start of file
#define	Fseek0(fp,n)		{if (fseek(fp, (long) (n), 0))	\
					Fail1("fseek %d", (int) n)}
#define	Fskip0(fp,n)		{if (fseek(fp, (long) (n), 1))	\
					Fail1("fseek1 %d", (int) n)}
#define	Fclose0(fp)		{if (fclose(fp))	Fail("fclose")}
#define	Fread0(p,s,n,fp)	{Size_t ZF = Fread(p, s, n, fp); \
				if (ZF != (Size_t) (n)) \
				Fail3("fread found %zu bytes, not %zu at %ld", \
					ZF, (Size_t) (n), ftell(fp))}
#define	Fwrite0(p,s,n,fp)	{if (((Size_t) (n)) != Fwrite(p, s, n, fp)) Fail("fwrite")}


#if defined(Userand48)
#	define Rand01		drand48()
#	define Seeder(s)	srand48(s)
#elif defined(Userand) // for PC
#	define Rand01		(rand() / RAND_MAX)
#	define Seeder(s)	srand(s)
#else
#	define Rand01		(random() / 2147483647.0)
#	define Seeder(s)	srandom(s)
#endif

// The remaining macro definitions should be environment independent

#define Sqr(x)		((x) * (x))
#define Odd(n)		((n)/2 != ((n)+1)/2)
#define Even(n)		((n)/2 == ((n)+1)/2)
#ifndef Abs
#	define Abs(x)	(((x) >= 0) ? (x) : (-(x)))
#endif
#ifndef Max
#	define Max(x,y)	(((x) > (y)) ? (x) : (y))
#endif
#ifndef Min
#	define Min(x,y)	(((x) < (y)) ? (x) : (y))
#endif

// Kernighan & Ritchie, p45 2nd ed., say float to int just truncates.
// So for nonnegative arguments we can just cast to an integer
// to avoid a function call.  But Floor_pos will be wrong if x < 0!
#define Floor_pos(x)	((int) (x))
// the usual floor function works for both positive and negative values
#define Floor(x)	floor(x)
#define Floorf(x)	floorf(x)
#define Ceil(x)		ceil(x)
#define Ceilf(x)	ceilf(x)

#ifndef Rint // round to nearest integer, ala rint()
#	define Rint(x)	( (int) ((x) >= 0 ? ((x) + 0.5) : ((x) - 0.5)) )
	// formerly ((int) ((x) + 0.5)), which is wrong for negatives!
#endif

#define Streq(a,b)	(((a) != NULL) && ((b) != NULL) && !strcmp(a,b))
#define Streqn(a,b,n)	(((a) != NULL) && ((b) != NULL) && !strncmp(a,b,n))

#define	Alloc0pure(ptr,type,n,s) { \
	if ( !((ptr) = (type *) calloc((Size_t) (n), (Size_t) (s))) ) \
		Fail2("alloc %zu %zu", (Size_t) (n), (Size_t) (s)) }

#define Free0pure(p)	free((void *) (p));

#ifdef CountAlloc

#define	Alloc0(ptr,type,n,s)	{ \
	if ( !((ptr) = (type *) io_mem_alloc((Size_t) (n), (Size_t) (s), \
		__FILE__, __LINE__, #ptr)) ) \
		Fail2("io_mem_alloc(%zu,%zu)", (Size_t) (n), (Size_t) (s)) }
#define Free(p)		{ \
		if (!io_mem_free((void *) (p), __FILE__, __LINE__, #p)) { \
			Warn2("bug free %s %d", __FILE__, __LINE__) \
			Exit("bug") } }
#define Free0(p)	{ \
		Call(io_mem_free, ((void *) (p), __FILE__, __LINE__, #p)) \
		p = NULL; }

#define PrintAlloc	io_mem_print(__FILE__, __LINE__, 1);
#define MemUsed		io_mem_usage(__FILE__, __LINE__, Chat);
#define MemInfo(p)	io_mem_info((cvoid *) p, __FILE__, __LINE__, #p);
#define MemChat(chat)	Call(io_mem_chat, (chat))

#else

#define	Alloc0		Alloc0pure
#define Free(p)		free((void *) (p))
#define Free0		Free0pure

#define PrintAlloc	{}
#define MemUsed		{}
#define MemInfo(p)	{}
#define MemChat(chat)	{}

#endif

#ifdef __CUDACC__
#ifndef jf_use_typeof
#define jf_use_typeof
#endif
#endif // __CUDACC__

#ifdef jf_use_typeof // gcc provides typeof() extension (sometimes)
#define jf_typeof(p) typeof(p)
#else
#define jf_typeof(p) void
#endif

#ifndef Mem0
#define Mem0(p,n)	Alloc0(p, jf_typeof(*(p)), n, sizeof(*(p)))
#define Mem0void(p,n,s)	Alloc0(p, jf_typeof(*(p)), n, s)
#define Mem0pure(p,n)	Alloc0pure(p, jf_typeof(*(p)), n, sizeof(*(p)))
#endif

// macros for making code sections with static variable work for threads
// Mutex_init declares a *static* mutex that can be used repeatedly
// Mutex_destroy need not be called except possibly at program exit
#if defined(Use_thread)
#	define Mutex_init \
		static pthread_mutex_t jf_mutex_reent = PTHREAD_MUTEX_INITIALIZER;
#	define Mutex_lock \
		if (pthread_mutex_lock( &jf_mutex_reent )) \
			Fail_nolock("pthread_mutex_lock() failed")
#	define Mutex_unlock \
		if (pthread_mutex_unlock( &jf_mutex_reent )) \
			Fail_nolock("pthread_mutex_unlock() failed")
#	define Mutex_destroy \
		if (pthread_mutex_destroy( &jf_mutex_reent )) \
			Fail_nolock("pthread_mutex_destroy() failed")
#	define Ok_unlock Mutex_unlock Ok
#else
#	define Mutex_init
#	define Mutex_lock
#	define Mutex_unlock
#	define Mutex_destroy
#	define Ok_unlock Ok
#endif

#endif // Defs_macro
// def,inline.h
// Simple vector operations done in-line.
// Copyright 1994-12, Jeff Fessler, University of Michigan

#ifndef DefInline
#define DefInline

// #define whileZn(start)	register int Zn = start; for (; Zn; --Zn)
#define whileZn(start)	for (size_t Zn = start; Zn-- > 0; )

// #define Register	register
#define Register	{}

// v op c(onstant)
#define Inline1vector(type1, v1, op, type2, c, n) \
{				\
Register type1 *Z1 = v1;	\
Register Const type2 Zc = (type2) c;	\
	whileZn(n)		\
		*Z1++ op Zc;	\
}

// v1 op v2
#define Inline2vector(type1, v1, op, type2, v2, n) \
{				\
Register type1 *Z1 = v1;	\
Register Const type2 *Z2 = v2;	\
	whileZn(n) {		\
		*Z1++ op (*Z2);	\
	++Z2; }			\
}




// Vector with scalar

// v *= c(onstant)
#define VectScale(type, v, c, n) \
	Inline1vector(type, v, *=, double, c, n)

// v += c(onstant)
#define VectInc(type, v, c, n) \
	Inline1vector(type, v, +=, double, c, n)

// v = c(onstant)
#define VectSet(type, v, c, n) \
	Inline1vector(type, v, =, type, c, n)

// v = Min(v,c)
#define VectMin(type, v, c, n) \
{				\
Register type *Zv = v;		\
	whileZn(n) {		\
		if (*Zv > c) *Zv = c; \
		++Zv; }		\
}


// v = Max(v,c)
#define VectMax(type, v, c, n) \
{			\
Register type *Zv = v;	\
	whileZn(n) {	\
		if (*Zv < c) *Zv = c; \
		++Zv; } \
}

// v = Max(v,0)
#define VectNonneg(type, v, n)	VectMax(type, v, 0, n)

// clamp vector to [min max]
#define VectClamp(type, v, vmin, vmax, n) \
{					\
Register type *Zv = v;			\
Register Const type Zmin = vmin;	\
Register Const type Zmax = vmax;	\
	whileZn(n) {				\
		if (*Zv > Zmax) *Zv = Zmax;	\
		else				\
		if (*Zv < Zmin) *Zv = Zmin;	\
		++Zv;				\
	}					\
}


// s = sum(v)
#define VectAccum(type, s, v, n) \
{					\
Register double Zsum = 0.;		\
Register Const type *Zv = v;		\
	whileZn(n)			\
		Zsum += *(Zv++);	\
	s = Zsum;			\
}


// s = sum(|v|)
#define VectNorm1(type, s, v, n) \
{					\
Register Const type *Zv = v;		\
Register double Zsum = 0;		\
	whileZn(n) {			\
		Zsum += Abs(*Zv);	\
		++Zv;			\
	} s = Zsum;			\
}


// s = <v,v>
#define VectNorm2(type, s, v, n) \
{					\
Register Const type *Zv = v;	\
Register double Zsum = 0;		\
	whileZn(n) {			\
		Zsum += *Zv * *Zv;	\
		++Zv;			\
	} s = Zsum;			\
}


// s = <v,v>_w
#define VectNorm2wtd(type, s, v, w, n) \
{					\
Register Const type *Zv = v;	\
Register Const type *Zw = w;	\
Register double Zsum = 0;			\
	whileZn(n) {				\
		Zsum += *Zv * *Zv * *Zw;	\
		++Zv; ++Zw;			\
	} s = Zsum;				\
}


// [v; v; ...; v] repeated nrep times
#define VectRepeat(type, v, npoint, nrep) \
{				\
Register type *Zv = v;		\
Const size_t Zn = npoint;	\
Const size_t Zr = nrep;		\
	for (size_t Zi=1; Zi < Zr; Zi++) \
		Bcopy(Zv, Zv+Zi*Zn, Zn)	\
}


// Vector with Vector

// v1 = v2
#define VectCopy(type1, v1, type2, v2, n) \
	Inline2vector(type1, v1, =, type2, v2, n)

// r *= v
#define VectMul(type, r, v, n) \
	Inline2vector(type, r, *=, type, v, n)

// r /= v
#define VectDiv(type, r, v, n) \
	Inline2vector(type, r, /=, type, v, n)

// r += v
#define VectAdd(type, r, v, n) \
	Inline2vector(type, r, +=, type, v, n)

// r += v
#define VectAdd2(type1, r, type2, v, n) \
	Inline2vector(type1, r, +=, type2, v, n)

// r += const * v
#define VectAddScale(type, r, c, v, n) \
	Inline2vector(type, r, += c *, type, v, n)

// r = s * v
#define VectCopyScale(type, r, s, v, n) \
	Inline2vector(type, r, = s *, type, v, n)

// r -= v
#define VectSub(type, r, v, n) \
	Inline2vector(type, r, -=, type, v, n)

// r += v1 + v2
#define VectInc2(type, r, v1, v2, n) \
{					\
Register type	*Zr = r;		\
Register Const type *Z1 = v1;		\
Register Const type *Z2 = v2;		\
	whileZn(n)			\
		*Zr++ += *(Z1++) + *(Z2++);	\
}

// r = r / v or 0 if v=0
#define VectDiv0(type, r, v, n) \
{				\
Register type *Zr = r;		\
Register Const type *Zv = v;	\
	whileZn(n) {		\
		if (*Zv)	\
			*Zr++ /= *Zv;	\
		else {			\
			*Zr++ = 0;	\
		}			\
		++Zv; }		\
}


// s = <v1, v2>
#define VectInprod3d(s, type1, v1, type2, v2, n1, n2, n3) \
{					\
Register Const type1 *Z1 = v1;		\
Register Const type2 *Z2 = v2;		\
Const size_t Zn2 = n2;		\
Register Const size_t Zn1 = n1;	\
double Zsum = 0.;		\
	whileZn(n3) {		\
		double Zsum2 = 0.;	\
		for (size_t Zi2=0; Zi2 < Zn2; ++Zi2) {		\
			Register double Zsum1 = 0.;		\
			for (size_t Zi1=0; Zi1 < Zn1; ++Zi1)	\
				Zsum1 += *(Z1++) * *(Z2++);	\
			Zsum2 += Zsum1;				\
		}		\
		Zsum += Zsum2;	\
	}			\
	s = Zsum;		\
}

#define VectInprod2d(s, type1, v1, type2, v2, n1, n2) \
	VectInprod3d(s, type1, v1, type2, v2, n1, n2, 1)

// s = <v1, v2>
#define VectInprod(type, s, v1, v2, n) \
{					\
Register Const type *Z1 = v1;		\
Register Const type *Z2 = v2;		\
Register double Zsum = 0.;		\
	whileZn(n)			\
		Zsum += *(Z1++) * *(Z2++);	\
	s = Zsum;			\
}


// s = <v1, v2.^2>
#define VectInprod2(type, s, v1, v2, n) \
{					\
Register Const type *Z1 = v1;		\
Register Const type *Z2 = v2;		\
Register double Zsum = 0.;		\
	whileZn(n) {			\
		Zsum += *(Z1++) * *Z2 * *Z2;	\
		++Z2;			\
	} s = Zsum;			\
}


// s = <v1, v2.^2>
// Supposes v2 is mostly zeros
#define VectInprod2_0(type, s, v1, v2, n) \
{					\
Register Const type *Z1 = v1;		\
Register Const type *Z2 = v2;		\
Register double Zsum = 0.;		\
	whileZn(n) {			\
		if (*Z2)		\
			Zsum += *Z1 * *Z2 * *Z2;	\
		++Z1; ++Z2;		\
	} s = Zsum;			\
}


// min and max of vector
#define VectMinMax(type1, vmin, vmax, type2, v, n) \
{						\
Register Const type2 *Zv = v;		\
Register type1 Zmin = (type1) *Zv;	\
Register type1 Zmax = (type1) *Zv++;	\
	whileZn(n-1) {				\
		if (*Zv > Zmax)	Zmax = *Zv;	\
		else				\
		if (*Zv < Zmin)	Zmin = *Zv;	\
		++Zv;				\
	}					\
	vmin = Zmin;				\
	vmax = Zmax;				\
}


// min and max of vector, and min > 0
#define VectMinMax0(type1, this_type_max, vmin, vmin0, vmax, type2, v, n) \
{						\
Register Const type2 *Zv = v;			\
Register type1 Zmin = (type1) *Zv;		\
Register type1 Zmin0 = (type1) (Zmin > 0 ? Zmin : this_type_max);	\
Register type1 Zmax = (type1) *Zv++;		\
	whileZn(n-1) {				\
		if (*Zv > Zmax)	Zmax = *Zv;	\
		else				\
		if (*Zv < Zmin)	Zmin = *Zv;	\
		if (*Zv > 0 && *Zv < Zmin0)	Zmin0 = *Zv;	\
		++Zv;				\
	}					\
	vmin = Zmin;				\
	vmin0 = Zmin0;				\
	vmax = Zmax;				\
}


// for debugging

#define VectInfo(type, v, n, msg) { \
	double Ztmin, Ztmax, Ztmin0, ZZsum; \
	VectMinMax0(type, MAX##type, Ztmin, Ztmin0, Ztmax, type, v, n) \
	VectAccum(type, ZZsum, v, n) \
	Note("%s [%d] min=%g min0=%g max=%g sum=%g", \
		msg, n, Ztmin, Ztmin0, Ztmax, ZZsum) }

// #undef Register
// #undef whileZn

#endif // DefInline
// def,alloc.h
// Copyright 1997-4, Jeff Fessler, University of Michigan

#ifndef DefAlloc
#define DefAlloc

extern void *io_mem_alloc(cSize_t n, cSize_t s, cchar *file, cint line, cchar *name);
extern sof io_mem_free(void *p, cchar *file, cint line, cchar *name);
extern sof io_mem_print(cchar *file, cint line, cint chat);
extern sof io_mem_usage(cchar *file, cint line, cint chat);
extern sof io_mem_info(cvoid *p, cchar *file, cint line, cchar *name);
extern sof io_mem_chat(cint chat);

#endif // DefAlloc
// def,proto.h
// Prototypes not provided by some old OS versions
//
// Copyright 2000-5, Jeff Fessler, University of Michigan

#ifndef DefProto
#define DefProto

#ifdef Need_proto_strstr
	extern char *strstr(char *, char *);
#endif

#ifdef Need_proto_random
	extern void srandom(unsigned);
	extern long random(void);
#endif

#ifdef Need_proto_lgamma
	extern double lgamma(cdouble);
#endif

#ifdef Need_proto_erf
	extern double erf(cdouble);
#endif

#endif // DefProto
// def,simd.h
// macros for SIMD instructions (SSE)
// based on src/contrib/ndimtro/SSE.h by Ned Dmitro
// Copyright 2010-03-08, Jeff Fessler, University of Michigan

#ifndef _def_simd_h
#define _def_simd_h

// Notes
// When macro arguments are preceded by a #, the compiler (cpp)
// "stringifies" them into " $Value ".
// It also concatenates any adjacent string literals.

// ATT syntax, unlike Intel syntax is COMMAND SOURCE DEST
// assembler immediate values must be preceded by $


// Put the 32-bit data "val" (e.g., float) in all 4 spots in XMM register NUM
// movss moves single precision floating point to lowest word of XMM register
// shufps selects parts of source and dest to move to dest

#define simd_fill4(val, NUM) { \
	__asm__( \
	"movss %0, %%xmm" #NUM \
	"; shufps $0x00, %%xmm" #NUM ", %%xmm" #NUM ";" : : "m" (val) ); \
	}


// The following either move data Aligned (a) or Unaligned (u)
// either "in" to XMM register NUM from memory, or "out"
// val is the starting value of the 4-float array used, e.g., data[0]

#define simd_mov_a4_in(val, NUM) { \
	__asm__("movaps %0, %%xmm" #NUM ";" : : "m" (val) ); \
	}

#define simd_mov_a4_out(NUM, val) { \
	__asm__("movaps %%xmm" #NUM ", %0;" : "=m" (val) : ); \
	}

#define simd_mov_u4_in(val, NUM) { \
	__asm__("movups %0, %%xmm" #NUM ";" : : "m" (val) ); \
	}

#define simd_mov_u4_out(NUM, val) { \
	__asm__("movups %%xmm" #NUM ", %0;" : "=m" (val) : ); \
	}


// Perform an operation on XMM NUM1 and XMM NUM2 and leave result in XMM NUM 2
// XMM2 = XMM2 op XMM1

#define simd_add4(NUM1, NUM2) { \
	__asm__("addps %xmm" #NUM1 ", %xmm" #NUM2 ";" ); \
	}

#define simd_mul4(NUM1, NUM2) { \
	__asm__("mulps %xmm" #NUM1 ", %xmm" #NUM2 ";" ); \
	}

#define simd_div4(NUM1, NUM2) { \
	__asm__("divps %xmm" #NUM1 ", %xmm" #NUM2 ";" ); \
	}

#define simd_sub4(NUM1, NUM2) { \
	__asm__("subps %xmm" #NUM1 ", %xmm" #NUM2 ";" ); \
	}


#if 0 // hide

#define shufSSE(immediate, source, dest) {\
	__asm__("shufps $0b" #immediate ", %xmm" #source ", %xmm" #dest ";" );\
	}

#define cvttSSE(source, dest) {\
	__asm__("cvttps2pi %%xmm" #source ", %%mm" #dest ";" : : : "st", "st(1)", "st(2)", "st(3)", "st(4)", "st(5)", "st(6)", "st(7)");\
	}

#define cvtSSE(source, dest) {\
	__asm__("cvtps2pi %%xmm" #source ", %%mm" #dest ";" : : : "st", "st(1)", "st(2)", "st(3)", "st(4)", "st(5)", "st(6)", "st(7)");\
	}

// The following preform the operation on XMM NUM1 and XMM NUM2 and
// place the result in XMM NUM 2 ie: XMM2 = XMM2 op XMM1

#define andSSE(NUM1, NUM2) {\
	__asm__("andps %xmm" #NUM1 ", %xmm" #NUM2 ";" ); \
	}

#define andnSSE(NUM1, NUM2) {\
	__asm__("andnps %xmm" #NUM1 ", %xmm" #NUM2 ";" ); \
	}

#define orSSE(NUM1, NUM2) {\
	__asm__("orps %xmm" #NUM1 ", %xmm" #NUM2 ";" ); \
	}

#define xorSSE(NUM1, NUM2) {\
	__asm__("xorps %xmm" #NUM1 ", %xmm" #NUM2 ";" ); \
	}

/* The following prefetches the line of data that contains the byte
 val into non-temporal cache storage
*/
#define prefetchSSE(val) { \
	__asm__("prefetchnta %0;" : : "m" (val) );\
	}

#define movqOutMMX(NUM, dest) {\
	__asm__("movq %%mm" #NUM ", %0;" : "=m" (dest) : );\
	}

#define movqInMMX(source, NUM) {\
	__asm__("movq %0, %%mm" #NUM ";" : : "m" (source) : "st(" #NUM ")" );\
	}

#define adddMMX(NUM1, NUM2) {\
	__asm__("paddd %mm" #NUM1 ", %mm" #NUM2 ";");\
	}

#define emmsMMX() {\
	__asm__("emms;");\
	}

#endif // hide

#endif // _def_simd_h
