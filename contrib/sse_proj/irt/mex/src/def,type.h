/*
* def,type.h
* type definitions
*
* Copyright 1995-2, Jeff Fessler, The University of Michigan
*/

#ifndef Defs_type
#define Defs_type

/* uint will come from sys/types.h */
#if !defined(NO_BOOL) && !defined(mat_h)	/* thank you MathWorks */
typedef int		bool;	/* True or False, Success Or Failure	*/
#endif
typedef short		word;	/* two byte integers			*/
typedef unsigned short	uword;	/* two byte integers, unsigned		*/
typedef unsigned char	byte;
typedef double	TypeCalc;	/* fastest floating point type		*/

typedef Const char	cchar;
typedef Const byte	cbyte;
typedef Const short	cshort;
typedef Const int	cint;
typedef Const unsigned int	cuint;
typedef Const float	cfloat;
typedef Const double	cdouble;
typedef Const void	cvoid;
typedef Const int	cbool;

typedef cchar		*pcchar;
typedef const pcchar	*pcpcchar;	/* for argv */
typedef cfloat		*pcfloat;

#endif /* Defs_type */
