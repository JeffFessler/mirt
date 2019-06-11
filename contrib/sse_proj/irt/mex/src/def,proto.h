/*
* def,proto.h
* Prototypes
*
* Copyright 2000-5, Jeff Fessler, The University of Michigan
*/

#ifndef DefProto
#define DefProto

#ifdef Need_proto_strstr
	extern char *strstr(char *, char *);
#endif

#ifdef Need_proto_random
	extern void	srandom(unsigned);
	extern long	random(void);
#endif

#ifdef Need_proto_lgamma
	extern double lgamma(cdouble);
#endif

#ifdef Need_proto_erf
	extern double erf(cdouble);
#endif

#endif	/* DefProto */
