/*
* def,alloc.h
* Copyright 1997-4, Jeff Fessler, The University of Michigan
*/

#ifndef DefAlloc
#define DefAlloc

/* trick: thanks to Mathworks and NO_BOOL, use "int" here rather than bool */
extern void	*io_mem_alloc(cuint n, cuint s, cchar *, cint, cchar *);
extern bool	io_mem_free(void *p, cchar *, cint, cchar *);
extern bool	io_mem_print(cchar *, cint, cint);
extern bool	io_mem_usage(cchar *, cint, cint);
extern bool	io_mem_info(cvoid *, cchar *, cint, cchar *);


#endif /* DefAlloc */
