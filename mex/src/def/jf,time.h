// jf,time.h
// Copyright 2009-5-22, Jeff Fessler, University of Michigan

#ifndef jf_time_h
#define jf_time_h

#include "defs-env.h"

#include <sys/time.h>

typedef struct timeval jf_timeval;

// extern jf_timeval *jf_time_alloc(void);
extern sof jf_time_tic(jf_timeval *);
extern sof jf_time_toc(Const jf_timeval *ptr, cchar *arg);
// extern sof jf_time_free(jf_timeval *ptr);

extern sof jf_time(cchar *arg);

#endif // jf_time_h
