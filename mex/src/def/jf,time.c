// jf,time.c
// measure wall times
// Copyright 2006-5-1, Jeff Fessler, University of Michigan

#include "jf,time.h"
// #include <sys/time.h>


#if 0
// jf_time_alloc()
jf_timeval *jf_time_alloc(void)
{
	jf_timeval *ptr;
	Mem0(ptr, 1)
	return ptr;
}
#endif


// jf_time_tic()
// start clock
sof jf_time_tic(jf_timeval *ptr)
{
	if (gettimeofday(ptr, NULL /* struct timezone */))
		Fail("gettimeofday() failed")
	Ok
}


// jf_time_toc()
sof jf_time_toc(Const jf_timeval *ptr, cchar *arg)
{
	struct timeval t2;
	if (gettimeofday(&t2, NULL /* struct timezone */))
		Fail("gettimeofday() failed")
	time_t ds = t2.tv_sec - ptr->tv_sec;
	int du = t2.tv_usec - ptr->tv_usec;
	Note2("wall time %s: %g", arg, ds + du * 1e-6)
	Ok
}


#if 0
// jf_time_free()
sof jf_time_free(jf_timeval *ptr)
{
	Free0(ptr)
	Ok
}
#endif


// jf_time()
sof jf_time(cchar *arg)
{
	static struct timeval t1;
	static truf t1_set = False;

	if (Streq(arg, "start"))
	{
		if (gettimeofday(&t1, NULL /* struct timezone */))
			Fail("gettimeofday() failed")
		t1_set = True;
	}

	else if (Streqn(arg, "report", 6))
	{
		if (!t1_set)
			Fail("bug: must first 'start'")
		struct timeval t2;
		if (gettimeofday(&t2, NULL /* struct timezone */))
			Fail("gettimeofday() failed")
#if 0
		Note4("t2 = %d %d, t1 = %d %d",
			(int) t2.tv_sec, (int) t1.tv_sec,
			(int) t2.tv_usec, (int) t1.tv_usec)
#endif
		time_t ds = t2.tv_sec - t1.tv_sec;
		int du = t2.tv_usec - t1.tv_usec;
//		Note3("wall time: %d,%d = %g", (int) ds, (int) du, ds + du * 1e-6)
		Note2("wall time%s: %g", arg+6, ds + du * 1e-6)
	}

	else
		Fail1("unknown arg %s", arg)

	Ok
}
