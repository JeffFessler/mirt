/*
* def,thread.h
*
* Must include pthread.h *before* my definitions, because pthread.h
* sets a __need_timespec flag which affects declarations in time.h
*
* Copyright 2001-8-22, Jeff Fessler, The University of Michigan
*/

/*
* this was my attempt to make it work on solaris,
* which does not seem to have the __need_timespec flag.
* i gave up since my dual-processor solaris machines are slow anyway
*/
#if 0
typedef struct timespec {
	long int dummy1;
	long int dummy2;
}
#endif

#include <pthread.h>

#define Thread2(func, t1, t2) \
	{ \
		t1.ok = False; \
		t2.ok = False; \
		{ \
		pthread_t thread_id1, thread_id2; \
		if (pthread_create(&thread_id1, NULL, func, (void*) &t1)) \
			Fail("error creating thread1") \
		if (pthread_create(&thread_id2, NULL, func, (void*) &t2)) \
			Fail("error creating thread2") \
		pthread_join(thread_id1, NULL); \
		pthread_join(thread_id2, NULL); \
		} \
		if (!t1.ok || !t2.ok) \
			Fail("thread failed") \
	}
