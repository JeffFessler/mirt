// jf,thread1.c
// generic routines for multi-threading
// for linux, now using "Portable Linux Processor Affinity (PLPA)"
// http://www.open-mpi.org/software/plpa/overview.php
// Copyright 2006-5-3, Jeff Fessler, University of Michigan


#define jf_nthread_big 160 // avoid alloc/free if nthread <= this value

#ifdef Use_thread

#include <stdarg.h> // for va_start in jf_thread_print() below

# ifdef Use_nptl // linux: native posix thread library
#  ifndef _GNU_SOURCE
#   define _GNU_SOURCE
#  endif
#  include <sched.h>
#  include <nptl/pthread.h>

# elif Use_plpa
#  include <plpa.h>
#  include <pthread.h>

# else // mac os
#  include <pthread.h>
# endif

#endif // Use_thread

#include "jf,thread1.h"

#ifdef Use_ncore_sysctl // old way for mac
#include <sys/sysctl.h>
#endif

#include <unistd.h> // for sysconf() to get ncore

#ifdef Use_aff_mac1 // affinity for mac leopard (and above?)
#include <mach/thread_policy.h>
#include <mach/mach.h>
#endif


#ifdef Use_thread

// jf_thread_print()
// thread-safe version of "fprintf" that keeps thread messages distinct
// called by Note() Warn() Fail()
void jf_thread_print(FILE *stream,
cchar *how, // usually "Note" or "WARN" or "FAIL"
cchar *file, // usually __FILE_
cint line, // usually __LINE__
cchar *format, ...)
{
	Mutex_init
	Mutex_lock
	(void) fprintf(stream, "%s %s %d: ", how, file, line);
	va_list args;
	va_start(args, format);
	(void) vfprintf(stream, format, args);
	(void) fprintf(stream, "\n");
	(void) fflush(stream);
	va_end(args);
	Mutex_unlock
}

#endif // Use_thread


// jf_thread1_node()
sof jf_thread1_node(
jf_thread1_node_mode mode, // jf_thread1_node_mode_set || jf_thread1_node_mode_get
cint inode, // mpi node index 0, 1, ..., nnode-1, used for jf_thread1_node_mode_set
int *p_node, // return address for previously set inode, used for jf_thread1_node_mode_get
cint chat)
{
	(void) chat;
	static int inode_static = 0;
	if (mode == jf_thread1_node_mode_set)
	{
		inode_static = inode;
	}
	else if (mode == jf_thread1_node_mode_get)
	{
		if (!p_node)
			Fail("p_node required for _get")
		*p_node = inode_static;
	}
		
	Ok
}


// jf_thread1_ncore()
// return number of available cores, if possible
// caution: returns 0 on failure
int jf_thread1_ncore(cint nwant)
{
	int ncore = 0;

#ifdef Use_ncore_sysctl // old way for mac
        int mib[2] = {CTL_HW, HW_NCPU};
        size_t len;
        len = sizeof(ncore);
        sysctl(mib, 2, &ncore, &len, NULL, 0);
	if (ncore <= 0)
		Fail1("sysctl returned %d", ncore)
#endif

#if 1 // this works on mac osx and linux
	ncore = sysconf( _SC_NPROCESSORS_ONLN );
	if (ncore == -1)
		Fail("sysconf() failed")
#endif

	if (nwant == -1) // user wants as many as possible
	{
		if (ncore)
			return ncore;
		else
			Fail("cannot determine # of cores")
	}

	else if (nwant == 0) // user would like many, but will accept 1
	{
		if (ncore)
			return ncore;
		else
		{
			Warn("cannot determine # of cores, defaulting to 1")
			return 1;
		}
	}

	else
	{
		if (ncore && ncore < nwant)
			Fail2("want %d cores but have only %d", nwant, ncore)
	}

        return ncore;
}


// pthread_attr_setaffinity_np()
// fake temporary routine for setting affinity.
// needed only if gcc compiler cannot find the real one.
// probably superceded by plpa library
#ifdef Provide_setaffinity
void pthread_attr_setaffinity_np(void)
{
	static int warned = 0;
	if (!warned)
	{
		Note("calling dummy setaffinity")
		warned = 1;
	}
}
#endif // Provide_setaffinity


#ifdef Use_aff_mac1

// jf_thread1_setaffinity_mac()
// for mac (leopard and above?) we set affinity after starting the thread
// but before doing any work!?
// apple thread affinity help, but does not mention pthread:
// http://developer.apple.com/releasenotes/Performance/RN-AffinityAPI/
// for example see:
// http://www.opensource.apple.com/darwinsource/projects/other/xnu-1228.3.13/tools/tests/affinity/sets.c
static sof jf_thread1_setaffinity_mac(cint ithread)
// Const jf_thread1_affinity *aff, // affinity control
{
	thread_extended_policy_data_t epolicy;
	epolicy.timeshare = FALSE;
	kern_return_t ret = thread_policy_set(
		mach_thread_self(), THREAD_EXTENDED_POLICY,
		(thread_policy_t) &epolicy, THREAD_EXTENDED_POLICY_COUNT);
	if (ret != KERN_SUCCESS)
		Fail1("thread_policy_set returned %d", ret)

	thread_affinity_policy_data_t apolicy;
	apolicy.affinity_tag = ithread + 1; // set affinity tag

	ret = thread_policy_set(
		mach_thread_self(), THREAD_EXTENDED_POLICY,
		(thread_policy_t) &apolicy, THREAD_EXTENDED_POLICY_COUNT);
	if (ret != KERN_SUCCESS)
		Fail1("thread_policy_set returned %d", ret)

	Ok
}

#endif // Use_aff_mac1


// jf_thread1_affinity_check()
sof jf_thread1_affinity_check(cint chat)
{
#if Use_nptl
	if (chat)
		Note("using nptl, so affinity should work")
#elif Use_plpa
	Call(PLPA_PROBE_OK == plpa_api_probe, ())
	if (chat)
		Note("using plpa, and affinity probe ok")
#elif Use_aff_mac1
	if (chat)
		Note("using mac affinity sets, which i hope works")
#else
	if (chat)
		Warn("affinity check called without support, disregarding")
#endif
	Ok
}


#ifdef Use_thread

// jf_thread1_setaffinity_attr()
// match thread to a given cpu
static sof jf_thread1_setaffinity_attr(
pthread_attr_t *attr,
cint ithread,
Const jf_thread1_affinity *aff, // affinity control
cint chat)
{
	if (!aff || aff->type == jf_thread1_affinity_none) Ok

#if Use_nptl || Use_plpa
	int affinity = ithread; // usual affinity
	if (aff->type == jf_thread1_affinity_mod && aff->nmod)
		affinity = ithread % aff->nmod;
	if (aff->type == jf_thread1_affinity_list && aff->list)
		affinity = aff->list[ithread];
#endif

#if Use_nptl
	{
	cpu_set_t cs;
	size_t cpu_set_size = sizeof(cs);
	__CPU_ZERO(&cs);
	__CPU_SET(affinity, &cs);
	pthread_attr_setaffinity_np(attr, cpu_set_size, &cs);
	if (chat) Note2("set affinity for thread %d to %d", ithread, affinity)
	}

#elif Use_plpa
	(void) attr;
	(void) chat;
	{
	plpa_cpu_set_t cs;
	size_t cpu_set_size = sizeof(cs);
	int ret;
	PLPA_CPU_ZERO(&cs);
	PLPA_CPU_SET(affinity, &cs);
	ret = plpa_sched_setaffinity(0, cpu_set_size, &cs);
	if (ret)
		Fail2("plpa_sched_setaffinity(affinity=%d) returned %d\n"
			"Perhaps you tried to use more threads than cores??",
			affinity, ret)
	}

#elif Use_aff_mac1
	// mac doesn't use attr to set affinity
	(void) attr;
	(void) chat;
	if (ithread == 0 && aff->type != jf_thread1_affinity_try)
		Fail("mac version supports only basic affinity support")

#else
	(void) attr;
	(void) chat;
	if (ithread == 0 && aff->type != jf_thread1_affinity_try)
		Warn("affinity support requested but not enabled!")
#endif
	Ok
}

#endif // Use_thread


// jf_thread1_glue()
// interface routine for threads
static void *jf_thread1_glue(void *in)
{
	jf_thread1_s *pt = (jf_thread1_s *) in;

// 2009-6-2 found that PRTS increases every time we call this!?
#ifdef Use_aff_mac1
	if (!jf_thread1_setaffinity_mac(pt->id))
	{
		pt->ok = sof_failure;
		return NULL;
	}
#endif // Use_aff_mac1
//	(void) jf_thread1_setaffinity_mac;

	pt->ok = (pt->init)(pt->ps, pt->id, pt->nthread);

//	pthread_exit((void*) in); // 2009-6-2 per llnl example
	return NULL; // ""
}


// jf_thread1_tops()
// top-level interface to threaded operations
// trick: only one of "ps" or "pps" should be used!
sof jf_thread1_tops(
jf_thread1_init_t fun_init, // required user function
jf_thread1_wrap_t fun_wrap, // optional user function
void *ps, // pointer to data structure used by threads
void **pps, // [nthread] pointers to structures ""
cint nthread, // # threads
Const jf_thread1_affinity *aff, // affinity control
cint chat)
{
	jf_thread1_s *pt;
	jf_thread1_s pt_pre[jf_nthread_big];

	if (nthread > jf_nthread_big)
	{
		Warn1("allocating space for %d threads", nthread)
		Mem0pure(pt, nthread)
	}
	else
		pt = pt_pre;

	if (ps && pps) Fail("only one of 'ps' and 'pps' may be non-null")
	if (!ps && !pps) Fail("one of 'ps' and 'pps' must be non-null")

	for (int it=0; it < nthread; ++it)
	{
		int inode = 0; // for mpi
		Call(jf_thread1_node, (jf_thread1_node_mode_get, 0, &inode, Chat))

		pt[it].init = fun_init;
		pt[it].ok = sof_failure;
		pt[it].id = it + inode * nthread; // 2012-09-27 for mpi
		pt[it].nthread = nthread;
		if (ps)
			pt[it].ps = ps; // all threads get same structure!
		else
			pt[it].ps = pps[it]; // each thread gets its own
	}

	if (nthread == 1) // to support non-threaded compiles
	{
		jf_thread1_glue(pt+0);
		if (!pt[0].ok)
			Fail("single thread failed")

		if (fun_wrap)
			Warn("fun_wrap unused for nthread=1")
	}

#ifdef Use_thread
	else
	{
		pthread_t *pid;
		pthread_t pid_pre[jf_nthread_big];

		if (nthread > jf_nthread_big)
			Mem0pure(pid, nthread)
		else
			pid = pid_pre;

		pthread_attr_t attr_, *p_attr = &attr_;
		pthread_attr_init(p_attr);
		pthread_attr_setdetachstate(p_attr, PTHREAD_CREATE_JOINABLE);

		for (int it=0; it < nthread; ++it)
		{
			// match thread to a given cpu, if requested
			Call(jf_thread1_setaffinity_attr, (p_attr, it, aff, chat))

			if (pthread_create(pid+it, p_attr, jf_thread1_glue,
				(void *) (pt+it)))
				Fail1("error creating thread %d", it)
		}

		if (pthread_attr_destroy(p_attr))
			Fail("pthread_attr_destroy()")

		for (int it=0; it < nthread; ++it)
		{
			if (pthread_join(pid[it], NULL))
				Fail1("pthread_join %d failed", it)
			if (!pt[it].ok)
				Fail1("thread %d failed", it)
		}

		if (fun_wrap)
			Call(fun_wrap, (pt, nthread))

		if (nthread > jf_nthread_big)
			Free0pure(pid)
	}
#else

	else
		Fail1("threads %d not done", nthread)

	(void) fun_wrap;
	(void) aff;
	(void) chat;
//	(void) jf_thread1_setaffinity_attr;

#endif

	if (nthread > jf_nthread_big)
		Free0pure(pt)

	Ok
}


// jf_thread1_top()
// simpler top-level interface to threaded operations
sof jf_thread1_top(
jf_thread1_init_t fun_init, // required user function
jf_thread1_wrap_t fun_wrap, // optional user function
void *ps, // pointer to data structure passed to threads
cint nthread, // # threads
cint chat)
{
	jf_thread1_affinity aff_, *aff = &aff_;
	aff->type = jf_thread1_affinity_try;
	Call(jf_thread1_tops,
		(fun_init, fun_wrap, ps, NULL, nthread, aff, Chat))
	Ok
}
