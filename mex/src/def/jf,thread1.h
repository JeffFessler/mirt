// jf,thread1.h
// generic declarations for jf multi-threading interface.
// these do not rely on any OS specific calls so they are portable.
// Copyright 2006-5-3, Jeff Fessler, University of Michigan

#ifndef jfDefThread
#define jfDefThread

#include "defs-env.h"

// user init routine type.
// user provides this routine to be called to start each thread
typedef sof jf_thread1_init_t(void *, cint id, cint nthread);

// single structure passed to pthread_create()
typedef struct {
	jf_thread1_init_t *init; // init function pointer
	void *ps; // pointer to structure with all data etc.
	int id; // 0:(nthread-1)
	int nthread; // # threads
	sof ok; // 1 for success, 0 for failure
} jf_thread1_s;

// wrap-up routine type
// user provides this optional routine for clean up after all threads exit.
typedef sof jf_thread1_wrap_t(jf_thread1_s *, cint nthread);

// control thread affinity (if possible)
typedef enum {
	jf_thread1_affinity_none, // do not set affinity
	jf_thread1_affinity_try, // use: ithread (if possible)
	jf_thread1_affinity_mod, // use: ithread % nmod
	jf_thread1_affinity_list, // use: list[ithread]
} jf_thread1_affinity_type;

// help for mpi 
typedef enum {
	jf_thread1_node_mode_dummy,
	jf_thread1_node_mode_set,
	jf_thread1_node_mode_get,
} jf_thread1_node_mode;



// thread affinity control arguments
typedef struct {
	jf_thread1_affinity_type type;
	int nmod;
	int *list;
} jf_thread1_affinity;


// jf,thread1.c

extern sof jf_thread1_node(jf_thread1_node_mode mode, cint inode, int *p_node, cint chat);

extern int jf_thread1_ncore(cint nwant);

// user calls this routine to check if affinity control available
extern sof jf_thread1_affinity_check(cint chat);

// user calls this top routine to start a threaded operation
extern sof jf_thread1_tops(
jf_thread1_init_t, // required user function
jf_thread1_wrap_t, // optional user function
void *ps, // pointer to structure with relevant data needed by threads
void **pps, // [nthread] pointers to structures ""
cint nthread, // # threads
Const jf_thread1_affinity *paff, // affinity control
cint chat);

// simple interface
extern sof jf_thread1_top(
jf_thread1_init_t, // required user function
jf_thread1_wrap_t, // optional user function
void *ps, // pointer to structure with relevant data needed by threads
cint nthread,
cint chat);

// jf,ranges.c
extern sof jf_thread1_ranges(
int *p_start,
int *p_inc,
int *p_end, // [start, end)
cint i_start,
cint i_inc,
cint i_end, // [start, end)
cint ithread,
cint nthread,
cint inode,
cint nnode,
cint chat);

extern sof jf_thread1_mask_init(
byte *mask_int, // [nn] out: integers 0 or 1 ... nthread
cbyte *mask_bin, // [nn] in: binary 0|1, possibly null
cint nn, // usually nx * ny
cint nthread,
ctruf chat);

#endif // jfDefThread
