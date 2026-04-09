/*
* cbct,def.h
* cone-beam CT definitions
* Copyright 2008-10-09, Jeff Fessler, University of Michigan
*/
#ifndef jf_cbct_def_h
#define jf_cbct_def_h

#include "defs-env.h"

// initial version used "double" for scalars
#ifndef rscalar
#define rscalar float
#endif
#define crscalar Const rscalar

#define Isinf(x)	isinf(x)

// image geometry
typedef struct {
	int	nx;	// image dimensions
	int	ny;
	int	nz;
	rscalar dx;	// pixel size (can be negative)
	rscalar dy;		// can be negative to cause flip
	rscalar dz;		// cannot be negative
	rscalar offset_x;	// center offset in pixels (usually 0)
	rscalar offset_y;
	rscalar offset_z;
	byte	*mask2;		// [nx ny] 2D support mask
				// 0 or 1 ... nthread
} cbct_ig;

// cone-beam geometry
typedef struct {
	rscalar dso;	// distance from source to isocenter
	rscalar dsd;	// distance from source to detector
	rscalar dfs;	// distance from detector focal point to source
			// 0 for 3rd-gen CT, infinity for flat detector
	int	ns;	// # detector channels per row
	int	nt;	// # detector rows (along axial direction)
	rscalar ds;	// horizontal (transaxial) ray spacing
	rscalar dt;	// vertical (axial) ray spacing
	rscalar offset_s;	// channel offset [pixels]
	rscalar offset_t;	// vertical offset on detector [pixels]
} cbct_cg;

// work space (one for each thread)
typedef struct {
	float *view; // ns*nt
} cbct_work;


// cbct,mask2.c
extern jool cbct_mask_init(
byte *mask_int,	// [nx ny] out: integers 0 or 1 ... nthread
cbyte *mask_bin, // [nx ny] in: binary, possibly null
cint nx, cint ny,
cint nthread,
cjool chat);

typedef enum {
	cbct_back_error,
	cbct_back_zero, // zero before back-projection
	cbct_back_inc, // incrementing back-projection (no zeroing first)
} cbct_back_init;

// cbct,pd1,back.c
// #define cbct_pd1_back1_args

typedef jool cbct_pd1_back1_type(
float *image,	// [nz nx ny] <- trick!
cint	nx,
cint	ny,
cint	nz,
crscalar dx,
crscalar dy,		// can be negative to cause flip
crscalar dz,
crscalar offset_x,	// center offset in pixels (usually 0)
crscalar offset_y,
crscalar offset_z_shift,	// offset_z - zshifts[ia]
cbyte	*mask2,		// [nx ny] 2D support mask: 0, 1, ..., nthread
cbyte	mask_id,	// 1 ... nthread
crscalar dso,	// distance from source to isocenter
crscalar dsd,	// distance from source to detector
crscalar dfs,	// distance from focal point to source (0 or inf)
cint	ns,
cint	nt,
crscalar ds,		// horizontal ray spacing
crscalar dt,		// vertical ray spacing
crscalar beta,		// source angle [radians]
crscalar offset_s,	// channel offset [pixels]
crscalar offset_t,	// vertical offset on detector [pixels]
cfloat	*proj);		// [nt ns] <- trick! projection view at angle beta

extern cbct_pd1_back1_type cbct_pd1_back1;
extern cbct_pd1_back1_type cbct_nn1_back1;


// cbct,pd1,back,t.c
extern jool cbct_pd1_back_t(
float *image, // [nz nx ny] <- trick!
const cbct_ig *ig,
const cbct_cg *cg,
cbct_work *cw,
cint na, // # of views
cfloat *proj, // [nt ns na] <- trick! projection views
cfloat *beta, // [na] source angles [radians]
cfloat *offset_s, // [na] possibly null
cfloat *offset_t, // [na] possibly null
cfloat *zshifts, // [na]
cint nthread, // # of threads
cint iblock, // for OS, which block
cint nblock, // for OS, #blocks
cjool compact, // 0 (for now)
cjool is_nn1,
cbct_back_init,
cfloat scale,
cint chat);


// cbct,pd1,proj.c
typedef jool cbct_pd1_proj1_type(
cfloat *image, // [nz, nx, ny] <- trick!
cint	nx,
cint	ny,
cint	nz,
crscalar dx,
crscalar dy,		// can be negative to cause flip
crscalar dz,
crscalar offset_x,	// center offset in pixels (usually 0)
crscalar offset_y,
crscalar offset_z_shift,	// offset_z - zshifts[ia]
cbyte	*mask2,		// [nx ny] 2D support mask: 0, 1, ..., nthread
crscalar dso,	// distance from source to isocenter
crscalar dsd,	// distance from source to detector
crscalar dfs,	// distance from focal point to source (0 or inf)
cint	ns,
cint	nt,
crscalar ds,		// horizontal ray spacing
crscalar dt,		// vertical ray spacing
crscalar beta,		// source angle [radians]
crscalar offset_s,	// channel offset [pixels]
crscalar offset_t,	// vertical offset on detector [pixels]
float	*proj);		// [nt ns] <- trick! projection view at angle beta

extern cbct_pd1_proj1_type cbct_pd1_proj1;
extern cbct_pd1_proj1_type cbct_nn1_proj1;

// cbct,pd1,proj,t.c
extern jool cbct_pd1_proj_t(
cfloat *image, // [nz nx ny] <- trick!
const cbct_ig *ig,
const cbct_cg *cg,
cbct_work *cw,
cint na, // # of views
float *proj, // [nt ns na] <- trick! projection views
cfloat *beta, // [na] source angles [radians]
cfloat *offset_s, // [na] possibly null
cfloat *offset_t, // [na] possibly null
cfloat *zshifts, // [na]
cint nthread, // # of threads
cint iblock, // for OS, which block
cint nblock, // for OS, #blocks
cjool compact, // 0 (for now)
cjool is_nn1,
cfloat scale,
cint chat);

// cbct,work.c
extern cbct_work *cbct_work_alloc(
const cbct_ig *ig,
const cbct_cg *cg,
cint nthread);

extern jool cbct_work_free(cbct_work *cw, cint nthread);

#endif // jf_cbct_def_h
