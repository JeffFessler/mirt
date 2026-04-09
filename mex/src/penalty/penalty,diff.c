/*
* penalty,diff.c
* compute differences between neighboring pixels for roughness penalties
* Copyright 2004-05-14, Jeff Fessler, The University of Michigan
*/
#include "defs-misc.h"


/*
* penalty_diff_dxyz()
* From the scalar offset value, determine the dx,dy,dz neighbor offsets.
* These are typically -1,0,+1 for the usual nearest neighbors.
* akin to matlab's ind2sub
*/
static sof penalty_diff_dxyz(
int	*p_dx,
int	*p_dy,
int	*p_dz,
cint	nx,
cint	ny,
cint	nz,
cint	offset)
{
	/* Upper bound on dx etc. is (almost) half the image size. */
	cint hx = Max(floor(nx/2), 1);
	cint hy = Max(floor(ny/2), 1);
	cint hz = Max(floor(nz/2), 1);

	/* trick: these next lines allow dx,dy,dz to be (somewhat) negative */
	cint dz = floor((offset + hx-1 + (hy-1)*ny + nz*nx*ny) / (nx*ny)) - nz;
	cint dy = floor((offset - dz*nx*ny + hx-1 + ny*nx) / nx) - ny;
	cint dx = offset - dz*nx*ny - dy*nx;

	if (offset < 0) Fail("need offset > 0")

	if (offset != dx + dy*nx + dz*nx*ny) Fail("offset bug")

	/* check offset bounds */
	if (	abs(dx) >= hx
	||	abs(dy) >= hy
	||	abs(dz) >= hz )
		Fail4("dx=%d dy=%d dz=%d for offset=%d", dx, dy, dz, offset)

	/* warn if the user is asking for "unusual" neighbors */
	if (	(dx != 1 && dx != -1 && dx != 0)
	||	(dy != 1 && dy != -1 && dy != 0)
	||	(dz != 1 && dz != -1 && dz != 0) )
		Warn4("dx=%d dy=%d dz=%d for offset=%d", dx, dy, dz, offset)

	*p_dx = dx;
	*p_dy = dy;
	*p_dz = dz;
	Ok
}


/*
* penalty_diff1_wk_tight_offset()
* out[ii] = kappa[ii] * kappa[ii-offset] / distance^distance_power
*/
static sof penalty_diff1_wk_tight_offset(
float	*po,	/* [[nx,ny,nz]] */
cfloat	*kappa,	/* [[nx,ny,nz]] */
cint	nx,
cint	ny,
cint	nz,
cint	offset,
cdouble	distance_power) /* if 1, yields "usual" 1/dis |x_j - x_k|^2,
		or use 2 for | (x_j - x_k) / dis |^2 which may be preferable */
{
	int ix, iy, iz;
	int dx, dy, dz;
	int mx, my, mz;
	float distance, distance_factor;

	Call(penalty_diff_dxyz, (&dx, &dy, &dz, nx, ny, nz, offset))

	mz = Min(nz, nz+dz);
	my = Min(ny, ny+dy);
	mx = Min(nx, nx+dx);
	distance = sqrt(dx*dx + dy*dy + dz*dz);
	distance_factor = 1. / pow(distance, distance_power);

	Bzero(po, nx*ny*nz)
	for (iz=Max(dz,0); iz < mz; ++iz)
		for (iy=Max(dy,0); iy < my; ++iy)
			for (ix=Max(dx,0); ix < mx; ++ix) {
				cint ii = ix + iy*nx + iz*nx*ny;
				po[ii] = kappa[ii] * kappa[ii-offset]
					* distance_factor;
			}
	Ok
}


/*
* penalty_diff2_wk_tight_offset()
* out[ii] = kappa[ii] * sqrt(kappa[ii-offset] * kappa[ii+offset])
* / distance^distance_power, for appropriate (non-boundary) ii
*/
static sof penalty_diff2_wk_tight_offset(
float	*po,	/* [[nx,ny,nz]] */
cfloat	*kappa,	/* [[nx,ny,nz]] */
cint	nx,
cint	ny,
cint	nz,
cint	offset,
cdouble	distance_power) /* if 1, yields "usual" 1/dis | . |^2,
		or use 2 for | . / dis |^2 which may be preferable */
{
	int ix, iy, iz;
	int dx, dy, dz;
	int mx, my, mz;
	float distance, distance_factor;

	Call(penalty_diff_dxyz, (&dx, &dy, &dz, nx, ny, nz, offset))
	dx = Abs(dx); /* make positive since both directions anyway */
	dy = Abs(dy);
	dz = Abs(dz);

	mx = nx-dx;
	my = ny-dy;
	mz = nz-dz;
	distance = sqrt(dx*dx + dy*dy + dz*dz);
	distance_factor = 1. / pow(distance, distance_power);

	Bzero(po, nx*ny*nz)
	for (iz=dz; iz < mz; ++iz)
		for (iy=dy; iy < my; ++iy)
			for (ix=dx; ix < mx; ++ix) {
				cint ii = ix + iy*nx + iz*nx*ny;
				po[ii] = distance_factor * kappa[ii]
				* sqrt(kappa[ii-offset] * kappa[ii+offset]);
			}
	Ok
}


/*
* penalty_diff_wk_tight()
*/
sof penalty_diff_wk_tight(
float	*po,	/* [[dim_i],noffset], each mask */
cfloat	*kappa,	/* [[dim_i]] "kappa" array or mask */
cint	*dim_i,	/* [nod_i] */
cint	nod_i,	/* # of dim of mask */
cint	*offset, /* [noffset] */
cint	noffset,
cdouble	distance_power,
cint	order)
{
	int ioff;
	cint nx = Max(dim_i[0], 1);
	cint ny = nod_i > 1 ? dim_i[1] : 1;
	cint nz = nod_i > 2 ? dim_i[2] : 1;
	cint nn = nx * ny * nz;
	if (nod_i > 3) Fail(">3D not done")
	for (ioff=0; ioff < noffset; ++ioff)
		if (order == 1)
			penalty_diff1_wk_tight_offset(po+ioff*nn, kappa,
			nx, ny, nz, offset[ioff], distance_power);
		else
			penalty_diff2_wk_tight_offset(po+ioff*nn, kappa,
			nx, ny, nz, offset[ioff], distance_power);
	Ok
}



/*
* penalty_diff1_wk_leak_offset()
* out[ii] = sqrt(kappa[ii] * kappa[ii-offset] * distance^distance_power)
*/
static sof penalty_diff1_wk_leak_offset(
float	*po,	/* [[nx,ny,nz]] */
cfloat	*kappa,	/* [[nx,ny,nz]] */
cint	nx,
cint	ny,
cint	nz,
cint	offset,
cdouble	distance_power) /* if 1, yields "usual" 1/dis |x_j - x_k|^2,
		or use 2 for | (x_j - x_k) / dis |^2 which may be preferable */
{
	int ix, iy, iz;
	int dx, dy, dz;
	int mx, my, mz;
	float distance, distance_factor;

	Call(penalty_diff_dxyz, (&dx, &dy, &dz, nx, ny, nz, offset))

	mz = Min(nz, nz+dz);
	my = Min(ny, ny+dy);
	mx = Min(nx, nx+dx);
	distance = sqrt(dx*dx + dy*dy + dz*dz);
	distance_factor = 1. / pow(distance, distance_power);

	Bzero(po, nx*ny*nz)
	for (iz=Max(dz,0); iz < mz; ++iz)
		for (iy=Max(dy,0); iy < my; ++iy)
			for (ix=Max(dx,0); ix < mx; ++ix) {
				cint ii = ix + iy*nx + iz*nx*ny;
				cdouble k0 = kappa[ii];
				cdouble k1 = kappa[ii-offset];
				double kk = 0.;
				if (k0 && k1)
					kk = k0 * k1;
				else if (k0)
					kk = k0 * k0;
				else if (k1)
					kk = k1 * k1;
				po[ii] = kk * distance_factor;
			}
	Ok
}


/*
* penalty_diff_wk_leak()
*/
sof penalty_diff_wk_leak(
float	*po,	/* [[dim_i],noffset], each mask */
cfloat	*kappa,	/* [[dim_i]] "kappa" array or mask */
cint	*dim_i,	/* [nod_i] */
cint	nod_i,	/* # of dim of mask */
cint	*offset, /* [noffset] */
cint	noffset,
cdouble	distance_power,
cint	order)
{
	int ioff;
	cint nx = Max(dim_i[0], 1);
	cint ny = nod_i > 1 ? dim_i[1] : 1;
	cint nz = nod_i > 2 ? dim_i[2] : 1;
	cint nn = nx * ny * nz;
	if (nod_i > 3) Fail(">3D not done")

	if (kappa[0])
		Warn("'leak' option is questionable for mask with non-zero outer border.")

	for (ioff=0; ioff < noffset; ++ioff)
		if (order == 1)
			penalty_diff1_wk_leak_offset(po+ioff*nn, kappa,
				nx, ny, nz, offset[ioff], distance_power);
		else
			Fail("not done")
	Ok
}


/*
* penalty_diff1_forw1_offset()
* out[ii] = image[ii] - image[ii-offset] for ii=offset to nn-1
*/
static void penalty_diff1_forw1_offset(
float	*po,	/* [[nn]] */
cfloat	*pi,	/* [[nn]] */
cint	offset,	/* caution: must be positive! */
cint	nn)
{
	register int ii;
	register cfloat *pd = pi;
	Bzero(po, offset)
	pi += offset;
	po += offset;
	for (ii=nn-offset; ii; --ii)
		*po++ = *pi++ - *pd++;
}


#if 0
/*
* penalty_diff1_forw1_anti()
* out[ii] = image[ii] - image[ii+offset] for ii=0 to nn-offset-1
* anti-causal version
*/
typedef enum {anti_copy_end, anti_zero_end} anti_end_type;
static void penalty_diff1_forw1_anti(
float	*po,	/* [[nn]] */
cfloat	*pi,	/* [[nn]] */
cint	offset,	/* caution: must be positive! */
cint	nn,
Const anti_end_type end_type)
{
	register int ii;
	register cfloat *pd = pi + offset;
	for (ii=nn-offset; ii; --ii)
		*po++ = *pi++ - *pd++;
	if (end_type == anti_copy_end)
		Bcopy(pi, po, offset)
	else
		Bzero(po, offset)
}
#endif


/*
* penalty_diff2_forw1_offset()
* 2nd finite differences.  for ii=offset to nn-offset-1:
* out[ii] = 2*image[ii] - image[ii+offset] - image[ii-offset]
*/
static void penalty_diff2_forw1_offset(
float	*po,	/* [[nn]] */
cfloat	*pi,	/* [[nn]] */
cint	offset,	/* caution: must be positive! */
cint	nn)
{
	register int ii;
	register cfloat *pl = pi;
	register cfloat *pr = pi + 2*offset;

	/* handle bad cases gracefully */
	if (offset <= 0 || 2*offset > nn) {
		Bzero(po, nn)
		return;
	}

	Bzero(po, offset)
	pi += offset;
	po += offset;
	for (ii=nn-2*offset; ii; --ii)
		*po++ = 2 * *pi++ - *pl++ - *pr++;
	Bzero(po, offset)
}


/*
* penalty_diff1_forw2_offset()
* c_{kj}^2 = 1, so
* out[ii] = image[ii] + image[ii-offset] for ii=offset to nn-offset-1
*/
static void penalty_diff1_forw2_offset(
float	*po,	/* [[nn]] */
cfloat	*pi,	/* [[nn]] */
cint	offset,
cint	nn)
{
	register int ii;
	register cfloat *pd = pi;
	pi += offset;
	po += offset;
	for (ii=nn-offset; ii; --ii)
		*po++ = *pi++ + *pd++;
}


/*
* penalty_diff1_back1_offset()
* adjoint of:
* image[ii] - image[ii-offset] for ii=offset to nn-1
*/
static void penalty_diff1_back1_offset(
float	*po,	/* [[nn]] */
cfloat	*pi,	/* [[nn]] */
cint	offset,
cint	nn)
{
	register int ii;
	register float *pd = po;
	pi += offset;
	po += offset;
	for (ii=nn-offset; ii; --ii) {
		*po++ += *pi;
		*pd++ -= *pi++;
	}
}


/*
* penalty_diff2_back1_offset()
* 2nd finite differences.  for ii=offset to nn-offset-1:
* adjoint of: 2*image[ii] - image[ii+offset] - image[ii-offset]
*/
static void penalty_diff2_back1_offset(
float	*po,	/* [[nn]] */
cfloat	*pi,	/* [[nn]] */
cint	offset,	/* caution: must be positive! */
cint	nn)
{
	register int ii;
	register float *pl = po;
	register float *pr = po + 2*offset;

	if (offset <= 0 || 2*offset > nn)
		return; /* handle bad cases gracefully */

	pi += offset;
	po += offset;
	for (ii=nn-2*offset; ii; --ii) {
		register cfloat value = *pi++;
		*pl++ -= value;
		*pr++ -= value;
		*po++ += 2 * value;
	}
}


/*
* penalty_diff1_back2_offset()
* c_{kj}^2 = 1, so
* \sum_k wk c_{kj}^2
*/
static void penalty_diff1_back2_offset(
float	*po,	/* [[nn]] */
cfloat	*pi,	/* [[nn]] wk */
cint	offset,
cint	nn)
{
	register int ii;
	register float *pd = po;
	pi += offset;
	po += offset;
	for (ii=nn-offset; ii; --ii) {
		*po++ += *pi;
		*pd++ += *pi++;
	}
}


/*
* penalty_diff_forw()
*/
sof penalty_diff_forw(
float	*po,	/* [[nn],noffset,nrep], each set of differences */
cfloat	*pi,	/* [[nn],nrep] */
cint	nn,
cint	*offset, /* [noffset] */
cint	noffset,
cint	nrep,
cint	order,	/* 1st or 2nd order differences */
cint	power)	/* 1 for ckj or 2 for ckj^2 */
{
	int ioff, irep;
	void (*fun)(float *, cfloat *, cint, cint);

	if (order == 1 && power == 1)
		fun = penalty_diff1_forw1_offset;
	else if (order == 2 && power == 1)
		fun = penalty_diff2_forw1_offset;
	else if (order == 1 && power == 2)
		fun = penalty_diff1_forw2_offset;
	else
		Fail2("unknown order=%d power=%d combination", order, power)

	for (irep=0; irep < nrep; ++irep, pi += nn)
		for (ioff=0; ioff < noffset; ++ioff, po += nn) {
			if (offset[ioff] > 0)
				fun(po, pi, offset[ioff], nn);
			else if (offset[ioff] < 0)
				Fail1("bad offset %d", offset[ioff])
			else /* offset = 0 */
				Bcopy(pi, po, nn)
		}
	Ok
}


/*
* penalty_diff_back()
*/
sof penalty_diff_back(
float	*po,	/* [[nn]] */
cfloat	*pi,	/* [[nn],noffset], each set of differences  */
cint	nn,
cint	*offset, /* [noffset] */
cint	noffset,
cint	nrep,
cint	order,
cint	power)
{
	int ioff, irep;
	void (*fun)(float *, cfloat *, cint, cint);

	if (order == 1 && power == 1)
		fun = penalty_diff1_back1_offset;
	else if (order == 2 && power == 1)
		fun = penalty_diff2_back1_offset;
	else if (order == 1 && power == 2)
		fun = penalty_diff1_back2_offset;
	else
		Fail2("unknown order=%d power=%d combination", order, power)

	Bzero(po, nn * nrep)
	for (irep=0; irep < nrep; ++irep, po += nn)
		for (ioff=0; ioff < noffset; ++ioff, pi += nn) {
			if (offset[ioff] > 0)
				fun(po, pi, offset[ioff], nn);
			else if (offset[ioff] < 0)
				Fail1("bad offset %d", offset[ioff])
			else /* offset = 0 */
				VectAdd(float, po, pi, nn)
		}
	Ok
}
