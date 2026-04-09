// cbct,mask2.c
// mask operations

#include "cbct,def.h"


// cbct_mask_iy_init()
// ithread responsible for iy_start[ithread] <= iy < iy_end[ithread]
sof cbct_mask_iy_init(
int *iy_start, // [nthread]
int *iy_end, // [nthread]
cbyte *mask2, // [nx ny]
cint nx,
cint ny,
cint nthread,
cint chat)
{
	for (int ithread=0; ithread < nthread; ++ithread)
	{
		iy_start[ithread] = -1;
		iy_end[ithread] = -1;
	}

	for (int iy=0; iy < ny; ++iy)
	{
		for (int ix=0; ix < nx; ++ix, ++mask2)
		{
			cint ithread = *mask2 - 1;
			if (ithread == -1)
				continue;
			if (iy_start[ithread] == -1)
				iy_start[ithread] = iy;
			iy_end[ithread] = iy+1;
		}
	}

	if (chat)
	{
		for (int ithread=0; ithread < nthread; ++ithread)
			Note3("ithread=%d iy range [%d,%d)", ithread,
				iy_start[ithread], iy_end[ithread])
	}

	Ok
}


// cbct_mask_init()
// initialize mask with values for each thread in it
// input and output masks can be the same pointer!
sof cbct_mask_init(
byte *mask_int,	// [nx ny] out: integers 0 or 1 ... nthread
cbyte *mask_bin, // [nx ny] in: binary, possibly null
cint nx, cint ny,
cint nthread,
ctruf chat)
{
	if (nthread > 250) Fail("too many threads")

	long sum = 0; // # of x-y pixels in 2d mask
	if (!mask_bin)
		sum = nx * ny;
	else
		for (int jj=0; jj < nx * ny; ++jj)
			sum += mask_bin[jj] != 0;

	cint nj_per_thread = (int) Ceilf(sum / (float) nthread);
	if (chat > 99)
		Note3("%ld xy pixels / %d threads = %d each",
			sum, nthread, nj_per_thread)

	int jj=0;
	for (int it=1; it <= nthread; ++it)
	{
		int nj_here = Min(sum, it * nj_per_thread) - (it-1) * nj_per_thread;
		nj_here = Max(nj_here, 0); // fixed 2011-10-12

		if (chat > 99) Note2("thread %d nj %d", it-1, nj_here)

		while (nj_here)
		{
			if (!mask_bin || mask_bin[jj])
			{
				mask_int[jj] = it;
				--nj_here;
			}
			else
				mask_int[jj] = 0;
			++jj;
		}
	}
	// if (jj != nx * ny) Fail("bug")
	// Iwrite2byte("mask-int.fld", mask_int, nx, ny)
	Ok
}
