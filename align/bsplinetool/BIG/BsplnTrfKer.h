/*******************************************************************************

	Trivial extension of BsplnTrf.h in registration package at
	http://bigwww.epfl.ch/algorithms.html so that any kernel 
	can be used in inverse filtering (or interpolation), 
	with zero end condition

        Coded by Se Young Chun, Jan 19, 2007, The Univ. of Michigan

*******************************************************************************/

extern	int	inverseBsplineKernelFinite(
			double		*inPtr,
			double		*outPtr,
			long		nx,
			long		ny,
			long		nz,
			double		*kernel,
			long		nk);

extern	int	inverseBsplineKernel3Finite(
			double		*inPtr,
			double		*outPtr,
			long		nx,
			long		ny,
			long		nz,
			double		*kernelx,
			long		nkx,
			double		*kernely,
			long		nky,
			double		*kernelz,
			long		nkz);
