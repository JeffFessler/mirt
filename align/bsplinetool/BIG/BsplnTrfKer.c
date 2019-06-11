/*******************************************************************************
        Trivial extension of BsplnTrf.h in registration package at
        http://bigwww.epfl.ch/algorithms.html so that any kernel
        can be used in inverse filtering (or interpolation),
        with mirror end condition

        Coded by Se Young Chun, Jan 19, 2007, The Univ. of Michigan

*******************************************************************************/

#include		<stddef.h>
#include		<stdlib.h>
#include		<string.h>
#include		<math.h>
#include		"BsplnTrfKer.h"

/* you need to include registration package from 
   http://bigwww.epfl.ch/thevenaz/registration/    */
#include		"reg/phil.h"
#include		"reg/convolve.h"
#include		"reg/getPut.h"

int	inverseBsplineKernelFinite (
		double	*inPtr,
		double	*outPtr,
		long	nx,
		long	ny,
		long	nz,
		double 	*kernel,
		long	nk) 
{
	double	*inData, *outData;
	long	x, y, z;

	outData = outPtr;
	for (z = 0L; (z < nz); z++)
		for (y = 0L; (y < ny); inPtr += (ptrdiff_t)nx, 
			outData += (ptrdiff_t)nx, y++)
			firConvolveFinite(inPtr, outData, nx, 0L, kernel, nk);

	if (ny > 1L) {
	 	inData = (double *)malloc((size_t)ny * sizeof(double));
	  	if (inData == (double *)NULL) {
			message("ERROR - Unable to allocate y inData");
			return(ERROR);
		}
		outData = (double *)malloc((size_t)ny * sizeof(double));
		if (outData == (double *)NULL) {
			free(inData);
			message("ERROR - Unable to allocate y outData");
			return(ERROR);
		}
		for (z = 0L; (z < nz); z++)
			for (x = 0L; (x < nx); x++) {
			  	getyD2D(outPtr, x, 0L, z, nx, ny, inData, ny);
			  	firConvolveFinite(inData, outData, ny, 0L, 
					kernel, nk);
			  	putyD2D(outPtr, x, 0L, z, nx, ny, outData, ny);
			}

		free(outData);
		free(inData);
	}

	if (nz > 1L) {
		inData = (double *)malloc((size_t)nz * sizeof(double));
		if (inData == (double *)NULL) {
			message("ERROR - Unable to allocate z inData");
			return(ERROR);
		}
		outData = (double *)malloc((size_t)nz * sizeof(double));
		if (outData == (double *)NULL) {
			free(inData);
			message("ERROR - Unable to allocate z outData");
			return(ERROR);
		}
		for (y = 0L; (y < ny); y++)
			for (x = 0L; (x < nx); x++) {
			  	getzD2D(outPtr, x, y, 0L, nx, ny, inData, nz);
			  	firConvolveFinite(inData, outData, nz, 0L, 
					kernel, nk);
			  	putzD2D(outPtr, x, y, 0L, nx, ny, outData, nz);
			}

		free(outData);
		free(inData);
	}

	return(!ERROR);

} /* End of inverseBsplineKernelFinite */



int	inverseBsplineKernel3Finite (
		double	*inPtr,
		double	*outPtr,
		long	nx,
		long	ny,
		long	nz,
		double 	*kernelx,
		long	nkx, 
		double 	*kernely,
		long	nky, 
		double 	*kernelz,
		long	nkz) 
{
	double	*inData, *outData;
	long	x, y, z;

	outData = outPtr;
	for (z = 0L; (z < nz); z++)
		for (y = 0L; (y < ny); inPtr += (ptrdiff_t)nx, 
			outData += (ptrdiff_t)nx, y++)
			firConvolveFinite(inPtr, outData, nx, 0L, kernelx, nkx);

	if (ny > 1L) {
		inData = (double *)malloc((size_t)ny * sizeof(double));
		if (inData == (double *)NULL) {
			message("ERROR - Unable to allocate y inData");
			return(ERROR);
		}
		outData = (double *)malloc((size_t)ny * sizeof(double));
		if (outData == (double *)NULL) {
			free(inData);
			message("ERROR - Unable to allocate y outData");
			return(ERROR);
		}
		for (z = 0L; (z < nz); z++)
			for (x = 0L; (x < nx); x++) {
			  	getyD2D(outPtr, x, 0L, z, nx, ny, inData, ny);
			  	firConvolveFinite(inData, outData, ny, 0L, 
					kernely, nky);
			  	putyD2D(outPtr, x, 0L, z, nx, ny, outData, ny);
			}

		free(outData);
		free(inData);
	}
	if (nz > 1L) {
		inData = (double *)malloc((size_t)nz * sizeof(double));
		if (inData == (double *)NULL) {
			message("ERROR - Unable to allocate z inData");
			return(ERROR);
		}
		outData = (double *)malloc((size_t)nz * sizeof(double));
		if (outData == (double *)NULL) {
			free(inData);
			message("ERROR - Unable to allocate z outData");
			return(ERROR);
		}
		for (y = 0L; (y < ny); y++)
			for (x = 0L; (x < nx); x++) {
			  	getzD2D(outPtr, x, y, 0L, nx, ny, inData, nz);
			  	firConvolveFinite(inData, outData, nz, 0L, 
					kernelz, nkz);
			  	putzD2D(outPtr, x, y, 0L, nx, ny, outData, nz);
			}

		free(outData);
		free(inData);
	}

	return(!ERROR);

} /* End of inverseBsplineKernel3Finite */
