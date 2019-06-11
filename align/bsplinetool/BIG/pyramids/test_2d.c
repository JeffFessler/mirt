/* ----------------------------------------------------------------------------
	Filename:  	test_2d.c
	
	Author:		Daniel Sage
				Swiss Federal Institute of Technology - Lausanne
				Biomedical Imaging Group
				EPFL/DMT/IOA, BM-Ecublens, CH-1015 Lausanne, Switzerland
				
	Date:		20 November 1999
	
	Purpose:	Test program for the pyramid routines,
				The program build a 2D signal, reduce the signal, expand the 
				signal then compute the error with the original signal.

	Convention:	g[ng] notes the reduce filter
				h[nh] notes the expansion filter
						 			 
---------------------------------------------------------------------------- */

/* --- System includes --- */
#include <stdio.h>
#include <stdlib.h>

/* --- Private includes --- */
#include "configs.h"
#include "pyramidtools.h"

/* --- Defines --- */
#define MAXF 200L		/* Maximum size of the filter */
#define NX 4L			/* Size of the image (X Axis), should be even */
#define NY 4L			/* Size of the Image (Y Axis), should be even */
							
#define SPLINE			"Spline"			/* Spline filter (l2-norm) */
#define SPLINE_L2		"Spline L2"			/* Spline filter (L2-norm) */
#define SPLINE_CENT		"Centered Spline"	/* Centered Spline filter (l2-norm) */
#define SPLINE_CENT_L2	"Centered Spline L2"/* Centered Spline filter (L2-norm) */

/* --- Main --- */
int main(void)
{
float	InputSignal[NX*NY];			/* Input signal */
float 	ReducedSignal[NX/2L*NY/2L];	/* Reduced signal */
float 	ExpandedSignal[NX*NY];		/* Expanded signal */
double  g[MAXF];					/* Coefficients of the reduce filter */
long  	ng;							/* Number of coefficients of the reduce filter */
double	h[MAXF];					/* Coefficients of the expansion filter */
long 	nh;							/* Number of coefficients of the expansion filter */
long 	i, j;						/* Index loop */
short	IsCentered;					/* Equal TRUE if the filter is a centered spline, FALSE otherwise */

	/* Get the filter coefficients for the Spline (order = 3) filter*/
	printf("Filter: " SPLINE "\n");
	if (GetPyramidFilter( SPLINE, 3, g, &ng, h, &nh, &IsCentered) == ERROR) {
		printf("Unable to load the filter coeffiients");
		return 1;
	}
	printf("Size of the reduce filter: %ld\n", ng);
	printf("Size of the expand filter: %ld\n", nh);

	/* Creation of a 1D input */
	InputSignal[0] =  0.0;
	InputSignal[1] =  1.0;
	InputSignal[2] =  2.0;
	InputSignal[3] =  3.0;
	InputSignal[4] =  1.0;
	InputSignal[5] =  2.0;
	InputSignal[6] =  3.0;
	InputSignal[7] =  4.0;
	InputSignal[8] =  2.0;
	InputSignal[9] =  3.0;
	InputSignal[10] =  4.0;
	InputSignal[11] =  5.0;
	InputSignal[12] =  3.0;
	InputSignal[13] =  4.0;
	InputSignal[14] =  3.0;
	InputSignal[15] =  2.0;
	
	printf("Input : \n"); 
	for (j=0L; j<NY; j++) { 
		for (i=0L; i<NX; i++) 
			printf("%2.3f ", InputSignal[i+j*NX]);
		printf("\n");
	}
	
	/* Reducing */
	Reduce_2D(InputSignal, NX, NY, ReducedSignal, g, ng, IsCentered);

	printf("Reduced: \n"); 
	for (j=0L; j<NY/2L; j++) { 
		for (i=0L; i<NX/2L; i++) 
			printf("%2.3f ", ReducedSignal[i+j*NX/2L]);
		printf("\n");
	}
	
	/* Expanding */
	Expand_2D(ReducedSignal, NX/2L, NY/2L, ExpandedSignal, h, nh, IsCentered);

	printf("Expanded : \n"); 
	for (j=0L; j<NY; j++) { 
		for (i=0L; i<NX; i++) 
			printf("%2.3f ", ExpandedSignal[i+j*NX]);
		printf("\n");
	}

	/* Computing the error */
	Expand_2D(ReducedSignal, NX/2L, NY/2L, ExpandedSignal, h, nh, IsCentered);

	printf("Error : \n"); 
	for (j=0L; j<NY; j++) { 
		for (i=0L; i<NX; i++) 
			printf("%2.3f ", ExpandedSignal[i+j*NX] - InputSignal[i+j*NX]);
		printf("\n");
	}
	
	
	return 0;
}

