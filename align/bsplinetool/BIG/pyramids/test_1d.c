/* ----------------------------------------------------------------------------
	Filename:  	test_1d.c
	
	Author:		Daniel Sage
				Swiss Federal Institute of Technology - Lausanne
				Biomedical Imaging Group
				EPFL/DMT/IOA, BM-Ecublens, CH-1015 Lausanne, Switzerland
				
	Date:		20 November 1999
	
	Purpose:	Test program for the pyramid routines.
				The program build a 1D signal, reduce the signal, expand the 
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
#define LENGTH 10L		/* Length of the input and output 1D signal, should be even */
							
#define SPLINE			"Spline"			/* Spline filter (l2-norm) */
#define SPLINE_L2		"Spline L2"			/* Spline filter (L2-norm) */
#define SPLINE_CENT		"Centered Spline"	/* Centered Spline filter (l2-norm) */
#define SPLINE_CENT_L2	"Centered Spline L2"/* Centered Spline filter (L2-norm) */

/* --- Main --- */
int main(void)
{
double	InputSignal[LENGTH];		/* Input signal */
double 	ReducedSignal[LENGTH/2L];	/* Reduced signal */
double 	ExpandedSignal[LENGTH];		/* Expanded signal */
double  g[MAXF];					/* Coefficients of the reduce filter */
long  	ng;							/* Number of coefficients of the reduce filter */
double	h[MAXF];					/* Coefficients of the expansion filter */
long 	nh;							/* Number of coefficients of the expansion filter */
long 	i;							/* Index loop */
short	IsCentered;					/* Equal TRUE if the filter is a centered spline, FALSE otherwise */

	/* Get the filter coefficients for the Spline (order = 3) filter*/
	printf("Filter: " SPLINE "\n");
	if (GetPyramidFilter( SPLINE_CENT, 3L, g, &ng, h, &nh, &IsCentered) == ERROR) {
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
	InputSignal[4] =  2.0;
	InputSignal[5] =  1.0;
	InputSignal[6] =  0.0;
	InputSignal[7] =  -2.0;
	InputSignal[8] =  -4.0;
	InputSignal[9] =  -6.0;

	printf("Input   : "); 
	for (i=0L; i<LENGTH; i++) 
		printf("%2.3f ", InputSignal[i]);
	printf("\n");

	/* Reducing */
	Reduce_1D(InputSignal, LENGTH, ReducedSignal, g, ng, IsCentered);

	printf("Reduced : "); 
	for (i=0L; i<LENGTH/2L; i++) 
		printf("%2.3f ", ReducedSignal[i]);
	printf("\n");
	
	/* Expanding */
	Expand_1D(ReducedSignal, LENGTH/2L, ExpandedSignal, h, nh, IsCentered);

	printf("Expanded: "); 
	for (i=0L; i<LENGTH; i++) 
		printf("%2.3f ", ExpandedSignal[i]);
	printf("\n");

	/* Computing the error */
	printf("Error   : "); 
	for (i=0L; i<LENGTH; i++) 
		printf("%2.3f ", ExpandedSignal[i] - InputSignal[i]);
	printf("\n");

	return 0;
}

