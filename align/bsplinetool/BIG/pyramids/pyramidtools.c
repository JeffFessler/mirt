/* ----------------------------------------------------------------------------
	Filename:  	pyramidtools.h
	
	Author:		Daniel Sage
				Swiss Federal Institute of Technology - Lausanne
				Biomedical Imaging Group
				EPFL/DMT/IOA, BM-Ecublens, CH-1015 Lausanne, Switzerland
				
	Date:		17 March 1999
	
	Purpose:	Basic functions reduce and expand by factors of 2 in 1D 
				signal or in 2D signal.
				Includes the standard pyramid and the centered pyramid.
				
	References:
	
	[1]	M. Unser, "Splines: a perfect fit for signal and image processing," 
		IEEE Signal Processing Magazine, 1999.
			
	[2]	M. Unser, A. Aldroubi and M. Eden, "B-spline signal processing: 
		Part II‹efficient design and applications," IEEE Trans. Signal Processing, 
		vol. 41, no. 2, pp. 834-848, February 1993.
			
	[3]	M. Unser, A. Aldroubi and M. Eden, "The L2 polynomial spline pyramid," 
		IEEE Trans. Pattern Anal. Mach. Intell., vol. 15, no. 4, pp. 364-379, April 1993.
			
	[4]	P. Brigger, F. Müller, K. Illgner and M. Unser, "Centered pyramids," 
		IEEE Trans. Image Processing, vol. 6, no. 9, pp. 1254-1264, September 1999.
			
	[5]	P.J. Burt and E.H. Adelson, "The Laplacian pyramid as a compact code," 
		IEEE Trans. Commun., vol. COM-31, no. 4, pp. 337-345, April 1983.

			 
---------------------------------------------------------------------------- */

/* --- System includes --- */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

/* --- Private includes --- */
#include "configs.h"
#include "messagedisplay.h"
#include "pyramidtools.h"
#include "pyramidfilters.h"

/* ------------------------------------------------------------------------- */
/* Declaration of static procedures                                          */
/* ------------------------------------------------------------------------- */
static void ExpandStandard_1D(	
				double In[], long NxIn,
				double Out[],
				double w[], long nw
				);
						
static void ReduceStandard_1D(	
				double In[], long NxIn,
				double Out[],
				double w[], long nw
				);
						
static void ExpandCentered_1D(	
				double In[], long NxIn,
 				double Out[],
 				double w[], long nw
 				);

static void ReduceCentered_1D(	
				double In[], long NxIn,
 				double Out[],
 				double w[], long nw
 				);

static void GetRow( 
				float *Image, long Nx, long Ny,
				long RowNb,
				double *Row, long RowSize
				);
				
static void GetColumn( 
				float *Image, long Nx, long Ny,
				long ColumnNb,
				double *Column, long ColumnSize
				);
				
static void PutRow( 
				float *Image, long Nx, long Ny,
				long RowNb,
				double *Row, long RowSize
				);
				
static void PutColumn( 
				float *Image, long Nx, long Ny,
				long ColumnNb, 
				double *Column, long ColumnSize
				);
			
/* ----------------------------------------------------------------------------
	
	Function: 
		GetPyramidFilter
		
	Purpose:
		Get the coefficients of the filter (reduce and expand filter)
		Return the coefficients in g[ng] and in h[nh]
		
	Convention:
		g[ng] for the reduce filter
		h[nh] for the expansion filter
		
	Parameters:
		Filter is the name of the filter
		
		Order is the order for the filters based on splines
			For the "Spline" filter, Order is 0, 1, 2 or 3
			For the "Spline L2" filter, Order is 0, 1, 3 or 5
			For the "Centered Spline" filter, Order is 0, 1, 2, 3 or 4
			For the "Centered Spline L2" filter, Order is 0, 1, 2, 3 or 4
			
		IsCentered is a return value indicates if the filter is a centered filter
			TRUE if it is a centered filter
			FALSE if it is not a centered filter 
			
---------------------------------------------------------------------------- */
extern int GetPyramidFilter(
					char *Filter, 				
					long Order, 				
					double g[], long *ng,
					double h[], long *nh,
					short *IsCentered)		
{

	ng[0] = -1L;
	nh[0] = -1L;
	*IsCentered = FALSE;
		
	if ( !strcmp(Filter, "Spline"))	{	
		PyramidFilterSplinel2(g, ng, h, nh, Order); 
		*IsCentered = FALSE;	
	}
		
	if ( !strcmp(Filter, "Spline L2")) {
		PyramidFilterSplineL2(g, ng, h, nh, Order);
		*IsCentered = FALSE;	
	}
		
	if ( !strcmp(Filter, "Centered Spline")) {	
		PyramidFilterCentered(g, ng, h, nh, Order); 
		*IsCentered = TRUE;	
	}
		
 	if ( !strcmp(Filter, "Centered Spline L2"))	{
		PyramidFilterCenteredL2(g, ng, h, nh, Order); 
		*IsCentered = TRUE;	
	}
       
	if ( ng[0] == -1L && nh[0] == -1L) {
        MessageDisplay( "This familly filters is unknown");
        return(ERROR);
    }
    return( !ERROR);
    
}

/* ----------------------------------------------------------------------------
	
	Function: 
		Reduce_2D
	
	Purpose: 
 		Reduces an image by a factor of two in each dimension.
 		
	Note: 
 		Expects the output array (Out) to be allocated.

	Parameters:
		Input image:  	In[NxIn*NyIn]
		Output image: 	Out[NxIn/2*NyIn/2]
		Filter:			g[ng] coefficients of the filter
		
---------------------------------------------------------------------------- */
extern int Reduce_2D(	
				float *In, long NxIn, long NyIn,
				float *Out,
				double g[], long ng,
				short IsCentered
				)
{
float	*Tmp;
double 	*InBuffer;		/* Input buffer to 1D process */ 
double	*OutBuffer;		/* Output buffer to 1D process */ 
long	kx, ky;
long 	NxOut;
long 	NyOut;

	/* --- Define dimension of the output --- */
	NxOut = NxIn/2L;
	if (NxOut < 1L) NxOut = 1L;
	
	NyOut = NyIn/2L;
	if (NyOut < 1L) NyOut = 1L;

	/* --- Allocate a temporary image --- */
	Tmp = (float *)malloc((size_t)(NxOut*NyIn*(long)sizeof(float)));
	if (Tmp == (float *)NULL) {
		MessageDisplay("Unable to allocate memory");
		return(ERROR);
	}
	
	/* --- X processing --- */
	if (NxIn > 1L) {
		InBuffer = (double *)malloc((size_t)(NxIn*(long)sizeof(double)));
		if (InBuffer == (double *)NULL) {
			free(Tmp);
			MessageDisplay("Unable to allocate memory");
			return(ERROR);
		}
		OutBuffer = (double *)malloc((size_t)(NxOut*(long)sizeof(double)));
		if (OutBuffer == (double *)NULL) {
			free(Tmp);
			free(InBuffer);
			MessageDisplay("Unable to allocate memory");
			return(ERROR);
		}
		for (ky=0L; ky<NyIn; ky++) {
			GetRow( In, NxIn, NyIn, ky, InBuffer, NxIn); 
			Reduce_1D(InBuffer, NxIn, OutBuffer, g, ng, IsCentered); 
	    	PutRow( Tmp, NxOut, NyIn, ky, OutBuffer, NxOut);
    	}
		free(InBuffer);
		free(OutBuffer);
    }
    else
	   	memcpy( Tmp, In, (size_t)(NyIn*(long)sizeof(float)));
    	
	/* --- Y processing --- */
	if (NyIn > 1L) {
		InBuffer = (double *)malloc((size_t)(NyIn*(long)sizeof(double)));
		if (InBuffer == (double *)NULL) {
			MessageDisplay("Unable to allocate memory");
			return(ERROR);
		}
		OutBuffer = (double *)malloc((size_t)(NyOut*(long)sizeof(double)));
		if (OutBuffer == (double *)NULL) {
			free(InBuffer);
			MessageDisplay("Unable to allocate memory");
			return(ERROR);
		}
	    for (kx=0L; kx<NxOut; kx++) {
			GetColumn( Tmp, NxOut, NyIn, kx, InBuffer, NyIn); 
			Reduce_1D(InBuffer, NyIn, OutBuffer, g, ng,IsCentered);
	 	   	PutColumn( Out, NxOut, NyOut, kx, OutBuffer, NyOut); 
   		}
		free(InBuffer);
		free(OutBuffer);
	}
	else
	   	memcpy( Out, Tmp, (size_t)(NxOut*(long)sizeof(float)));
	
	/* --- Free the temporary image --- */
	free(Tmp);
	
	return (!ERROR);
}

/* ----------------------------------------------------------------------------
	
	Function: 
		Expand_2D
	
	Purpose: 
 		Expands an image by a factor of two in each dimension.
 		
	Note: 
 		Expects the output array (Out) to be allocated.

	Parameters:
		Input volume:  	In[NxIn,NyIn]
		Output voulme: 	Out[NxIn*2,NyIn*2]
		Filter coef:	h[nh]
		
---------------------------------------------------------------------------- */
extern int Expand_2D(	
				float *In, long NxIn, long NyIn,
				float *Out,
				double h[], long nh,
				short IsCentered
				)
{
double  *InBuffer; 		/* Input buffer to 1D process */ 
double  *OutBuffer;		/* Output buffer to 1D process */ 
long 	kx, ky;
long 	NxOut;
long 	NyOut;
	
	if (NxIn <= 1L) 
		NxOut = 1L; 
	else 
		NxOut = NxIn*2L;
		
	if (NyIn <= 1L) 
		NyOut = 1L; 
	else 
		NyOut = NyIn*2L;
 
	/* --- X processing --- */
	if (NxIn > 1L) {
		InBuffer = (double *)malloc((size_t)(NxIn*(long)sizeof(double)));
		if (InBuffer == (double *)NULL) {
			MessageDisplay("Unable to allocate memory");
			return(ERROR);
		}
		OutBuffer = (double *)malloc((size_t)(NxOut*(long)sizeof(double)));
		if (OutBuffer == (double *)NULL) {
			free(InBuffer);
			MessageDisplay("Unable to allocate memory");
			return(ERROR);
		}
		for (ky=0L; ky<NyIn; ky++) {	
			GetRow( In, NxIn, NyIn, ky, InBuffer, NxIn);
			Expand_1D(InBuffer, NxIn, OutBuffer , h, nh, IsCentered);
			PutRow( Out, NxOut, NyOut, ky, OutBuffer, NxOut);  
		}
		free(InBuffer);
		free(OutBuffer);
	}
	else
	   	memcpy( Out, In, (size_t)(NxIn*NyIn*(long)sizeof(float)));


	/* --- Y processing --- */
	if (NyIn > 1L) {
		InBuffer = (double *)malloc((size_t)(NyIn*(long)sizeof(double)));
		if (InBuffer == (double *)NULL) {
			MessageDisplay("Unable to allocate memory");
			return(ERROR);
		}
		OutBuffer = (double *)malloc((size_t)(NyOut*(long)sizeof(double)));
		if (OutBuffer == (double *)NULL) {
			free(InBuffer);
			MessageDisplay("Unable to allocate memory");
			return(ERROR);
		} 
		for (kx=0L; kx<NxOut; kx++) {
			GetColumn( Out, NxOut, NyOut, kx, InBuffer, NyIn); 
			Expand_1D(InBuffer, NyIn, OutBuffer, h, nh, IsCentered); 
    		PutColumn( Out, NxOut, NyOut, kx, OutBuffer, NyOut); 
		}
		free(InBuffer);
		free(OutBuffer);
	}
	 
	return (ERROR);
	
}

/* ----------------------------------------------------------------------------

	Function:	
		Reduce_1D
	
	Purpose:	
		Router function to call ReduceStandard_1D or ReduceCentered_1D
	
---------------------------------------------------------------------------- */
extern void Reduce_1D(	
				double In[], long NxIn,
				double Out[], 
				double g[], long ng,
				short IsCentered
				)
{

	if (IsCentered)
		ReduceCentered_1D( In, NxIn, Out, g, ng);
	else
		ReduceStandard_1D( In, NxIn, Out, g, ng);
}
/* ----------------------------------------------------------------------------

	Function:	
		Expand_1D
	
	Purpose:	
		Router function to call ExpandStandard_1D or ExpandCentered_1D
	
---------------------------------------------------------------------------- */
extern void Expand_1D(	
				double In[], long NxIn,
				double Out[],
				double h[], long nh,
				short IsCentered
				)
{

	if (IsCentered)
		ExpandCentered_1D( In, NxIn, Out, h, nh);
	else
		ExpandStandard_1D( In, NxIn, Out, h, nh);
}


/* ----------------------------------------------------------------------------

	Function:	
		ReduceStandard_1D
	
	Purpose:	
		Basic function to reduce a 1D signal
	
	Parameters:
		In[NxIn] is the input signal	(NxIn should be greater than 2 and even)
		Out[NxIn/2] is the output signal
		g[ng] is an array that contains the coefficients of the filter
		
	Author:  
		Michael Unser, NIH, BEIP, June 1994
		Daniel Sage, EPFL, Biomedical Imaging Group, April 1999
				
---------------------------------------------------------------------------- */
static void ReduceStandard_1D(	
				double In[], long NxIn,
				double Out[],
				double g[], long ng)
{
long k, i, i1, i2;
long kk, kn, nred, n;
	
	nred = NxIn/2L;
	n  = nred*2L;
	kn = n-1L;			/* kn=n-2; DS Modified */ 
	     
	if (ng<2L) {     	/* length filter < 2 */
		for (kk=0L; kk<nred; kk++) {
			k  = 2L*kk;
			i2 = k+1L;
			if (i2 > n-1L) 
				i2 = kn-i2;
			Out[kk] = (In[k]+In[i2])/2.;
		}
	}
	     
	else {
		for (kk=0L; kk<nred; kk++) {
			k = 2L*kk;
			Out[kk] = In[k]*g[0];
			for (i=1L; i<ng; i++) {
				i1 = k-i;
				i2 = k+i;
				if (i1<0L) {
					i1 = (-i1) % kn;
					if (i1 > n-1) 
						i1=kn-i1;
				}
				if (i2 > n-1L) {
					i2 = i2 % kn;
					if (i2 > n-1L) 
						i2=kn-i2;
				}
				Out[kk] = Out[kk] + g[i]*(In[i1]+In[i2]);  
			}
		}
	}
	
}
	
/* ----------------------------------------------------------------------------

	Function:	
		ExpandStandard_1D
	
	Purpose:	
		Basic function to expand a 1D signal
	
	Parameters:
		In[NxIn] is the input signal	(NxIn should be greater than 1)
		Out[NxIn*2] is the output signal
		w[nw] is an array that contains the coefficients of the filter
		
	Author:  
		Michael Unser, NIH, BEIP, June 1994
		Daniel Sage, EPFL, Biomedical Imaging Group, April 1999
				
---------------------------------------------------------------------------- */
static void ExpandStandard_1D(	
				double In[], long NxIn,
				double Out[], 
				double h[], long nh)
{
long k, j, i, i1, i2;
long kn, nexp, n;
	
	nexp = NxIn*2L;
	n = NxIn;
	kn = n-1;
    
	if (nh < 2L) {	
		for (i=0L; i<NxIn; i++) {   
			j = i*2L;
			Out[j] = In[i];
			Out[j+1] = In[i];
		}
	}
	     
	else {        
		for (i=0L; i<nexp; i++) {    
			Out[i] = 0.0; 
			for (k=(i % 2L) ;k<nh; k+=2L) {  
				i1 = (i-k)/2L;
				if (i1 < 0L) {   
					i1 = (-i1) % kn;
					if (i1 > kn) 
						i1=kn-i1;
				}
				Out[i] = Out[i] + h[k]*In[i1];  
			}
			 
			for (k=2L-(i % 2L); k<nh; k+=2L) { 
				i2 = (i+k)/2L;
				if (i2 > kn) {
					i2 = i2 % kn;
				    i2 = kn-i2;
				    if (i2 > kn) 
				    	i2 = kn-i2;
				}
				Out[i] = Out[i] + h[k]*In[i2]; 
			}
		}
	}

}	

/* ----------------------------------------------------------------------------

	Function:	ReduceCentered_1D
		
	Purpose:	Reduces an image by a factor of two
				The reduced image grid is between the finer grid

	Parameters:
		In[NxIn] is the input signal	(NxIn should be greater than 2 and even)
		Out[NxIn/2] is the output signal
		g[ng] is an array that contains the coefficients of the filter		
				
	Author:		
		Michael Unser, NIH, BEIP, June 1994
		Patrick Brigger, NIH, BEIP,	May 1996, modified
		Daniel Sage, EPFL, Biomedical Imaging Group, April 1999

---------------------------------------------------------------------------- */
extern void ReduceCentered_1D(	double In[], long NxIn,
 								double Out[],
 								double g[], long ng)
{
double 	*y_tmp;
long 	k, i, i1, i2;
long	kk, kn, nred, n;
	 
	nred = NxIn/2L;
	n = nred*2L;
	kn = 2L*n;
	
	/* --- Allocated memory for a temporary buffer --- */
	y_tmp = (double *)malloc( (size_t)(n*(long)sizeof(double)));
	if ( y_tmp == (double *)NULL) {
		MessageDisplay("Out of memory in reduce_centered!");
		return;
	}
	
	/* --- Apply the symmetric filter to all the coefficients --- */
	for (k=0L; k<n; k++) {
		y_tmp[k] = In[k]*g[0 ];
		for (i=1L; i<ng; i++) {
			i1 = k-i;
			i2 = k+i;
			if (i1 < 0L) {
				i1 = (2L*n-1-i1) % kn;
				if (i1 >= n) 
					i1 = kn-i1-1L;
			}
			if (i2 > (n-1L)) {
				i2 = i2 % kn;
				if (i2 >= n) i2 = kn-i2-1L;
			}
			y_tmp[k] = y_tmp[k] + g[i]*(In[i1]+In[i2]);
		}
	}
	
	/* --- Now apply the Haar and perform downsampling --- */
	for(kk=0L; kk<nred; kk++) {
		k = 2L*kk;
		Out[kk] = (y_tmp[k] + y_tmp[k+1])/2.;
	}

	/* --- Free allocated memory --- */
	free(y_tmp);
	
}

/* ----------------------------------------------------------------------------

	Function:	
		ExpandCentered_1D
	
	Purpose:	
		Expands an 1D signal by a factor of two using the specified 
		FIR interpolation filter. The expansion is based on 
		a subsampled grid that is centered with respect to the
		finer grid

	Parameters:
		In[NxIn] is the input signal	(NxIn should be greater than 1)
		Out[NxIn"*2] is the output signal
		h[nh] is an array that contains the coefficients of the filter		
				
	Author:		
		Michael Unser, NIH, BEIP, June 1994
		Patrick Brigger, NIH, BEIP,	May 1996, modified
		Daniel Sage, EPFL, Biomedical Imaging Group, April 1999

---------------------------------------------------------------------------- */
extern void ExpandCentered_1D(	double In[], long NxIn,
								double Out[],
								double h[], long nh)
{
long k, i, j, i1, k0, i2;
long kk, kn, nexp, n;

	nexp = NxIn*2L;
	n = NxIn;
	kn = 2L*n;
	k0 = (nh/2L)*2-1L;
	     
	for (i=0L; i<NxIn; i++) {
		j = i*2L;
		Out[j] = In[i]*h[0];
		for (k=2L; k<nh; k=k+2L) {
			i1 = i-k/2L;
			i2 = i+k/2L;
			if (i1 < 0L) {				/* Provide the correct border conditions. */
				i1 = (2L*n-1-i1) % kn;	/* --> pseudo mirror image                */
				if (i1 >= n) 
					i1=kn-i1-1L;
			}
			if (i2 >= n) {				
				i2= (i2) % kn;
				if (i2 >= n) 
					i2=kn-i2-1L;
			}
			Out[j] = Out[j] + h[k]*(In[i1]+In[i2]);
		}
		Out[j+1] = 0.;  
		for (k=-k0; k<nh; k=k+2L) {
			kk=labs(k);    			/* filter coeff. are shifted with respect to above. */
			i1 = i+(k+1L)/2L;
			if (i1 < 0L) {
				i1 = (2*n-1-i1) % kn;
				if (i1 > n-1) 
					i1 = kn-i1-1L;
			}
			if (i1 >= n) {				
				i1 = (i1) % kn;
				if (i1 >= n) 
					i1=kn-i1-1;
			}
			Out[j+1L] = Out[j+1L] + h[kk]*In[i1];
		}
	}
	
	/* Now apply the Haar[-x] and  */
	for (j=nexp-1L; j>0L; j--)
		Out[j] = (Out[j] + Out[j-1])/2.0;
	Out[0] /= 2.0;
		
}

/* ----------------------------------------------------------------------------

	Function:	
		GetRow
		
	Purpose:	
		Get a row from an image

	Parameters:
		Image[Nx*Ny] is the input image
		RowNb is the number of the row
		Row[RowSize] is the output signal
				
---------------------------------------------------------------------------- */
static void GetRow( 
				float *Image, long Nx, long Ny,
				long RowNb,
				double *Row, long RowSize
				)
{
int i;
int BaseIndex;

	(void) Ny;
	BaseIndex = RowNb*Nx;
	for (i=0L; i<RowSize; i++) 
		Row[i] = (double)Image[BaseIndex+i];
}

/* ----------------------------------------------------------------------------

	Function:	
		GetColumn
		
	Purpose:	
		Get a column from an image

	Parameters:
		Image[Nx*Ny] is the input image
		ColumnNb is the number of the row
		Column[ColumnSize] is the output signal
				
---------------------------------------------------------------------------- */
static void GetColumn( 
				float *Image, long Nx, long Ny,
				long ColumnNb,
				double *Column, long ColumnSize
				)
{
int j;
int Index;

	(void) Ny;
	Index = ColumnNb;
	for (j=0L; j<ColumnSize; j++) {
		Column[j] = (double)Image[Index];
		Index += Nx;
	}
}

/* ----------------------------------------------------------------------------

	Function:	
		PutRow
		
	Purpose:	
		Put a row to an image

	Parameters:
		Image[Nx*Ny] is the input image
		RowNb is the number of the row
		Row[RowSize] is the output signal
				
---------------------------------------------------------------------------- */
static void PutRow( 
				float *Image, long Nx, long Ny,
				long RowNb,
				double *Row, long RowSize
				)
{
int i;
int BaseIndex;

	(void) Ny;
	BaseIndex = RowNb*Nx;
	for (i=0L; i<RowSize; i++) 
		Image[BaseIndex+i] = (float)Row[i];
}

/* ----------------------------------------------------------------------------

	Function:	
		PutColumn
		
	Purpose:	
		Put a column to an image

	Parameters:
		Image[Nx*Ny] is the input image
		ColumnNb is the number of the row
		Column[ColumnSize] is the output signal 
				
---------------------------------------------------------------------------- */
static void PutColumn( 
				float *Image, long Nx, long Ny,
				long ColumnNb, 
				double *Column, long ColumnSize
				)
{
int j;
int Index;

	(void) Ny;
	Index = ColumnNb;
	for (j=0L; j<ColumnSize; j++) {
		Image[Index] = (float)Column[j];
		Index += Nx;
	}
}
