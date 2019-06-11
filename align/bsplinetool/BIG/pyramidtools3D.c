/*******************************************************************************        Trivial extension of pyramidtools.h in spline pyramids package at
        http://bigwww.epfl.ch/algorithms.html to 3D reduce and expand

        Coded by Se Young Chun, May 25, 2007, The Univ. of Michigan

*******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include <stddef.h>

/* available at http://bigwww.epfl.ch/sage/pyramids/ */
#include "pyramids/configs.h"
#include "pyramids/messagedisplay.h"
#include "pyramids/pyramidtools.h"
#include "pyramids/pyramidfilters.h"

#include "pyramidtools3D.h"


static void 	GetColumn3D (
                	float *Image, long x, long z,         
                	double Line[],     
                	long Width, long Height   
                );

static void 	GetRow3D (
                	float *Image, long y, long z,     
                	double Line[], 
                	long Width, long Height
                );

static void 	GetSlice3D (
                	float *Image, long x, long y,      
                	double Line[], 
                	long Width, long Height, long Slice
                );

static void 	PutColumn3D (
                	float *Image, long x, long z,          
                	double Line[],      
                	long Width, long Height         
                );

static void 	PutRow3D (
                	float *Image, long y, long z,          
                	double Line[],    
                	long Width, long Height
                );

static void 	PutSlice3D (
                	float *Image, long x, long y,          
                	double Line[],    
                	long Width, long Height, long Slice
                );



extern int 	Reduce_3D (	
			float *In, long NxIn, long NyIn, long NzIn,
			float *Out,
			double g[], long ng,
			short IsCentered
		)
{
	float	*Tmp, *Tmp2;
	double 	*InBuffer;
	double	*OutBuffer;
	long	kx, ky, kz;
	long 	NxOut;
	long 	NyOut;
	long 	NzOut;

	/* --- Define dimension of the output --- */
	NxOut = NxIn/2L;
	if (NxOut < 1L) NxOut = 1L;
	
	NyOut = NyIn/2L;
	if (NyOut < 1L) NyOut = 1L;

	NzOut = NzIn/2L;
	if (NzOut < 1L) NzOut = 1L;

	/* --- Allocate a temporary image --- */
	Tmp = (float *)malloc((size_t)(NxOut*NyIn*NzIn*(long)sizeof(float)));
	if (Tmp == (float *)NULL) {
		MessageDisplay("Unable to allocate memory");
		return(ERROR);
	}

	Tmp2 = (float *)malloc((size_t)(NxOut*NyOut*NzIn*(long)sizeof(float)));
	if (Tmp2 == (float *)NULL) {
		free(Tmp);
		MessageDisplay("Unable to allocate memory");
		return(ERROR);
	}
	
	/* --- X processing --- */
	if (NxIn > 1L) {
		InBuffer=(double*)malloc((size_t)(NxIn*(long)sizeof(double)));
		if (InBuffer == (double *)NULL) {
			free(Tmp); free(Tmp2);
			MessageDisplay("Unable to allocate memory");
			return(ERROR);
		}
		OutBuffer=(double*)malloc((size_t)(NxOut*(long)sizeof(double)));
		if (OutBuffer == (double *)NULL) {
			free(Tmp); free(Tmp2); free(InBuffer);
			MessageDisplay("Unable to allocate memory");
			return(ERROR);
		}
		for (kz=0L; kz<NzIn; kz++) {
			for (ky=0L; ky<NyIn; ky++) {
				GetRow3D(In, ky, kz, InBuffer, NxIn, NyIn);
				Reduce_1D(InBuffer, NxIn, OutBuffer, g, ng, 
					IsCentered); 
 				PutRow3D(Tmp, ky, kz, OutBuffer, NxOut, NyIn);
    			}
		}
		free(InBuffer);
		free(OutBuffer);
    	}
    	else
	   	memcpy( Tmp, In, (size_t)(NyIn*NzIn*(long)sizeof(float)));
    	
	/* --- Y processing --- */
	if (NyIn > 1L) {
		InBuffer=(double*)malloc((size_t)(NyIn*(long)sizeof(double)));
		if (InBuffer == (double *)NULL) {
			MessageDisplay("Unable to allocate memory");
			free(Tmp); free(Tmp2);
			return(ERROR);
		}
		OutBuffer=(double*)malloc((size_t)(NyOut*(long)sizeof(double)));
		if (OutBuffer == (double *)NULL) {
			free(InBuffer); free(Tmp); free(Tmp2);
			MessageDisplay("Unable to allocate memory");
			return(ERROR);
		}
		for (kz=0L; kz<NzIn; kz++) {
	    		for (kx=0L; kx<NxOut; kx++) {
				GetColumn3D(Tmp, kx, kz, InBuffer, NxOut, NyIn);
				Reduce_1D(InBuffer, NyIn, OutBuffer, g, ng,
					IsCentered);
	 	   		PutColumn3D(Tmp2,kx,kz,OutBuffer,NxOut,NyOut);
			}
   		}
		free(InBuffer);
		free(OutBuffer);
	}
	else
	   	memcpy( Tmp2, Tmp, (size_t)(NxOut*NzIn*(long)sizeof(float)));
	
	/* --- Z processing --- */
	if (NzIn > 1L) {
		InBuffer=(double*)malloc((size_t)(NzIn*(long)sizeof(double)));
		if (InBuffer == (double *)NULL) {
			MessageDisplay("Unable to allocate memory");
			free(Tmp); free(Tmp2);
			return(ERROR);
		}
		OutBuffer=(double*)malloc((size_t)(NzOut*(long)sizeof(double)));
		if (OutBuffer == (double *)NULL) {
			free(InBuffer); free(Tmp); free(Tmp2);
			MessageDisplay("Unable to allocate memory");
			return(ERROR);
		}
	    	for (ky=0L; ky<NyOut; ky++) {
	    		for (kx=0L; kx<NxOut; kx++) {
				GetSlice3D(Tmp2,kx,ky,InBuffer,NxOut,NyOut
					,NzIn);
				Reduce_1D(InBuffer, NzIn, OutBuffer, g, ng
					,IsCentered);
	 	   		PutSlice3D(Out,kx,ky,OutBuffer,NxOut,NyOut
					,NzOut);
   			}
   		}
		free(InBuffer);
		free(OutBuffer);
	}
	else
	   	memcpy( Out, Tmp2, (size_t)(NxOut*NyOut*(long)sizeof(float)));
	
	/* --- Free the temporary image --- */
	free(Tmp); free(Tmp2);
	
	return (!ERROR);
}



extern int 	Expand_3D (	
			float *In, long NxIn, long NyIn, long NzIn,
			float *Out,
			double h[], long nh,
			short IsCentered
		)
{
	float	*Tmp, *Tmp2;
	double  *InBuffer; 	/* Input buffer to 1D process */ 
	double  *OutBuffer;	/* Output buffer to 1D process */ 
	long 	kx, ky, kz;
	long 	NxOut;
	long 	NyOut;
	long 	NzOut;
	
	if (NxIn <= 1L) 
		NxOut = 1L; 
	else 
		NxOut = NxIn*2L;
		
	if (NyIn <= 1L) 
		NyOut = 1L; 
	else 
		NyOut = NyIn*2L;
 
	if (NzIn <= 1L) 
		NzOut = 1L; 
	else 
		NzOut = NzIn*2L;
 
	/* --- Allocate a temporary image --- */
	Tmp = (float *)malloc((size_t)(NxOut*NyIn*NzIn*(long)sizeof(float)));
	if (Tmp == (float *)NULL) {
		MessageDisplay("Unable to allocate memory");
		return(ERROR);
	}

	Tmp2 = (float *)malloc((size_t)(NxOut*NyOut*NzIn*(long)sizeof(float)));
	if (Tmp2 == (float *)NULL) {
		free(Tmp);
		MessageDisplay("Unable to allocate memory");
		return(ERROR);
	}
	
	/* --- X processing --- */
	if (NxIn > 1L) {
		InBuffer=(double*)malloc((size_t)(NxIn*(long)sizeof(double)));
		if (InBuffer == (double *)NULL) {
			MessageDisplay("Unable to allocate memory");
			free(Tmp); free(Tmp2);
			return(ERROR);
		}
		OutBuffer=(double*)malloc((size_t)(NxOut*(long)sizeof(double)));
		if (OutBuffer == (double *)NULL) {
			free(InBuffer); free(Tmp); free(Tmp2);
			MessageDisplay("Unable to allocate memory");
			return(ERROR);
		}
		for (kz=0L; kz<NzIn; kz++) {
			for (ky=0L; ky<NyIn; ky++) {
				GetRow3D(In, ky, kz, InBuffer, NxIn, NyIn);
				Expand_1D(InBuffer, NxIn, OutBuffer, h, nh
					, IsCentered);
 				PutRow3D(Tmp, ky, kz, OutBuffer, NxOut, NyIn);
    			}
		}
		free(InBuffer);
		free(OutBuffer);
	}
	else
	   	memcpy(Tmp, In, (size_t)(NxIn*NyIn*NzIn*(long)sizeof(float)));

	/* --- Y processing --- */
	if (NyIn > 1L) {
		InBuffer=(double*)malloc((size_t)(NyIn*(long)sizeof(double)));
		if (InBuffer == (double *)NULL) {
			MessageDisplay("Unable to allocate memory");
			free(Tmp); free(Tmp2);
			return(ERROR);
		}
		OutBuffer=(double*)malloc((size_t)(NyOut*(long)sizeof(double)));
		if (OutBuffer == (double *)NULL) {
			free(InBuffer); free(Tmp); free(Tmp2);
			MessageDisplay("Unable to allocate memory");
			return(ERROR);
		} 
		for (kz=0L; kz<NzIn; kz++) {
	    		for (kx=0L; kx<NxOut; kx++) {
				GetColumn3D(Tmp,kx,kz,InBuffer,NxOut,NyIn);
				Expand_1D(InBuffer, NyIn, OutBuffer, h, nh
					, IsCentered); 
	 	   		PutColumn3D(Tmp2,kx,kz,OutBuffer,NxOut,NyOut); 
			}
   		}
		free(InBuffer);
		free(OutBuffer);
	}
	else
	   	memcpy(Tmp2, Tmp, (size_t)(NxOut*NyIn*NzIn*(long)sizeof(float)));
	 
	/* --- Z processing --- */
	if (NzIn > 1L) {
		InBuffer=(double*)malloc((size_t)(NzIn*(long)sizeof(double)));
		if (InBuffer == (double *)NULL) {
			MessageDisplay("Unable to allocate memory");
			free(Tmp); free(Tmp2);
			return(ERROR);
		}
		OutBuffer=(double*)malloc((size_t)(NzOut*(long)sizeof(double)));
		if (OutBuffer == (double *)NULL) {
			free(InBuffer); free(Tmp); free(Tmp2);
			MessageDisplay("Unable to allocate memory");
			return(ERROR);
		} 
	    	for (ky=0L; ky<NyOut; ky++) {
	    		for (kx=0L; kx<NxOut; kx++) {
				GetSlice3D(Tmp2, kx, ky, InBuffer, NxOut, NyOut
					, NzIn);
				Expand_1D(InBuffer, NzIn, OutBuffer, h, nh
					, IsCentered); 
	 	   		PutSlice3D(Out, kx, ky, OutBuffer, NxOut, NyOut
					, NzOut);
   			}
   		}
		free(InBuffer);
		free(OutBuffer);
	}
	else
	   	memcpy(Out,Tmp2,(size_t)(NxOut*NyOut*NzIn*(long)sizeof(float)));


	free(Tmp); free(Tmp2);

	return (ERROR);
	
}



static void 	GetColumn3D (
                	float   *Image,     
                        long    x,      
                        long    z,         
                        double  Line[],     
                        long    Width,     
                        long    Height   
                )
{ /* begin GetColumn3D */
        long    y;
        Image = Image + (ptrdiff_t)(x + z * Width * Height);
        for (y = 0L; y < Height; y++) {
                Line[y] = (double)*Image;
                Image += (ptrdiff_t)Width;
        }
} /* end GetColumn3D */



static void 	GetRow3D (
                	float   *Image,  
                        long    y,        
                        long    z,     
                        double  Line[], 
                        long    Width,   
                        long    Height
                )
{ /* begin GetRow3D */
        long    x;
        Image = Image + (ptrdiff_t)(y * Width + z * Width * Height);
        for (x = 0L; x < Width; x++) {
                Line[x] = (double)*Image++;
        }
} /* end GetRow3D */



static void 	GetSlice3D (
                	float   *Image,   
                        long    x,       
                        long    y,      
                        double  Line[], 
                        long    Width,    
                        long    Height,
                        long    Slice
                )
{ /* begin GetSlice3D */
        long    z;
        Image = Image + (ptrdiff_t)(x + y * Width);
        for (z = 0L; z < Slice; z++) {
                Line[z] = (double)*Image;
                Image += (ptrdiff_t)(Width * Height);
        }
} /* end GetSlice3D */



static void 	PutColumn3D (
                	float   *Image,   
                        long    x,         
                        long    z,          
                        double  Line[],      
                        long    Width,        
                        long    Height         
                )
{ /* begin PutColumn3D */
        long    y;
        Image = Image + (ptrdiff_t)(x + z * Width * Height);
        for (y = 0L; y < Height; y++) {
                *Image = (float)Line[y];
                Image += (ptrdiff_t)Width;
        }
} /* end PutColumn3D */



static void 	PutRow3D (
                	float   *Image,       
                        long    y,           
                        long    z,          
                        double  Line[],    
                        long    Width,   
                        long    Height
                )
{ /* begin PutRow3D */
        long    x;
        Image = Image + (ptrdiff_t)(y * Width + z * Width * Height);
        for (x = 0L; x < Width; x++) {
                *Image++ = (float)Line[x];
        }
} /* end PutRow3D */



static void 	PutSlice3D (
                	float   *Image,       
                        long    x,           
                        long    y,          
                        double  Line[],    
                        long    Width,    
                        long    Height,
                        long    Slice
                )
{ /* begin PutSlice3D */
        long    z;
        Image = Image + (ptrdiff_t)(x + y * Width);
        for (z = 0L; z < Slice; z++) {
                *Image = (float)Line[z] ;
                Image += (ptrdiff_t)(Width * Height);
        }
} /* end PutSlice3D */
