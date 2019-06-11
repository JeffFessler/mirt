/* ----------------------------------------------------------------------------
	Filename:  	pyramidtools.h
	
	Project:	Biomedical Imaging Library
 
	Author:		Daniel Sage
				Swiss Federal Institute of Technology - Lausanne
				Biomedical Imaging Group
				EPFL/DMT/IOA, BM-Ecublens, CH-1015 Lausanne, Switzerland
				
	Date:		17 March 1999
	
	Purpose:	Header associated to pyramidtools.c
			 
---------------------------------------------------------------------------- */

extern int GetPyramidFilter(
				char *Filter, 
				long Order, 
				double g[],long *ng,
				double h[],long *nh, 
				short *FlagCentered		
				);	
					
extern int Reduce_2D(	
				float *In, long NxIn, long NyIn,
				float *Out,
				double w[], long nw,
				short FlagCentered
				);
				
extern int Expand_2D(	
				float *In, long NxIn, long NyIn,
				float *Out,
				double w[], long nw,
				short FlagCentered
				);

extern void Reduce_1D(	
				double x[], long nx,
				double y[],
				double w[], long nw,
				short FlagCentered
				);
						
extern void Expand_1D(	
				double x[], long nx,
				double y[],
				double w[], long nw,
				short FlagCentered
				);
						
