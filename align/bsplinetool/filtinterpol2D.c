/*******************************************************************************
        Interpolation of grid locations of points using filtering technique
        (kernel) with zero boundary condition

        Created by Se Young Chun, Oct 26, 2006, the University of Michigan

*******************************************************************************/

#include	<math.h>
#include	<stddef.h>
#include	<stdio.h>
#include	<stdlib.h>

#include	"filtinterpol2D.h"
#include	"BIG/BsplnTrfKer.h"
 
extern void     FiltInterpolatedValue2DZero (
                        double  *Bout,  
                        double  *Bcoeff,
                        int     Width,  
                        int     Height, 
                        int     nx,    
                        int     offx, 
                        int     mx, 
                        int     ny,   
                        int     offy,
                        int     my, 
                        double  *kernelx, 
                        int     nkx,     
                        double  *kernely, 
                        int     nky      
		)
{ /* begin FiltInterpolatedValue2DZero */

	double	*Bin;
	int	idcoordx, idcoordy, nz = 1;

	Bin = (double*)calloc(nx*ny, sizeof(double));

	for (idcoordy = 0; idcoordy < Height; idcoordy++) {
		for (idcoordx = 0; idcoordx < Width; idcoordx++) {

			Bin[offx+idcoordx*mx + nx*(offy+idcoordy*my)] 
				= Bcoeff[idcoordx + Width*idcoordy];

		}
	}

	inverseBsplineKernel3Finite(Bin, Bout, (long)nx, (long)ny, (long)nz, 
		kernelx, (long)nkx, kernely, (long)nky, NULL, 0);

	free(Bin);

} /* end FiltInterpolatedValue2DZero */
