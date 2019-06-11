/*******************************************************************************

        Interpolation of grid locations of points using filtering technique
        (kernel) with zero boundary condition

        Created by Se Young Chun, Oct 20, 2006, the University of Michigan

*******************************************************************************/

#include	<math.h>
#include	<stddef.h>
#include	<stdio.h>
#include	<stdlib.h>

#include	"filtinterpol3D.h"
#include	"BIG/BsplnTrfKer.h"
 
extern void     FiltInterpolatedValue3DZero (
                        double  *Bout,   
                        double  *Bcoeff,
                        int     Width, 
                        int     Height,
                        int     Slice, 
                        int     nx,    
                        int     offx,  
                        int     mx,   
                        int     ny,  
                        int     offy,
                        int     my,  
                        int     nz,  
                        int     offz,
                        int     mz,  
                        double  *kernelx,
                        int     nkx,    
                        double  *kernely,
                        int     nky,    
                        double  *kernelz, 
                        int     nkz      
		)
{ /* begin FiltInterpolatedValue3DZero */

	double	*Bin;
	int	idcoordx, idcoordy, idcoordz;

	Bin = (double*)calloc(nx*ny*nz, sizeof(double));

	for (idcoordz = 0; idcoordz < Slice; idcoordz++) {
		for (idcoordy = 0; idcoordy < Height; idcoordy++) {
			for (idcoordx = 0; idcoordx < Width; idcoordx++) {

				Bin[offx+idcoordx*mx + nx*(offy+idcoordy*my) 
					+ nx*ny*(offz+idcoordz*mz)] 
					= Bcoeff[idcoordx + Width*idcoordy 
						+ Width*Height*idcoordz];

			}
		}
	}

	inverseBsplineKernel3Finite(Bin, Bout, (long)nx, (long)ny, (long)nz, 
		kernelx, (long)nkx, kernely, (long)nky, kernelz, (long)nkz);

	free(Bin);

} /* end FiltInterpolatedValue3DZero */
