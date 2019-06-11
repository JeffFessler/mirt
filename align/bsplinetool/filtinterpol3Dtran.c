/*******************************************************************************

        Transpose interpolation of grid locations of points using filtering   
      	technique (kernel) with zero boundary condition

        Created by Se Young Chun, Nov 2, 2006, the University of Michigan

*******************************************************************************/

#include	<math.h>
#include	<stddef.h>
#include	<stdio.h>
#include	<stdlib.h>

#include	"filtinterpol3Dtran.h"
#include	"BIG/BsplnTrfKer.h"
 
extern void     FiltInterpolatedValue3DTranZero (
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
{ /* begin FiltInterpolatedValue3DTranZero */

	double	*Btmp;
	int	idcoordx, idcoordy, idcoordz;

	Btmp = (double*)malloc(Width*Height*Slice*sizeof(double));

	inverseBsplineKernel3Finite(Bcoeff,Btmp,(long)Width,(long)Height,
	    (long)Slice,kernelx,(long)nkx,kernely,(long)nky,kernelz,(long)nkz);

	for (idcoordz = 0; idcoordz < nz; idcoordz++) {
		for (idcoordy = 0; idcoordy < ny; idcoordy++) {
			for (idcoordx = 0; idcoordx < nx; idcoordx++) {

				Bout[idcoordx + nx*idcoordy + nx*ny*idcoordz] = 
					Btmp[offx+idcoordx*mx 
					+ Width*(offy+idcoordy*my) 
					+ Height*Width*(offz+idcoordz*mz)];

			}
		}
	}

	free(Btmp);

} /* end FiltInterpolatedValue3DTranZero */
