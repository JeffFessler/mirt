/*******************************************************************************

        Transpose interpolation of grid locations of points using filtering 
        technique (kernel) with zero boundary condition

        Created by Se Young Chun, Nov 2, 2006, the University of Michigan

*******************************************************************************/
#include	<math.h>
#include	<stddef.h>
#include	<stdio.h>
#include	<stdlib.h>

#include	"filtinterpol2Dtran.h"
#include	"BIG/BsplnTrfKer.h"
 
extern void     FiltInterpolatedValue2DTranZero (
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
{ /* begin FiltInterpolatedValue2DTranZero */

	double	*Btmp;
	int	idcoordx, idcoordy;

	Btmp = (double*)malloc(Width*Height*sizeof(double));

	inverseBsplineKernel3Finite(Bcoeff, Btmp, (long)Width, (long)Height, 
		(long)1, kernelx, (long)nkx, kernely, (long)nky, NULL, 0);

	for (idcoordy = 0; idcoordy < ny; idcoordy++) {
		for (idcoordx = 0; idcoordx < nx; idcoordx++) {

			Bout[idcoordx + nx*idcoordy] = Btmp[offx+idcoordx*mx 
						+ Width*(offy+idcoordy*my)];

		}
	}

	free(Btmp);

} /* end FiltInterpolatedValue2DTranZero */
