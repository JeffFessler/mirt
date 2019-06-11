/*******************************************************************************

        Batch transpose interpolation for all locations of points
        with mirror and zero end condition

        Created by Se Young Chun, Jan 19, 2007, the University of Michigan

*******************************************************************************/

#include	<math.h>
#include	<stddef.h>
#include	<stdio.h>
#include	<stdlib.h>

#include	"batchinterpol2Dtran.h"
#include	"interpol2Dtran.h"
 
extern void     BatchInterpolatedValue2DTranMirr (
                        float   *Bout, 
                        float   *Bcoeff,
                        long    Width, 
                        long    Height,
                        double  nx,    
                        double  offx, 
                        double  mx,  
                        double  ny,    
                        double  offy, 
                        double  my,  
                        long    SplineDegree 
		)
{ /* begin BatchInterpolatedValue2DTranMirr */

	double	x, y;
	long	idx;
	long	idcoordx, idcoordy;

	for (idx = 0; idx < nx*ny; idx++) Bout[idx] = 0;

	idx = 0;	
	for (idcoordy = 0; idcoordy < Height; idcoordy++) {
		for (idcoordx = 0; idcoordx < Width; idcoordx++) {
			
			x = (idcoordx - offx) / mx;
			y = (idcoordy - offy) / my;

			InterpolatedValueTranMirr(Bout, nx, ny, 
				Bcoeff[idx++], x, y, SplineDegree);
		}
	}

} /* end BatchInterpolatedValue2DTranMirr */

extern void     BatchInterpolatedValue2DTranMirrWarp (
                        float   *Bout,      
                        float   *Bcoeff,   
                        long    Width,    
                        long    Height,  
                        double  nx,     
                        double  offx,  
                        double  mx,   
                        float   *WarpX,
                        double  ny,   
                        double  offy, 
                        double  my,   
                        float   *WarpY,
                        long    SplineDegree 
                )
{ /* begin BatchInterpolatedValue2DTranMirrWarp */

        double  x, y;
        long    idx;
        long    idcoordx, idcoordy;

        for (idx = 0; idx < nx*ny; idx++) Bout[idx] = 0;

        idx = 0;
        for (idcoordy = 0; idcoordy < Height; idcoordy++) {
                for (idcoordx = 0; idcoordx < Width; idcoordx++) {

                        x = (idcoordx + WarpX[idx] - offx) / mx;
                        y = (idcoordy + WarpY[idx] - offy) / my;

                        InterpolatedValueTranMirr(Bout, nx, ny, 
				Bcoeff[idx++], x, y, SplineDegree);
                }
        }
} /* end BatchInterpolatedValue2DTranMirrWarp */

extern void     BatchInterpolatedGradX2DTranMirr (
                        float   *Bout,      
                        float   *Bcoeff,   
                        long    Width,    
                        long    Height,  
                        double  nx,     
                        double  offx,  
                        double  mx,   
                        double  ny,  
                        double  offy,   
                        double  my,    
                        long    SplineDegree 
		)
{ /* begin BatchInterpolatedGradX2DTranMirr */

	double	x, y;
	long	idx;
	long	idcoordx, idcoordy;

	for (idx = 0; idx < nx*ny; idx++) Bout[idx] = 0;

	idx = 0;	
	for (idcoordy = 0; idcoordy < Height; idcoordy++) {
		for (idcoordx = 0; idcoordx < Width; idcoordx++) {
			
			x = (idcoordx - offx) / mx;
			y = (idcoordy - offy) / my;

			InterpolatedGradXTranMirr(Bout, nx, ny, 
				Bcoeff[idx++], x, y, SplineDegree);
		}
	}

	for (idx = 0; idx < nx*ny; idx++) Bout[idx] = Bout[idx] / mx;

} /* end BatchInterpolatedGradX2DTranMirr */

extern void     BatchInterpolatedGradX2DTranMirrWarp (
                        float   *Bout,      
                        float   *Bcoeff,   
                        long    Width,    
                        long    Height,  
                        double  nx,     
                        double  offx,  
                        double  mx,   
                        float   *WarpX,  
                        double  ny,     
                        double  offy,  
                        double  my,   
                        float   *WarpY, 
                        long    SplineDegree
                )
{ /* begin BatchInterpolatedGradX2DTranMirrWarp */

        double  x, y;
        long    idx;
        long    idcoordx, idcoordy;

        for (idx = 0; idx < nx*ny; idx++) Bout[idx] = 0;

        idx = 0;
        for (idcoordy = 0; idcoordy < Height; idcoordy++) {
                for (idcoordx = 0; idcoordx < Width; idcoordx++) {

                        x = (idcoordx + WarpX[idx] - offx) / mx;
                        y = (idcoordy + WarpY[idx] - offy) / my;

                        InterpolatedGradXTranMirr(Bout, nx, ny, 
				Bcoeff[idx++], x, y, SplineDegree);
                }
        }

        for (idx = 0; idx < nx*ny; idx++) Bout[idx] = Bout[idx] / mx;

} /* end BatchInterpolatedGradX2DTranMirrWarp */


extern void     BatchInterpolatedGradY2DTranMirr (
                        float   *Bout,      
                        float   *Bcoeff,   
                        long    Width,    
                        long    Height,  
                        double  nx,     
                        double  offx,  
                        double  mx,   
                        double  ny,   
                        double  offy, 
                        double  my,   
                        long    SplineDegree 
                )
{ /* begin BatchInterpolatedGradY2DTranMirr */

	double	x, y;
	long	idx;
	long	idcoordx, idcoordy;

	for (idx = 0; idx < nx*ny; idx++) Bout[idx] = 0;

	idx = 0;	
	for (idcoordy = 0; idcoordy < Height; idcoordy++) {
		for (idcoordx = 0; idcoordx < Width; idcoordx++) {
			
			x = (idcoordx - offx) / mx;
			y = (idcoordy - offy) / my;

			InterpolatedGradYTranMirr(Bout, nx, ny, 
				Bcoeff[idx++], x, y, SplineDegree);
		}
	}

	for (idx = 0; idx < nx*ny; idx++) Bout[idx] = Bout[idx] / my;

} /* end BatchInterpolatedGradY2DTranMirr */

extern void     BatchInterpolatedGradY2DTranMirrWarp (
                        float   *Bout,      
                        float   *Bcoeff,   
                        long    Width,    
                        long    Height,  
                        double  nx,     
                        double  offx,  
                        double  mx,  
                        float   *WarpX,  
                        double  ny,     
                        double  offy,  
                        double  my,   
                        float   *WarpY, 
                        long    SplineDegree
                )
{ /* begin BatchInterpolatedGradY2DTranMirrWarp */

        double  x, y;
        long    idx;
        long    idcoordx, idcoordy;

        for (idx = 0; idx < nx*ny; idx++) Bout[idx] = 0;

        idx = 0;
        for (idcoordy = 0; idcoordy < Height; idcoordy++) {
                for (idcoordx = 0; idcoordx < Width; idcoordx++) {

                        x = (idcoordx + WarpX[idx] - offx) / mx;
                        y = (idcoordy + WarpY[idx] - offy) / my;

                        InterpolatedGradYTranMirr(Bout, nx, ny, 
				Bcoeff[idx++], x, y, SplineDegree);
                }
        }

        for (idx = 0; idx < nx*ny; idx++) Bout[idx] = Bout[idx] / my;

} /* end BatchInterpolatedGradY2DTranMirrWarp */

extern void     BatchInterpolatedValue2DTranZero (
                        float   *Bout, 
                        float   *Bcoeff,
                        long    Width, 
                        long    Height,
                        double  nx,    
                        double  offx, 
                        double  mx,  
                        double  ny,    
                        double  offy, 
                        double  my,  
                        long    SplineDegree 
		)
{ /* begin BatchInterpolatedValue2DTranZero */

	double	x, y;
	long	idx;
	long	idcoordx, idcoordy;

	for (idx = 0; idx < nx*ny; idx++) Bout[idx] = 0;

	idx = 0;	
	for (idcoordy = 0; idcoordy < Height; idcoordy++) {
		for (idcoordx = 0; idcoordx < Width; idcoordx++) {
			
			x = (idcoordx - offx) / mx;
			y = (idcoordy - offy) / my;

			InterpolatedValueTranZero(Bout, nx, ny, 
				Bcoeff[idx++], x, y, SplineDegree);
		}
	}

} /* end BatchInterpolatedValue2DTranZero */

extern void     BatchInterpolatedValue2DTranZeroWarp (
                        float   *Bout,      
                        float   *Bcoeff,   
                        long    Width,    
                        long    Height,  
                        double  nx,     
                        double  offx,  
                        double  mx,   
                        float   *WarpX,
                        double  ny,   
                        double  offy, 
                        double  my,   
                        float   *WarpY,
                        long    SplineDegree 
                )
{ /* begin BatchInterpolatedValue2DTranZeroWarp */

        double  x, y;
        long    idx;
        long    idcoordx, idcoordy;

        for (idx = 0; idx < nx*ny; idx++) Bout[idx] = 0;

        idx = 0;
        for (idcoordy = 0; idcoordy < Height; idcoordy++) {
                for (idcoordx = 0; idcoordx < Width; idcoordx++) {

                        x = (idcoordx + WarpX[idx] - offx) / mx;
                        y = (idcoordy + WarpY[idx] - offy) / my;

                        InterpolatedValueTranZero(Bout, nx, ny, 
				Bcoeff[idx++], x, y, SplineDegree);
                }
        }
} /* end BatchInterpolatedValue2DTranZeroWarp */

extern void     BatchInterpolatedGradX2DTranZero (
                        float   *Bout,      
                        float   *Bcoeff,   
                        long    Width,    
                        long    Height,  
                        double  nx,     
                        double  offx,  
                        double  mx,   
                        double  ny,  
                        double  offy,   
                        double  my,    
                        long    SplineDegree 
		)
{ /* begin BatchInterpolatedGradX2DTranZero */

	double	x, y;
	long	idx;
	long	idcoordx, idcoordy;

	for (idx = 0; idx < nx*ny; idx++) Bout[idx] = 0;

	idx = 0;	
	for (idcoordy = 0; idcoordy < Height; idcoordy++) {
		for (idcoordx = 0; idcoordx < Width; idcoordx++) {
			
			x = (idcoordx - offx) / mx;
			y = (idcoordy - offy) / my;

			InterpolatedGradXTranZero(Bout, nx, ny, 
				Bcoeff[idx++], x, y, SplineDegree);
		}
	}

	for (idx = 0; idx < nx*ny; idx++) Bout[idx] = Bout[idx] / mx;

} /* end BatchInterpolatedGradX2DTranZero */

extern void     BatchInterpolatedGradX2DTranZeroWarp (
                        float   *Bout,      
                        float   *Bcoeff,   
                        long    Width,    
                        long    Height,  
                        double  nx,     
                        double  offx,  
                        double  mx,   
                        float   *WarpX,  
                        double  ny,     
                        double  offy,  
                        double  my,   
                        float   *WarpY, 
                        long    SplineDegree
                )
{ /* begin BatchInterpolatedGradX2DTranZeroWarp */

        double  x, y;
        long    idx;
        long    idcoordx, idcoordy;

        for (idx = 0; idx < nx*ny; idx++) Bout[idx] = 0;

        idx = 0;
        for (idcoordy = 0; idcoordy < Height; idcoordy++) {
                for (idcoordx = 0; idcoordx < Width; idcoordx++) {

                        x = (idcoordx + WarpX[idx] - offx) / mx;
                        y = (idcoordy + WarpY[idx] - offy) / my;

                        InterpolatedGradXTranZero(Bout, nx, ny, 
				Bcoeff[idx++], x, y, SplineDegree);
                }
        }

        for (idx = 0; idx < nx*ny; idx++) Bout[idx] = Bout[idx] / mx;

} /* end BatchInterpolatedGradX2DTranZeroWarp */


extern void     BatchInterpolatedGradY2DTranZero (
                        float   *Bout,      
                        float   *Bcoeff,   
                        long    Width,    
                        long    Height,  
                        double  nx,     
                        double  offx,  
                        double  mx,   
                        double  ny,   
                        double  offy, 
                        double  my,   
                        long    SplineDegree 
                )
{ /* begin BatchInterpolatedGradY2DTranZero */

	double	x, y;
	long	idx;
	long	idcoordx, idcoordy;

	for (idx = 0; idx < nx*ny; idx++) Bout[idx] = 0;

	idx = 0;	
	for (idcoordy = 0; idcoordy < Height; idcoordy++) {
		for (idcoordx = 0; idcoordx < Width; idcoordx++) {
			
			x = (idcoordx - offx) / mx;
			y = (idcoordy - offy) / my;

			InterpolatedGradYTranZero(Bout, nx, ny, 
				Bcoeff[idx++], x, y, SplineDegree);
		}
	}

	for (idx = 0; idx < nx*ny; idx++) Bout[idx] = Bout[idx] / my;

} /* end BatchInterpolatedGradY2DTranZero */

extern void     BatchInterpolatedGradY2DTranZeroWarp (
                        float   *Bout,      
                        float   *Bcoeff,   
                        long    Width,    
                        long    Height,  
                        double  nx,     
                        double  offx,  
                        double  mx,  
                        float   *WarpX,  
                        double  ny,     
                        double  offy,  
                        double  my,   
                        float   *WarpY, 
                        long    SplineDegree
                )
{ /* begin BatchInterpolatedGradY2DTranZeroWarp */

        double  x, y;
        long    idx;
        long    idcoordx, idcoordy;

        for (idx = 0; idx < nx*ny; idx++) Bout[idx] = 0;

        idx = 0;
        for (idcoordy = 0; idcoordy < Height; idcoordy++) {
                for (idcoordx = 0; idcoordx < Width; idcoordx++) {

                        x = (idcoordx + WarpX[idx] - offx) / mx;
                        y = (idcoordy + WarpY[idx] - offy) / my;

                        InterpolatedGradYTranZero(Bout, nx, ny, 
				Bcoeff[idx++], x, y, SplineDegree);
                }
        }

        for (idx = 0; idx < nx*ny; idx++) Bout[idx] = Bout[idx] / my;

} /* end BatchInterpolatedGradY2DTranZeroWarp */
