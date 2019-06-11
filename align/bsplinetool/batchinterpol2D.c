/*******************************************************************************

        Batch interpolation for each location of points (+ deformation)
        using point interpolation functions with mirror and zero end condition

        Created by Se Young Chun, Sep 30, 2006, the University of Michigan

*******************************************************************************/

#include	<math.h>
#include	<stddef.h>
#include	<stdio.h>
#include	<stdlib.h>

#include	"batchinterpol2D.h"
#include	"interpol2D.h"
 
extern void     BatchInterpolatedValue2DMirr (
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
{ /* begin BatchInterpolatedValue2DMirr */

	double	x, y;
	long	idx;
	long	idcoordx, idcoordy;

	idx = 0;	
	for (idcoordy = 0; idcoordy < ny; idcoordy++) {
		for (idcoordx = 0; idcoordx < nx; idcoordx++) {
			
			x = (idcoordx - offx) / mx;
			y = (idcoordy - offy) / my;

			Bout[idx++] = InterpolatedValueMirr(Bcoeff, 
					Width, Height, x, y, SplineDegree);
		}
	}
} /* end BatchInterpolatedValue2DMirr */


extern void     BatchInterpolatedValue2DMirrWarp (
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
{ /* begin BatchInterpolatedValue2DMirrWarp */

        double  x, y;
        long    idx;
        long    idcoordx, idcoordy;

        idx = 0;
        for (idcoordy = 0; idcoordy < ny; idcoordy++) {
                for (idcoordx = 0; idcoordx < nx; idcoordx++) {

                        x = (idcoordx + WarpX[idx] - offx) / mx;
                        y = (idcoordy + WarpY[idx] - offy) / my;

                        Bout[idx++] = InterpolatedValueMirr(Bcoeff,
                                        Width, Height, x, y, SplineDegree);
                }
        }
} /* end BatchInterpolatedValue2DMirrWarp */


extern void     BatchInterpolatedGradX2DMirr (
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
{ /* begin BatchInterpolatedGradX2DMirr */

	double	x, y;
	long	idx;
	long	idcoordx, idcoordy;

	idx = 0;	
	for (idcoordy = 0; idcoordy < ny; idcoordy++) {
		for (idcoordx = 0; idcoordx < nx; idcoordx++) {
			
			x = (idcoordx - offx) / mx;
			y = (idcoordy - offy) / my;

			Bout[idx++] = InterpolatedGradXMirr(Bcoeff, 
					Width, Height, x, y, SplineDegree)/mx;
		}
	}
} /* end BatchInterpolatedGradX2DMirr */

extern void     BatchInterpolatedGradX2DMirrWarp (
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
{ /* begin BatchInterpolatedGradX2DMirrWarp */

        double  x, y;
        long    idx;
        long    idcoordx, idcoordy;

        idx = 0;
        for (idcoordy = 0; idcoordy < ny; idcoordy++) {
                for (idcoordx = 0; idcoordx < nx; idcoordx++) {

                        x = (idcoordx + WarpX[idx] - offx) / mx;
                        y = (idcoordy + WarpY[idx] - offy) / my;

                        Bout[idx++] = InterpolatedGradXMirr(Bcoeff,
                                        Width, Height, x, y, SplineDegree)/mx;
                }
        }
} /* end BatchInterpolatedGradX2DMirrWarp */


extern void     BatchInterpolatedGradY2DMirr (
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
{ /* begin BatchInterpolatedGradY2DMirr */

	double	x, y;
	long	idx;
	long	idcoordx, idcoordy;

	idx = 0;	
	for (idcoordy = 0; idcoordy < ny; idcoordy++) {
		for (idcoordx = 0; idcoordx < nx; idcoordx++) {
			
			x = (idcoordx - offx) / mx;
			y = (idcoordy - offy) / my;

			Bout[idx++] = InterpolatedGradYMirr(Bcoeff, 
					Width, Height, x, y, SplineDegree)/my;
		}
	}
} /* end BatchInterpolatedGradY2DMirr */

extern void     BatchInterpolatedGradY2DMirrWarp (
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
{ /* begin BatchInterpolatedGradY2DMirrWarp */

        double  x, y;
        long    idx;
        long    idcoordx, idcoordy;

        idx = 0;
        for (idcoordy = 0; idcoordy < ny; idcoordy++) {
                for (idcoordx = 0; idcoordx < nx; idcoordx++) {

                        x = (idcoordx + WarpX[idx] - offx) / mx;
                        y = (idcoordy + WarpY[idx] - offy) / my;

                        Bout[idx++] = InterpolatedGradYMirr(Bcoeff,
                                        Width, Height, x, y, SplineDegree)/my;
                }
        }
} /* end BatchInterpolatedGradY2DMirrWarp */

extern void     BatchInterpolatedValue2DZero (
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
{ /* begin BatchInterpolatedValue2DZero */

	double	x, y;
	long	idx;
	long	idcoordx, idcoordy;

	idx = 0;	
	for (idcoordy = 0; idcoordy < ny; idcoordy++) {
		for (idcoordx = 0; idcoordx < nx; idcoordx++) {
			
			x = (idcoordx - offx) / mx;
			y = (idcoordy - offy) / my;

			Bout[idx++] = InterpolatedValueZero(Bcoeff, 
					Width, Height, x, y, SplineDegree);
		}
	}
} /* end BatchInterpolatedValue2DZero */


extern void     BatchInterpolatedValue2DZeroWarp (
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
{ /* begin BatchInterpolatedValue2DZeroWarp */

        double  x, y;
        long    idx;
        long    idcoordx, idcoordy;

        idx = 0;
        for (idcoordy = 0; idcoordy < ny; idcoordy++) {
                for (idcoordx = 0; idcoordx < nx; idcoordx++) {

                        x = (idcoordx + WarpX[idx] - offx) / mx;
                        y = (idcoordy + WarpY[idx] - offy) / my;

                        Bout[idx++] = InterpolatedValueZero(Bcoeff,
                                        Width, Height, x, y, SplineDegree);
                }
        }
} /* end BatchInterpolatedValue2DZeroWarp */


extern void     BatchInterpolatedGradX2DZero (
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
{ /* begin BatchInterpolatedGradX2DZero */

	double	x, y;
	long	idx;
	long	idcoordx, idcoordy;

	idx = 0;	
	for (idcoordy = 0; idcoordy < ny; idcoordy++) {
		for (idcoordx = 0; idcoordx < nx; idcoordx++) {
			
			x = (idcoordx - offx) / mx;
			y = (idcoordy - offy) / my;

			Bout[idx++] = InterpolatedGradXZero(Bcoeff, 
					Width, Height, x, y, SplineDegree)/mx;
		}
	}
} /* end BatchInterpolatedGradX2DZero */

extern void     BatchInterpolatedGradX2DZeroWarp (
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
{ /* begin BatchInterpolatedGradX2DZeroWarp */

        double  x, y;
        long    idx;
        long    idcoordx, idcoordy;

        idx = 0;
        for (idcoordy = 0; idcoordy < ny; idcoordy++) {
                for (idcoordx = 0; idcoordx < nx; idcoordx++) {

                        x = (idcoordx + WarpX[idx] - offx) / mx;
                        y = (idcoordy + WarpY[idx] - offy) / my;

                        Bout[idx++] = InterpolatedGradXZero(Bcoeff,
                                        Width, Height, x, y, SplineDegree)/mx;
                }
        }
} /* end BatchInterpolatedGradX2DZeroWarp */


extern void     BatchInterpolatedGradY2DZero (
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
{ /* begin BatchInterpolatedGradY2DZero */

	double	x, y;
	long	idx;
	long	idcoordx, idcoordy;

	idx = 0;	
	for (idcoordy = 0; idcoordy < ny; idcoordy++) {
		for (idcoordx = 0; idcoordx < nx; idcoordx++) {
			
			x = (idcoordx - offx) / mx;
			y = (idcoordy - offy) / my;

			Bout[idx++] = InterpolatedGradYZero(Bcoeff, 
					Width, Height, x, y, SplineDegree)/my;
		}
	}
} /* end BatchInterpolatedGradY2DZero */

extern void     BatchInterpolatedGradY2DZeroWarp (
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
{ /* begin BatchInterpolatedGradY2DZeroWarp */

        double  x, y;
        long    idx;
        long    idcoordx, idcoordy;

        idx = 0;
        for (idcoordy = 0; idcoordy < ny; idcoordy++) {
                for (idcoordx = 0; idcoordx < nx; idcoordx++) {

                        x = (idcoordx + WarpX[idx] - offx) / mx;
                        y = (idcoordy + WarpY[idx] - offy) / my;

                        Bout[idx++] = InterpolatedGradYZero(Bcoeff,
                                        Width, Height, x, y, SplineDegree)/my;
                }
        }
} /* end BatchInterpolatedGradY2DZeroWarp */
