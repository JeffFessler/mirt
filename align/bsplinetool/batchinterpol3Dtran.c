/*******************************************************************************

        Batch transpose interpolation for all locations of points
        with mirror and zero end condition

        Created by Se Young Chun, Jan 19, 2007, the University of Michigan

*******************************************************************************/

#include	<math.h>
#include	<stddef.h>
#include	<stdio.h>
#include	<stdlib.h>

#include        "batchinterpol3Dtran.h"
#include	"interpol3Dtran.h"
 
extern void     BatchInterpolatedValue3DTranMirr (
                        float   *Bout,    
                        float   *Bcoeff, 
                        long    Width,  
                        long    Height,
                        long    Slice,
                        double  nx,          
                        double  offx,       
                        double  mx,        
                        double  ny,       
                        double  offy,    
                        double  my,     
                        double  nz,    
                        double  offz, 
                        double  mz,  
                        long    SplineDegree  
                )
{ /* begin BatchInterpolatedValue3DTranMirr */
	double	x, y, z;
	long	idx;
	long	idcoordx, idcoordy, idcoordz;

	for (idx = 0; idx < nx*ny*nz; idx++) Bout[idx] = 0;

	idx = 0;	
	for (idcoordz = 0; idcoordz < Slice; idcoordz++) {
		for (idcoordy = 0; idcoordy < Height; idcoordy++) {
			for (idcoordx = 0; idcoordx < Width; idcoordx++) {
				
				x = (idcoordx - offx) / mx;
				y = (idcoordy - offy) / my;
				z = (idcoordz - offz) / mz;

				InterpolatedValue3DTranMirr(Bout, nx, ny, nz,
 					Bcoeff[idx++], x, y, z, SplineDegree);
			}
		}
	}
} /* end BatchInterpolatedValue3DTranMirr */

extern void     BatchInterpolatedValue3DTranMirrWarp (
                        float   *Bout,       
                        float   *Bcoeff,    
                        long    Width,     
                        long    Height,   
                        long    Slice,   
                        double  nx,     
                        double  offx,  
                        double  mx,   
                        float   *WarpX,        
                        double  ny,           
                        double  offy,        
                        double  my,         
                        float   *WarpY,    
                        double  nz,       
                        double  offz,    
                        double  mz,     
                        float   *WarpZ,
                        long    SplineDegree 
                )
{ /* begin BatchInterpolatedValue3DTranMirrWarp */
        double  x, y, z;
        long    idx;
        long    idcoordx, idcoordy, idcoordz;

        for (idx = 0; idx < nx*ny*nz; idx++) Bout[idx] = 0;

        idx = 0;
        for (idcoordz = 0; idcoordz < Slice; idcoordz++) {
                for (idcoordy = 0; idcoordy < Height; idcoordy++) {
                        for (idcoordx = 0; idcoordx < Width; idcoordx++) {

                                x = (idcoordx + WarpX[idx] - offx) / mx;
                                y = (idcoordy + WarpY[idx] - offy) / my;
                                z = (idcoordz + WarpZ[idx] - offz) / mz;

                                InterpolatedValue3DTranMirr(Bout, nx, ny, nz,
 					Bcoeff[idx++], x, y, z, SplineDegree);
                        }
                }
        }
} /* end BatchInterpolatedValue3DTranMirrWarp */

extern void     BatchInterpolatedGradX3DTranMirr (
                        float   *Bout,       
                        float   *Bcoeff,    
                        long    Width,     
                        long    Height,   
                        long    Slice,   
                        double  nx,     
                        double  offx,  
                        double  mx,   
                        double  ny,  
                        double  offy,  
                        double  my,   
                        double  nz, 
                        double  offz,
                        double  mz, 
                        long    SplineDegree  
                )
{ /* begin BatchInterpolatedGradX3DTranMirr */

	double	x, y, z;
	long	idx;
	long	idcoordx, idcoordy, idcoordz;

	for (idx = 0; idx < nx*ny*nz; idx++) Bout[idx] = 0;

	idx = 0;	
	for (idcoordz = 0; idcoordz < Slice; idcoordz++) {
		for (idcoordy = 0; idcoordy < Height; idcoordy++) {
			for (idcoordx = 0; idcoordx < Width; idcoordx++) {
				
				x = (idcoordx - offx) / mx;
				y = (idcoordy - offy) / my;
				z = (idcoordz - offz) / mz;

				InterpolatedGradX3DTranMirr(Bout, nx, ny, nz,
 					Bcoeff[idx++], x, y, z, SplineDegree);
			}
		}
	}

	for (idx = 0; idx < nx*ny*nz; idx++) Bout[idx] = Bout[idx] / mx;

} /* end BatchInterpolatedGradX3DTranMirr */

extern void     BatchInterpolatedGradX3DTranMirrWarp (
                        float   *Bout,       
                        float   *Bcoeff,    
                        long    Width,     
                        long    Height,   
                        long    Slice,   
                        double  nx,     
                        double  offx,  
                        double  mx,   
                        float   *WarpX,   
                        double  ny,      
                        double  offy,   
                        double  my,    
                        float   *WarpY,  
                        double  nz,     
                        double  offz,  
                        double  mz,   
                        float   *WarpZ,   
                        long    SplineDegree 
                )
{ /* begin BatchInterpolatedGradX3DTranMirrWarp */

        double  x, y, z;
        long    idx;
        long    idcoordx, idcoordy, idcoordz;

        for (idx = 0; idx < nx*ny*nz; idx++) Bout[idx] = 0;

        idx = 0;
        for (idcoordz = 0; idcoordz < Slice; idcoordz++) {
                for (idcoordy = 0; idcoordy < Height; idcoordy++) {
                        for (idcoordx = 0; idcoordx < Width; idcoordx++) {

                                x = (idcoordx + WarpX[idx] - offx) / mx;
                                y = (idcoordy + WarpY[idx] - offy) / my;
                                z = (idcoordz + WarpZ[idx] - offz) / mz;

                                InterpolatedGradX3DTranMirr(Bout, nx, ny, nz,
 					Bcoeff[idx++], x, y, z, SplineDegree);
                        }
                }
        }

        for (idx = 0; idx < nx*ny*nz; idx++) Bout[idx] = Bout[idx] / mx;

} /* end BatchInterpolatedGradX3DTranMirrWarp */

extern void     BatchInterpolatedGradY3DTranMirr (
                        float   *Bout,      
                        float   *Bcoeff,   
                        long    Width,    
                        long    Height,  
                        long    Slice,  
                        double  nx,    
                        double  offx, 
                        double  mx,  
                        double  ny,         
                        double  offy,      
                        double  my,       
                        double  nz,      
                        double  offz,   
                        double  mz,    
                        long    SplineDegree 
                )
{ /* begin BatchInterpolatedGradY3DTranMirr */
	double	x, y, z;
	long	idx;
	long	idcoordx, idcoordy, idcoordz;

	for (idx = 0; idx < nx*ny*nz; idx++) Bout[idx] = 0;

	idx = 0;	
	for (idcoordz = 0; idcoordz < Slice; idcoordz++) {
		for (idcoordy = 0; idcoordy < Height; idcoordy++) {
			for (idcoordx = 0; idcoordx < Width; idcoordx++) {
				
				x = (idcoordx - offx) / mx;
				y = (idcoordy - offy) / my;
				z = (idcoordz - offz) / mz;

				InterpolatedGradY3DTranMirr(Bout, nx, ny, nz,
 					Bcoeff[idx++], x, y, z, SplineDegree);
			}
		}
	}

	for (idx = 0; idx < nx*ny*nz; idx++) Bout[idx] = Bout[idx] / my;

} /* end BatchInterpolatedGradY3DTranMirr */

extern void     BatchInterpolatedGradY3DTranMirrWarp (
                        float   *Bout,      
                        float   *Bcoeff,   
                        long    Width,    
                        long    Height,  
                        long    Slice,  
                        double  nx,    
                        double  offx, 
                        double  mx,  
                        float   *WarpX,   
                        double  ny,      
                        double  offy,   
                        double  my,    
                        float   *WarpY,    
                        double  nz,       
                        double  offz,    
                        double  mz,     
                        float   *WarpZ,   
                        long    SplineDegree 
                )
{ /* begin BatchInterpolatedGradY3DTranMirrWarp */
        double  x, y, z;
        long    idx;
        long    idcoordx, idcoordy, idcoordz;

        for (idx = 0; idx < nx*ny*nz; idx++) Bout[idx] = 0;

        idx = 0;
        for (idcoordz = 0; idcoordz < Slice; idcoordz++) {
                for (idcoordy = 0; idcoordy < Height; idcoordy++) {
                        for (idcoordx = 0; idcoordx < Width; idcoordx++) {

                                x = (idcoordx + WarpX[idx] - offx) / mx;
                                y = (idcoordy + WarpY[idx] - offy) / my;
                                z = (idcoordz + WarpZ[idx] - offz) / mz;

                                InterpolatedGradY3DTranMirr(Bout, nx, ny, nz,
 					Bcoeff[idx++], x, y, z, SplineDegree);
                        }
                }
        }

        for (idx = 0; idx < nx*ny*nz; idx++) Bout[idx] = Bout[idx] / my;

} /* end BatchInterpolatedGradY3DTranMirrWarp */


extern void     BatchInterpolatedGradZ3DTranMirr ( 
                        float   *Bout,      
                        float   *Bcoeff,   
                        long    Width,    
                        long    Height,  
                        long    Slice,  
                        double  nx,    
                        double  offx, 
                        double  mx,  
                        double  ny,       
                        double  offy,    
                        double  my,     
                        double  nz,    
                        double  offz, 
                        double  mz,  
                        long    SplineDegree 
                )
{ /* begin BatchInterpolatedGradZ3DTranMirr */
	double	x, y, z;
	long	idx;
	long	idcoordx, idcoordy, idcoordz;

	for (idx = 0; idx < nx*ny*nz; idx++) Bout[idx] = 0;

	idx = 0;	
	for (idcoordz = 0; idcoordz < Slice; idcoordz++) {
		for (idcoordy = 0; idcoordy < Height; idcoordy++) {
			for (idcoordx = 0; idcoordx < Width; idcoordx++) {
				
				x = (idcoordx - offx) / mx;
				y = (idcoordy - offy) / my;
				z = (idcoordz - offz) / mz;

				InterpolatedGradZ3DTranMirr(Bout, nx, ny, nz,
		 			Bcoeff[idx++], x, y, z, SplineDegree);
			}
		}
	}

	for (idx = 0; idx < nx*ny*nz; idx++) Bout[idx] = Bout[idx] / mz;

} /* end BatchInterpolatedGradZ3DTranMirr */

extern void     BatchInterpolatedGradZ3DTranMirrWarp (
                        float   *Bout,      
                        float   *Bcoeff,   
                        long    Width,    
                        long    Height,  
                        long    Slice,  
                        double  nx,    
                        double  offx, 
                        double  mx,  
                        float   *WarpX,    
                        double  ny,       
                        double  offy,    
                        double  my,     
                        float   *WarpY,
                        double  nz,   
                        double  offz,
                        double  mz, 
                        float   *WarpZ,    
                        long    SplineDegree 
                )
{ /* begin BatchInterpolatedGradZ3DTranMirrWarp */
        double  x, y, z;
        long    idx;
        long    idcoordx, idcoordy, idcoordz;

        for (idx = 0; idx < nx*ny*nz; idx++) Bout[idx] = 0;

        idx = 0;
        for (idcoordz = 0; idcoordz < Slice; idcoordz++) {
                for (idcoordy = 0; idcoordy < Height; idcoordy++) {
                        for (idcoordx = 0; idcoordx < Width; idcoordx++) {

                                x = (idcoordx + WarpX[idx] - offx) / mx;
                                y = (idcoordy + WarpY[idx] - offy) / my;
                                z = (idcoordz + WarpZ[idx] - offz) / mz;

                                InterpolatedGradZ3DTranMirr(Bout, nx, ny, nz,
				 	Bcoeff[idx++], x, y, z, SplineDegree);
                        }
                }
        }

        for (idx = 0; idx < nx*ny*nz; idx++) Bout[idx] = Bout[idx] / mz;

} /* end BatchInterpolatedGradZ3DTranMirrWarp */

extern void     BatchInterpolatedValue3DTranZero (
                        float   *Bout,    
                        float   *Bcoeff, 
                        long    Width,  
                        long    Height,
                        long    Slice,
                        double  nx,          
                        double  offx,       
                        double  mx,        
                        double  ny,       
                        double  offy,    
                        double  my,     
                        double  nz,    
                        double  offz, 
                        double  mz,  
                        long    SplineDegree  
                )
{ /* begin BatchInterpolatedValue3DTranZero */
	double	x, y, z;
	long	idx;
	long	idcoordx, idcoordy, idcoordz;

	for (idx = 0; idx < nx*ny*nz; idx++) Bout[idx] = 0;

	idx = 0;	
	for (idcoordz = 0; idcoordz < Slice; idcoordz++) {
		for (idcoordy = 0; idcoordy < Height; idcoordy++) {
			for (idcoordx = 0; idcoordx < Width; idcoordx++) {
				
				x = (idcoordx - offx) / mx;
				y = (idcoordy - offy) / my;
				z = (idcoordz - offz) / mz;

				InterpolatedValue3DTranZero(Bout, nx, ny, nz,
 					Bcoeff[idx++], x, y, z, SplineDegree);
			}
		}
	}
} /* end BatchInterpolatedValue3DTranZero */

extern void     BatchInterpolatedValue3DTranZeroWarp (
                        float   *Bout,       
                        float   *Bcoeff,    
                        long    Width,     
                        long    Height,   
                        long    Slice,   
                        double  nx,     
                        double  offx,  
                        double  mx,   
                        float   *WarpX,        
                        double  ny,           
                        double  offy,        
                        double  my,         
                        float   *WarpY,    
                        double  nz,       
                        double  offz,    
                        double  mz,     
                        float   *WarpZ,
                        long    SplineDegree 
                )
{ /* begin BatchInterpolatedValue3DTranZeroWarp */
        double  x, y, z;
        long    idx;
        long    idcoordx, idcoordy, idcoordz;

        for (idx = 0; idx < nx*ny*nz; idx++) Bout[idx] = 0;

        idx = 0;
        for (idcoordz = 0; idcoordz < Slice; idcoordz++) {
                for (idcoordy = 0; idcoordy < Height; idcoordy++) {
                        for (idcoordx = 0; idcoordx < Width; idcoordx++) {

                                x = (idcoordx + WarpX[idx] - offx) / mx;
                                y = (idcoordy + WarpY[idx] - offy) / my;
                                z = (idcoordz + WarpZ[idx] - offz) / mz;

                                InterpolatedValue3DTranZero(Bout, nx, ny, nz,
 					Bcoeff[idx++], x, y, z, SplineDegree);
                        }
                }
        }
} /* end BatchInterpolatedValue3DTranZeroWarp */

extern void     BatchInterpolatedGradX3DTranZero (
                        float   *Bout,       
                        float   *Bcoeff,    
                        long    Width,     
                        long    Height,   
                        long    Slice,   
                        double  nx,     
                        double  offx,  
                        double  mx,   
                        double  ny,  
                        double  offy,  
                        double  my,   
                        double  nz, 
                        double  offz,
                        double  mz, 
                        long    SplineDegree  
                )
{ /* begin BatchInterpolatedGradX3DTranZero */

	double	x, y, z;
	long	idx;
	long	idcoordx, idcoordy, idcoordz;

	for (idx = 0; idx < nx*ny*nz; idx++) Bout[idx] = 0;

	idx = 0;	
	for (idcoordz = 0; idcoordz < Slice; idcoordz++) {
		for (idcoordy = 0; idcoordy < Height; idcoordy++) {
			for (idcoordx = 0; idcoordx < Width; idcoordx++) {
				
				x = (idcoordx - offx) / mx;
				y = (idcoordy - offy) / my;
				z = (idcoordz - offz) / mz;

				InterpolatedGradX3DTranZero(Bout, nx, ny, nz,
 					Bcoeff[idx++], x, y, z, SplineDegree);
			}
		}
	}

	for (idx = 0; idx < nx*ny*nz; idx++) Bout[idx] = Bout[idx] / mx;

} /* end BatchInterpolatedGradX3DTranZero */

extern void     BatchInterpolatedGradX3DTranZeroWarp (
                        float   *Bout,       
                        float   *Bcoeff,    
                        long    Width,     
                        long    Height,   
                        long    Slice,   
                        double  nx,     
                        double  offx,  
                        double  mx,   
                        float   *WarpX,   
                        double  ny,      
                        double  offy,   
                        double  my,    
                        float   *WarpY,  
                        double  nz,     
                        double  offz,  
                        double  mz,   
                        float   *WarpZ,   
                        long    SplineDegree 
                )
{ /* begin BatchInterpolatedGradX3DTranZeroWarp */

        double  x, y, z;
        long    idx;
        long    idcoordx, idcoordy, idcoordz;

        for (idx = 0; idx < nx*ny*nz; idx++) Bout[idx] = 0;

        idx = 0;
        for (idcoordz = 0; idcoordz < Slice; idcoordz++) {
                for (idcoordy = 0; idcoordy < Height; idcoordy++) {
                        for (idcoordx = 0; idcoordx < Width; idcoordx++) {

                                x = (idcoordx + WarpX[idx] - offx) / mx;
                                y = (idcoordy + WarpY[idx] - offy) / my;
                                z = (idcoordz + WarpZ[idx] - offz) / mz;

                                InterpolatedGradX3DTranZero(Bout, nx, ny, nz,
 					Bcoeff[idx++], x, y, z, SplineDegree);
                        }
                }
        }

        for (idx = 0; idx < nx*ny*nz; idx++) Bout[idx] = Bout[idx] / mx;

} /* end BatchInterpolatedGradX3DTranZeroWarp */

extern void     BatchInterpolatedGradY3DTranZero (
                        float   *Bout,      
                        float   *Bcoeff,   
                        long    Width,    
                        long    Height,  
                        long    Slice,  
                        double  nx,    
                        double  offx, 
                        double  mx,  
                        double  ny,         
                        double  offy,      
                        double  my,       
                        double  nz,      
                        double  offz,   
                        double  mz,    
                        long    SplineDegree 
                )
{ /* begin BatchInterpolatedGradY3DTranZero */
	double	x, y, z;
	long	idx;
	long	idcoordx, idcoordy, idcoordz;

	for (idx = 0; idx < nx*ny*nz; idx++) Bout[idx] = 0;

	idx = 0;	
	for (idcoordz = 0; idcoordz < Slice; idcoordz++) {
		for (idcoordy = 0; idcoordy < Height; idcoordy++) {
			for (idcoordx = 0; idcoordx < Width; idcoordx++) {
				
				x = (idcoordx - offx) / mx;
				y = (idcoordy - offy) / my;
				z = (idcoordz - offz) / mz;

				InterpolatedGradY3DTranZero(Bout, nx, ny, nz,
 					Bcoeff[idx++], x, y, z, SplineDegree);
			}
		}
	}

	for (idx = 0; idx < nx*ny*nz; idx++) Bout[idx] = Bout[idx] / my;

} /* end BatchInterpolatedGradY3DTranZero */

extern void     BatchInterpolatedGradY3DTranZeroWarp (
                        float   *Bout,      
                        float   *Bcoeff,   
                        long    Width,    
                        long    Height,  
                        long    Slice,  
                        double  nx,    
                        double  offx, 
                        double  mx,  
                        float   *WarpX,   
                        double  ny,      
                        double  offy,   
                        double  my,    
                        float   *WarpY,    
                        double  nz,       
                        double  offz,    
                        double  mz,     
                        float   *WarpZ,   
                        long    SplineDegree 
                )
{ /* begin BatchInterpolatedGradY3DTranZeroWarp */
        double  x, y, z;
        long    idx;
        long    idcoordx, idcoordy, idcoordz;

        for (idx = 0; idx < nx*ny*nz; idx++) Bout[idx] = 0;

        idx = 0;
        for (idcoordz = 0; idcoordz < Slice; idcoordz++) {
                for (idcoordy = 0; idcoordy < Height; idcoordy++) {
                        for (idcoordx = 0; idcoordx < Width; idcoordx++) {

                                x = (idcoordx + WarpX[idx] - offx) / mx;
                                y = (idcoordy + WarpY[idx] - offy) / my;
                                z = (idcoordz + WarpZ[idx] - offz) / mz;

                                InterpolatedGradY3DTranZero(Bout, nx, ny, nz,
 					Bcoeff[idx++], x, y, z, SplineDegree);
                        }
                }
        }

        for (idx = 0; idx < nx*ny*nz; idx++) Bout[idx] = Bout[idx] / my;

} /* end BatchInterpolatedGradY3DTranZeroWarp */


extern void     BatchInterpolatedGradZ3DTranZero ( 
                        float   *Bout,      
                        float   *Bcoeff,   
                        long    Width,    
                        long    Height,  
                        long    Slice,  
                        double  nx,    
                        double  offx, 
                        double  mx,  
                        double  ny,       
                        double  offy,    
                        double  my,     
                        double  nz,    
                        double  offz, 
                        double  mz,  
                        long    SplineDegree 
                )
{ /* begin BatchInterpolatedGradZ3DTranZero */
	double	x, y, z;
	long	idx;
	long	idcoordx, idcoordy, idcoordz;

	for (idx = 0; idx < nx*ny*nz; idx++) Bout[idx] = 0;

	idx = 0;	
	for (idcoordz = 0; idcoordz < Slice; idcoordz++) {
		for (idcoordy = 0; idcoordy < Height; idcoordy++) {
			for (idcoordx = 0; idcoordx < Width; idcoordx++) {
				
				x = (idcoordx - offx) / mx;
				y = (idcoordy - offy) / my;
				z = (idcoordz - offz) / mz;

				InterpolatedGradZ3DTranZero(Bout, nx, ny, nz,
		 			Bcoeff[idx++], x, y, z, SplineDegree);
			}
		}
	}

	for (idx = 0; idx < nx*ny*nz; idx++) Bout[idx] = Bout[idx] / mz;

} /* end BatchInterpolatedGradZ3DTranZero */

extern void     BatchInterpolatedGradZ3DTranZeroWarp (
                        float   *Bout,      
                        float   *Bcoeff,   
                        long    Width,    
                        long    Height,  
                        long    Slice,  
                        double  nx,    
                        double  offx, 
                        double  mx,  
                        float   *WarpX,    
                        double  ny,       
                        double  offy,    
                        double  my,     
                        float   *WarpY,
                        double  nz,   
                        double  offz,
                        double  mz, 
                        float   *WarpZ,    
                        long    SplineDegree 
                )
{ /* begin BatchInterpolatedGradZ3DTranZeroWarp */
        double  x, y, z;
        long    idx;
        long    idcoordx, idcoordy, idcoordz;

        for (idx = 0; idx < nx*ny*nz; idx++) Bout[idx] = 0;

        idx = 0;
        for (idcoordz = 0; idcoordz < Slice; idcoordz++) {
                for (idcoordy = 0; idcoordy < Height; idcoordy++) {
                        for (idcoordx = 0; idcoordx < Width; idcoordx++) {

                                x = (idcoordx + WarpX[idx] - offx) / mx;
                                y = (idcoordy + WarpY[idx] - offy) / my;
                                z = (idcoordz + WarpZ[idx] - offz) / mz;

                                InterpolatedGradZ3DTranZero(Bout, nx, ny, nz,
				 	Bcoeff[idx++], x, y, z, SplineDegree);
                        }
                }
        }

        for (idx = 0; idx < nx*ny*nz; idx++) Bout[idx] = Bout[idx] / mz;

} /* end BatchInterpolatedGradZ3DTranZeroWarp */
