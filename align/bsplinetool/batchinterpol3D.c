/*******************************************************************************

	Batch interpolation for each location of points
        using point interpolation function, 3D case

        Created by Se Young Chun, Nov 28, 2006, the University of Michigan

*******************************************************************************/

#include	<math.h>
#include	<stddef.h>
#include	<stdio.h>
#include	<stdlib.h>

#include        "batchinterpol3D.h"
#include	"interpol3D.h"
 
extern void     BatchInterpolatedValue3DMirr (
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
{ /* begin BatchInterpolatedValue3DMirr */
	double	x, y, z;
	long	idx;
	long	idcoordx, idcoordy, idcoordz;

	idx = 0;	
	for (idcoordz = 0; idcoordz < nz; idcoordz++) {
		for (idcoordy = 0; idcoordy < ny; idcoordy++) {
			for (idcoordx = 0; idcoordx < nx; idcoordx++) {
				
				x = (idcoordx - offx) / mx;
				y = (idcoordy - offy) / my;
				z = (idcoordz - offz) / mz;

				Bout[idx++] = InterpolatedValue3DMirr(
						Bcoeff, Width, Height, Slice, 
						x, y, z, SplineDegree);
			}
		}
	}
} /* end BatchInterpolatedValue3DMirr */


extern void     BatchInterpolatedValue3DMirrWarp (
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
{ /* begin BatchInterpolatedValue3DMirrWarp */
        double  x, y, z;
        long    idx;
        long    idcoordx, idcoordy, idcoordz;

        idx = 0;
        for (idcoordz = 0; idcoordz < nz; idcoordz++) {
                for (idcoordy = 0; idcoordy < ny; idcoordy++) {
                        for (idcoordx = 0; idcoordx < nx; idcoordx++) {

                                x = (idcoordx + WarpX[idx] - offx) / mx;
                                y = (idcoordy + WarpY[idx] - offy) / my;
                                z = (idcoordz + WarpZ[idx] - offz) / mz;

                                Bout[idx++] = InterpolatedValue3DMirr(
						Bcoeff, Width, Height, Slice, 
						x, y, z, SplineDegree);
                        }
                }
        }
} /* end BatchInterpolatedValue3DMirrWarp */

extern void     BatchInterpolatedGradX3DMirr (
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
{ /* begin BatchInterpolatedGradX3DMirr */

	double	x, y, z;
	long	idx;
	long	idcoordx, idcoordy, idcoordz;

	idx = 0;	
	for (idcoordz = 0; idcoordz < nz; idcoordz++) {
		for (idcoordy = 0; idcoordy < ny; idcoordy++) {
			for (idcoordx = 0; idcoordx < nx; idcoordx++) {
				
				x = (idcoordx - offx) / mx;
				y = (idcoordy - offy) / my;
				z = (idcoordz - offz) / mz;

				Bout[idx++] = InterpolatedGradX3DMirr(
						Bcoeff, Width, Height, Slice, 
						x, y, z, SplineDegree) / mx;
			}
		}
	}

} /* end BatchInterpolatedGradX3DMirr */

extern void     BatchInterpolatedGradX3DMirrWarp (
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
{ /* begin BatchInterpolatedGradX3DMirrWarp */

        double  x, y, z;
        long    idx;
        long    idcoordx, idcoordy, idcoordz;

        idx = 0;
        for (idcoordz = 0; idcoordz < nz; idcoordz++) {
                for (idcoordy = 0; idcoordy < ny; idcoordy++) {
                        for (idcoordx = 0; idcoordx < nx; idcoordx++) {

                                x = (idcoordx + WarpX[idx] - offx) / mx;
                                y = (idcoordy + WarpY[idx] - offy) / my;
                                z = (idcoordz + WarpZ[idx] - offz) / mz;

                                Bout[idx++] = InterpolatedGradX3DMirr(
						Bcoeff, Width, Height, Slice, 
						x, y, z, SplineDegree) / mx;
                        }
                }
        }

} /* end BatchInterpolatedGradX3DMirrWarp */

extern void     BatchInterpolatedGradY3DMirr (
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
{ /* begin BatchInterpolatedGradY3DMirr */
	double	x, y, z;
	long	idx;
	long	idcoordx, idcoordy, idcoordz;

	idx = 0;	
	for (idcoordz = 0; idcoordz < nz; idcoordz++) {
		for (idcoordy = 0; idcoordy < ny; idcoordy++) {
			for (idcoordx = 0; idcoordx < nx; idcoordx++) {
				
				x = (idcoordx - offx) / mx;
				y = (idcoordy - offy) / my;
				z = (idcoordz - offz) / mz;

				Bout[idx++] = InterpolatedGradY3DMirr(
						Bcoeff, Width, Height, Slice, 
						x, y, z, SplineDegree) / my;
			}
		}
	}
} /* end BatchInterpolatedGradY3DMirr */


extern void     BatchInterpolatedGradY3DMirrWarp (
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
{ /* begin BatchInterpolatedGradY3DMirrWarp */
        double  x, y, z;
        long    idx;
        long    idcoordx, idcoordy, idcoordz;

        idx = 0;
        for (idcoordz = 0; idcoordz < nz; idcoordz++) {
                for (idcoordy = 0; idcoordy < ny; idcoordy++) {
                        for (idcoordx = 0; idcoordx < nx; idcoordx++) {

                                x = (idcoordx + WarpX[idx] - offx) / mx;
                                y = (idcoordy + WarpY[idx] - offy) / my;
                                z = (idcoordz + WarpZ[idx] - offz) / mz;

                                Bout[idx++] = InterpolatedGradY3DMirr(
						Bcoeff, Width, Height, Slice, 
						x, y, z, SplineDegree) / my;
                        }
                }
        }
} /* end BatchInterpolatedGradY3DMirrWarp */

extern void     BatchInterpolatedGradZ3DMirr ( 
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
{ /* begin BatchInterpolatedGradZ3DMirr */
	double	x, y, z;
	long	idx;
	long	idcoordx, idcoordy, idcoordz;

	idx = 0;	
	for (idcoordz = 0; idcoordz < nz; idcoordz++) {
		for (idcoordy = 0; idcoordy < ny; idcoordy++) {
			for (idcoordx = 0; idcoordx < nx; idcoordx++) {
				
				x = (idcoordx - offx) / mx;
				y = (idcoordy - offy) / my;
				z = (idcoordz - offz) / mz;

				Bout[idx++] = InterpolatedGradZ3DMirr(
						Bcoeff, Width, Height, Slice, 
						x, y, z, SplineDegree) / mz;
			}
		}
	}
} /* end BatchInterpolatedGradZ3DMirr */

extern void     BatchInterpolatedGradZ3DMirrWarp (
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
{ /* begin BatchInterpolatedGradZ3DMirrWarp */
        double  x, y, z;
        long    idx;
        long    idcoordx, idcoordy, idcoordz;

        idx = 0;
        for (idcoordz = 0; idcoordz < nz; idcoordz++) {
                for (idcoordy = 0; idcoordy < ny; idcoordy++) {
                        for (idcoordx = 0; idcoordx < nx; idcoordx++) {

                                x = (idcoordx + WarpX[idx] - offx) / mx;
                                y = (idcoordy + WarpY[idx] - offy) / my;
                                z = (idcoordz + WarpZ[idx] - offz) / mz;

                                Bout[idx++] = InterpolatedGradZ3DMirr(
						Bcoeff, Width, Height, Slice, 
						x, y, z, SplineDegree) / mz;
                        }
                }
        }
} /* end BatchInterpolatedGradZ3DMirrWarp */

extern void     BatchInterpolatedAll3DMirr (
                        float   *Bo, 
                        float   *Bogx, 
                        float   *Bogy, 
                        float   *Bogz, 
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
{ /* begin BatchInterpolatedAll3DMirr */
	double	x, y, z, gradx, grady, gradz;
	long	idx;
	long	idcoordx, idcoordy, idcoordz;

	idx = 0;	
	for (idcoordz = 0; idcoordz < nz; idcoordz++) {
		for (idcoordy = 0; idcoordy < ny; idcoordy++) {
			for (idcoordx = 0; idcoordx < nx; idcoordx++) {
				
				x = (idcoordx - offx) / mx;
				y = (idcoordy - offy) / my;
				z = (idcoordz - offz) / mz;

				Bo[idx] = InterpolatedAll3DMirr(
						&gradx, &grady, &gradz,
						Bcoeff, Width, Height, Slice, 
						x, y, z, SplineDegree);
				Bogx[idx] = gradx / mx;	
				Bogy[idx] = grady / my;	
				Bogz[idx++] = gradz / mz;	
			}
		}
	}
} /* end BatchInterpolatedAll3DMirr */

extern void     BatchInterpolatedAll3DMirrWarp (
                        float   *Bo, 
                        float   *Bogx, 
                        float   *Bogy, 
                        float   *Bogz, 
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
{ /* begin BatchInterpolatedAll3DMirrWarp */
        double  x, y, z, gradx, grady, gradz;
        long    idx;
        long    idcoordx, idcoordy, idcoordz;

        idx = 0;
        for (idcoordz = 0; idcoordz < nz; idcoordz++) {
                for (idcoordy = 0; idcoordy < ny; idcoordy++) {
                        for (idcoordx = 0; idcoordx < nx; idcoordx++) {

                                x = (idcoordx + WarpX[idx] - offx) / mx;
                                y = (idcoordy + WarpY[idx] - offy) / my;
                                z = (idcoordz + WarpZ[idx] - offz) / mz;

                                Bo[idx] = InterpolatedAll3DMirr(
						&gradx, &grady, &gradz,
						Bcoeff, Width, Height, Slice, 
						x, y, z, SplineDegree);
				Bogx[idx] = gradx / mx;	
				Bogy[idx] = grady / my;	
				Bogz[idx++] = gradz / mz;	
                        }
                }
        }
} /* end BatchInterpolatedAll3DMirrWarp */

extern void     BatchInterpolatedThree3DMirr (
                        float   *Bo1, 
                        float   *Bo1gx, 
                        float   *Bo1gy, 
                        float   *Bo1gz, 
                        float   *Bo2, 
                        float   *Bo2gx, 
                        float   *Bo2gy, 
                        float   *Bo2gz, 
                        float   *Bo3, 
                        float   *Bo3gx, 
                        float   *Bo3gy, 
                        float   *Bo3gz, 
                        float   *Bcoeff1, 
                        float   *Bcoeff2, 
                        float   *Bcoeff3, 
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
{ /* begin BatchInterpolatedThree3DMirr */
	double	x, y, z;
	double  out1, gradx1, grady1, gradz1;
	double  out2, gradx2, grady2, gradz2;
	double  out3, gradx3, grady3, gradz3;
	long	idx;
	long	idcoordx, idcoordy, idcoordz;

	idx = 0;	
	for (idcoordz = 0; idcoordz < nz; idcoordz++) {
		for (idcoordy = 0; idcoordy < ny; idcoordy++) {
			for (idcoordx = 0; idcoordx < nx; idcoordx++) {
				
				x = (idcoordx - offx) / mx;
				y = (idcoordy - offy) / my;
				z = (idcoordz - offz) / mz;

				InterpolatedThree3DMirr(
					&out1, &gradx1, &grady1, &gradz1,
					&out2, &gradx2, &grady2, &gradz2,
					&out3, &gradx3, &grady3, &gradz3,
					Bcoeff1, Bcoeff2, Bcoeff3, 
					Width, Height, Slice, 
					x, y, z, SplineDegree);
				Bo1[idx] = out1;	
				Bo1gx[idx] = gradx1 / mx;	
				Bo1gy[idx] = grady1 / my;	
				Bo1gz[idx] = gradz1 / mz;	
				Bo2[idx] = out2;	
				Bo2gx[idx] = gradx2 / mx;	
				Bo2gy[idx] = grady2 / my;	
				Bo2gz[idx] = gradz2 / mz;	
				Bo3[idx] = out3;	
				Bo3gx[idx] = gradx3 / mx;	
				Bo3gy[idx] = grady3 / my;	
				Bo3gz[idx++] = gradz3 / mz;	
			}
		}
	}
} /* end BatchInterpolatedThree3DMirr */

extern void     BatchInterpolatedValue3DZero (
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
{ /* begin BatchInterpolatedValue3DZero */
	double	x, y, z;
	long	idx;
	long	idcoordx, idcoordy, idcoordz;

	idx = 0;	
	for (idcoordz = 0; idcoordz < nz; idcoordz++) {
		for (idcoordy = 0; idcoordy < ny; idcoordy++) {
			for (idcoordx = 0; idcoordx < nx; idcoordx++) {
				
				x = (idcoordx - offx) / mx;
				y = (idcoordy - offy) / my;
				z = (idcoordz - offz) / mz;

				Bout[idx++] = InterpolatedValue3DZero(
						Bcoeff, Width, Height, Slice, 
						x, y, z, SplineDegree);
			}
		}
	}
} /* end BatchInterpolatedValue3DZero */


extern void     BatchInterpolatedValue3DZeroWarp (
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
{ /* begin BatchInterpolatedValue3DZeroWarp */
        double  x, y, z;
        long    idx;
        long    idcoordx, idcoordy, idcoordz;

        idx = 0;
        for (idcoordz = 0; idcoordz < nz; idcoordz++) {
                for (idcoordy = 0; idcoordy < ny; idcoordy++) {
                        for (idcoordx = 0; idcoordx < nx; idcoordx++) {

                                x = (idcoordx + WarpX[idx] - offx) / mx;
                                y = (idcoordy + WarpY[idx] - offy) / my;
                                z = (idcoordz + WarpZ[idx] - offz) / mz;

                                Bout[idx++] = InterpolatedValue3DZero(
						Bcoeff, Width, Height, Slice, 
						x, y, z, SplineDegree);
                        }
                }
        }
} /* end BatchInterpolatedValue3DZeroWarp */

extern void     BatchInterpolatedGradX3DZero (
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
{ /* begin BatchInterpolatedGradX3DZero */

	double	x, y, z;
	long	idx;
	long	idcoordx, idcoordy, idcoordz;

	idx = 0;	
	for (idcoordz = 0; idcoordz < nz; idcoordz++) {
		for (idcoordy = 0; idcoordy < ny; idcoordy++) {
			for (idcoordx = 0; idcoordx < nx; idcoordx++) {
				
				x = (idcoordx - offx) / mx;
				y = (idcoordy - offy) / my;
				z = (idcoordz - offz) / mz;

				Bout[idx++] = InterpolatedGradX3DZero(
						Bcoeff, Width, Height, Slice, 
						x, y, z, SplineDegree) / mx;
			}
		}
	}

} /* end BatchInterpolatedGradX3DZero */

extern void     BatchInterpolatedGradX3DZeroWarp (
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
{ /* begin BatchInterpolatedGradX3DZeroWarp */

        double  x, y, z;
        long    idx;
        long    idcoordx, idcoordy, idcoordz;

        idx = 0;
        for (idcoordz = 0; idcoordz < nz; idcoordz++) {
                for (idcoordy = 0; idcoordy < ny; idcoordy++) {
                        for (idcoordx = 0; idcoordx < nx; idcoordx++) {

                                x = (idcoordx + WarpX[idx] - offx) / mx;
                                y = (idcoordy + WarpY[idx] - offy) / my;
                                z = (idcoordz + WarpZ[idx] - offz) / mz;

                                Bout[idx++] = InterpolatedGradX3DZero(
						Bcoeff, Width, Height, Slice, 
						x, y, z, SplineDegree) / mx;
                        }
                }
        }

} /* end BatchInterpolatedGradX3DZeroWarp */

extern void     BatchInterpolatedGradY3DZero (
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
{ /* begin BatchInterpolatedGradY3DZero */
	double	x, y, z;
	long	idx;
	long	idcoordx, idcoordy, idcoordz;

	idx = 0;	
	for (idcoordz = 0; idcoordz < nz; idcoordz++) {
		for (idcoordy = 0; idcoordy < ny; idcoordy++) {
			for (idcoordx = 0; idcoordx < nx; idcoordx++) {
				
				x = (idcoordx - offx) / mx;
				y = (idcoordy - offy) / my;
				z = (idcoordz - offz) / mz;

				Bout[idx++] = InterpolatedGradY3DZero(
						Bcoeff, Width, Height, Slice, 
						x, y, z, SplineDegree) / my;
			}
		}
	}
} /* end BatchInterpolatedGradY3DZero */


extern void     BatchInterpolatedGradY3DZeroWarp (
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
{ /* begin BatchInterpolatedGradY3DZeroWarp */
        double  x, y, z;
        long    idx;
        long    idcoordx, idcoordy, idcoordz;

        idx = 0;
        for (idcoordz = 0; idcoordz < nz; idcoordz++) {
                for (idcoordy = 0; idcoordy < ny; idcoordy++) {
                        for (idcoordx = 0; idcoordx < nx; idcoordx++) {

                                x = (idcoordx + WarpX[idx] - offx) / mx;
                                y = (idcoordy + WarpY[idx] - offy) / my;
                                z = (idcoordz + WarpZ[idx] - offz) / mz;

                                Bout[idx++] = InterpolatedGradY3DZero(
						Bcoeff, Width, Height, Slice, 
						x, y, z, SplineDegree) / my;
                        }
                }
        }
} /* end BatchInterpolatedGradY3DZeroWarp */

extern void     BatchInterpolatedGradZ3DZero ( 
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
{ /* begin BatchInterpolatedGradZ3DZero */
	double	x, y, z;
	long	idx;
	long	idcoordx, idcoordy, idcoordz;

	idx = 0;	
	for (idcoordz = 0; idcoordz < nz; idcoordz++) {
		for (idcoordy = 0; idcoordy < ny; idcoordy++) {
			for (idcoordx = 0; idcoordx < nx; idcoordx++) {
				
				x = (idcoordx - offx) / mx;
				y = (idcoordy - offy) / my;
				z = (idcoordz - offz) / mz;

				Bout[idx++] = InterpolatedGradZ3DZero(
						Bcoeff, Width, Height, Slice, 
						x, y, z, SplineDegree) / mz;
			}
		}
	}
} /* end BatchInterpolatedGradZ3DZero */

extern void     BatchInterpolatedGradZ3DZeroWarp (
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
{ /* begin BatchInterpolatedGradZ3DZeroWarp */
        double  x, y, z;
        long    idx;
        long    idcoordx, idcoordy, idcoordz;

        idx = 0;
        for (idcoordz = 0; idcoordz < nz; idcoordz++) {
                for (idcoordy = 0; idcoordy < ny; idcoordy++) {
                        for (idcoordx = 0; idcoordx < nx; idcoordx++) {

                                x = (idcoordx + WarpX[idx] - offx) / mx;
                                y = (idcoordy + WarpY[idx] - offy) / my;
                                z = (idcoordz + WarpZ[idx] - offz) / mz;

                                Bout[idx++] = InterpolatedGradZ3DZero(
						Bcoeff, Width, Height, Slice, 
						x, y, z, SplineDegree) / mz;
                        }
                }
        }
} /* end BatchInterpolatedGradZ3DZeroWarp */

extern void     BatchInterpolatedAll3DZero (
                        float   *Bo, 
                        float   *Bogx, 
                        float   *Bogy, 
                        float   *Bogz, 
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
{ /* begin BatchInterpolatedAll3DZero */
	double	x, y, z, gradx, grady, gradz;
	long	idx;
	long	idcoordx, idcoordy, idcoordz;

	idx = 0;	
	for (idcoordz = 0; idcoordz < nz; idcoordz++) {
		for (idcoordy = 0; idcoordy < ny; idcoordy++) {
			for (idcoordx = 0; idcoordx < nx; idcoordx++) {
				
				x = (idcoordx - offx) / mx;
				y = (idcoordy - offy) / my;
				z = (idcoordz - offz) / mz;

				Bo[idx] = InterpolatedAll3DZero(
						&gradx, &grady, &gradz,
						Bcoeff, Width, Height, Slice, 
						x, y, z, SplineDegree);
				Bogx[idx] = gradx / mx;	
				Bogy[idx] = grady / my;	
				Bogz[idx++] = gradz / mz;	
			}
		}
	}
} /* end BatchInterpolatedAll3DZero */

extern void     BatchInterpolatedAll3DZeroWarp (
                        float   *Bo, 
                        float   *Bogx, 
                        float   *Bogy, 
                        float   *Bogz, 
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
{ /* begin BatchInterpolatedAll3DZeroWarp */
        double  x, y, z, gradx, grady, gradz;
        long    idx;
        long    idcoordx, idcoordy, idcoordz;

        idx = 0;
        for (idcoordz = 0; idcoordz < nz; idcoordz++) {
                for (idcoordy = 0; idcoordy < ny; idcoordy++) {
                        for (idcoordx = 0; idcoordx < nx; idcoordx++) {

                                x = (idcoordx + WarpX[idx] - offx) / mx;
                                y = (idcoordy + WarpY[idx] - offy) / my;
                                z = (idcoordz + WarpZ[idx] - offz) / mz;

                                Bo[idx] = InterpolatedAll3DZero(
						&gradx, &grady, &gradz,
						Bcoeff, Width, Height, Slice, 
						x, y, z, SplineDegree);
				Bogx[idx] = gradx / mx;	
				Bogy[idx] = grady /my;	
				Bogz[idx++] = gradz /mz;	
                        }
                }
        }
} /* end BatchInterpolatedAll3DZeroWarp */

extern void     BatchInterpolatedThree3DZero (
                        float   *Bo1, 
                        float   *Bo1gx, 
                        float   *Bo1gy, 
                        float   *Bo1gz, 
                        float   *Bo2, 
                        float   *Bo2gx, 
                        float   *Bo2gy, 
                        float   *Bo2gz, 
                        float   *Bo3, 
                        float   *Bo3gx, 
                        float   *Bo3gy, 
                        float   *Bo3gz, 
                        float   *Bcoeff1, 
                        float   *Bcoeff2, 
                        float   *Bcoeff3, 
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
                        float   *mask,   
                        long    SplineDegree  
                )
{ /* begin BatchInterpolatedThree3DZero */
	double	x, y, z;
	double  out1, gradx1, grady1, gradz1;
	double  out2, gradx2, grady2, gradz2;
	double  out3, gradx3, grady3, gradz3;
	long	idx;
	long	idcoordx, idcoordy, idcoordz;

	idx = 0;	
	for (idcoordz = 0; idcoordz < nz; idcoordz++) {
	   for (idcoordy = 0; idcoordy < ny; idcoordy++) {
	      for (idcoordx = 0; idcoordx < nx; idcoordx++) {
			
		 if (mask[idx] == 0) {			
			Bo1[idx] = 0;	
			Bo1gx[idx] = 0;	
			Bo1gy[idx] = 0;	
			Bo1gz[idx] = 0;	
			Bo2[idx] = 0;	
			Bo2gx[idx] = 0;	
			Bo2gy[idx] = 0;	
			Bo2gz[idx] = 0;	
			Bo3[idx] = 0;	
			Bo3gx[idx] = 0;	
			Bo3gy[idx] = 0;	
			Bo3gz[idx++] = 0;	
		 }
		 else {
			x = (idcoordx - offx) / mx;
			y = (idcoordy - offy) / my;
			z = (idcoordz - offz) / mz;

			InterpolatedThree3DZero(
				&out1, &gradx1, &grady1, &gradz1,
				&out2, &gradx2, &grady2, &gradz2,
				&out3, &gradx3, &grady3, &gradz3,
				Bcoeff1, Bcoeff2, Bcoeff3, 
				Width, Height, Slice, 
				x, y, z, SplineDegree);
			Bo1[idx] = out1;	
			Bo1gx[idx] = gradx1 / mx;	
			Bo1gy[idx] = grady1 / my;	
			Bo1gz[idx] = gradz1 / mz;	
			Bo2[idx] = out2;	
			Bo2gx[idx] = gradx2 / mx;	
			Bo2gy[idx] = grady2 / my;	
			Bo2gz[idx] = gradz2 / mz;	
			Bo3[idx] = out3;	
			Bo3gx[idx] = gradx3 / mx;	
			Bo3gy[idx] = grady3 / my;	
			Bo3gz[idx++] = gradz3 / mz;	
		 }
	      }
	   } 
	}
} /* end BatchInterpolatedThree3DZero */

