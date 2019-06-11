/*******************************************************************************

        Batch transpose interpolation for all locations of points
	with mirror and zero end condition

        Created by Se Young Chun, Jan 19, 2007, the University of Michigan

*******************************************************************************/

extern void 	BatchInterpolatedValue2DTranMirr (
  			float   *Bout,          /* output interpolated */
                	float   *Bcoeff,        /* input B-spline coeff array */
                	long    Width,          /* width of coeff array */
                	long    Height,         /* height of coeff array */
                	double  nx,             /* width of Bout */
                	double  offx,           /* offset of x direction */
                	double  mx,             /* magnif factor for x */
                	double  ny,             /* height of Bout */
                	double  offy,           /* offset of y direction */
                	double  my,             /* magnif factor for y */
        		long    SplineDegree    /* degree of spline model */
		);

extern void     BatchInterpolatedValue2DTranMirrWarp (
  			float   *Bout,          /* output interpolated */
                        float   *Bcoeff,        /* input B-spline coeff array */
                        long    Width,          /* width of coeff array */
                        long    Height,         /* height of coeff array */
                        double  nx,             /* width of Bout */
                        double  offx,           /* offset of x direction */
                        double  mx,             /* magnif factor for x */
                        float   *WarpX,         /* deformed value in x direct */
                        double  ny,             /* height of Bout */
                        double  offy,           /* offset of y direction */
                        double  my,             /* magnif factor for y */
                        float   *WarpY,         /* deformed value in y direct */
                        long    SplineDegree    /* degree of spline model */
                );

extern void 	BatchInterpolatedGradX2DTranMirr (
  			float   *Bout,          /* output interpolated */
                	float   *Bcoeff,        /* input B-spline coeff array */
                	long    Width,          /* width of coeff array */
                	long    Height,         /* height of coeff array */
                	double  nx,             /* width of Bout */
                	double  offx,           /* offset of x direction */
                	double  mx,             /* magnif factor for x */
                	double  ny,             /* height of Bout */
                	double  offy,           /* offset of y direction */
                	double  my,             /* magnif factor for y */
                	long    SplineDegree    /* degree of spline model */
		);

extern void     BatchInterpolatedGradX2DTranMirrWarp (
  			float   *Bout,          /* output interpolated */
                        float   *Bcoeff,        /* input B-spline coeff array */
                        long    Width,          /* width of coeff array */
                        long    Height,         /* height of coeff array */
                        double  nx,             /* width of Bout */
                        double  offx,           /* offset of x direction */
                        double  mx,             /* magnif factor for x */
                        float   *WarpX,         /* deformed value in x direct */
                        double  ny,             /* height of Bout */
                        double  offy,           /* offset of y direction */
                        double  my,             /* magnif factor for y */
                        float   *WarpY,         /* deformed value in y direct */
                        long    SplineDegree    /* degree of spline model */
                );

extern void 	BatchInterpolatedGradY2DTranMirr (
  			float   *Bout,          /* output interpolated */
                	float   *Bcoeff,        /* input B-spline coeff array */
                	long    Width,          /* width of coeff array */
                	long    Height,         /* height of coeff array */
                	double  nx,             /* width of Bout */
                	double  offx,           /* offset of x direction */
                	double  mx,             /* magnif factor for x */
                	double  ny,             /* height of Bout */
                	double  offy,           /* offset of y direction */
                	double  my,             /* magnif factor for y */
                	long    SplineDegree    /* degree of spline model */
		);

extern void     BatchInterpolatedGradY2DTranMirrWarp (
  			float   *Bout,          /* output interpolated */
                        float   *Bcoeff,        /* input B-spline coeff array */
                        long    Width,          /* width of coeff array */
                        long    Height,         /* height of coeff array */
                        double  nx,             /* width of Bout */
                        double  offx,           /* offset of x direction */
                        double  mx,             /* magnif factor for x */
                        float   *WarpX,         /* deformed value in x direct */
                        double  ny,             /* height of Bout */
                        double  offy,           /* offset of y direction */
                        double  my,             /* magnif factor for y */
                        float   *WarpY,         /* deformed value in y direct */
                        long    SplineDegree    /* degree of spline model */
                );

extern void 	BatchInterpolatedValue2DTranZero (
  			float   *Bout,          /* output interpolated */
                	float   *Bcoeff,        /* input B-spline coeff array */
                	long    Width,          /* width of coeff array */
                	long    Height,         /* height of coeff array */
                	double  nx,             /* width of Bout */
                	double  offx,           /* offset of x direction */
                	double  mx,             /* magnif factor for x */
                	double  ny,             /* height of Bout */
                	double  offy,           /* offset of y direction */
                	double  my,             /* magnif factor for y */
        		long    SplineDegree    /* degree of spline model */
		);

extern void     BatchInterpolatedValue2DTranZeroWarp (
  			float   *Bout,          /* output interpolated */
                        float   *Bcoeff,        /* input B-spline coeff array */
                        long    Width,          /* width of coeff array */
                        long    Height,         /* height of coeff array */
                        double  nx,             /* width of Bout */
                        double  offx,           /* offset of x direction */
                        double  mx,             /* magnif factor for x */
                        float   *WarpX,         /* deformed value in x direct */
                        double  ny,             /* height of Bout */
                        double  offy,           /* offset of y direction */
                        double  my,             /* magnif factor for y */
                        float   *WarpY,         /* deformed value in y direct */
                        long    SplineDegree    /* degree of spline model */
                );

extern void 	BatchInterpolatedGradX2DTranZero (
  			float   *Bout,          /* output interpolated */
                	float   *Bcoeff,        /* input B-spline coeff array */
                	long    Width,          /* width of coeff array */
                	long    Height,         /* height of coeff array */
                	double  nx,             /* width of Bout */
                	double  offx,           /* offset of x direction */
                	double  mx,             /* magnif factor for x */
                	double  ny,             /* height of Bout */
                	double  offy,           /* offset of y direction */
                	double  my,             /* magnif factor for y */
                	long    SplineDegree    /* degree of spline model */
		);

extern void     BatchInterpolatedGradX2DTranZeroWarp (
  			float   *Bout,          /* output interpolated */
                        float   *Bcoeff,        /* input B-spline coeff array */
                        long    Width,          /* width of coeff array */
                        long    Height,         /* height of coeff array */
                        double  nx,             /* width of Bout */
                        double  offx,           /* offset of x direction */
                        double  mx,             /* magnif factor for x */
                        float   *WarpX,         /* deformed value in x direct */
                        double  ny,             /* height of Bout */
                        double  offy,           /* offset of y direction */
                        double  my,             /* magnif factor for y */
                        float   *WarpY,         /* deformed value in y direct */
                        long    SplineDegree    /* degree of spline model */
                );

extern void 	BatchInterpolatedGradY2DTranZero (
  			float   *Bout,          /* output interpolated */
                	float   *Bcoeff,        /* input B-spline coeff array */
                	long    Width,          /* width of coeff array */
                	long    Height,         /* height of coeff array */
                	double  nx,             /* width of Bout */
                	double  offx,           /* offset of x direction */
                	double  mx,             /* magnif factor for x */
                	double  ny,             /* height of Bout */
                	double  offy,           /* offset of y direction */
                	double  my,             /* magnif factor for y */
                	long    SplineDegree    /* degree of spline model */
		);

extern void     BatchInterpolatedGradY2DTranZeroWarp (
  			float   *Bout,          /* output interpolated */
                        float   *Bcoeff,        /* input B-spline coeff array */
                        long    Width,          /* width of coeff array */
                        long    Height,         /* height of coeff array */
                        double  nx,             /* width of Bout */
                        double  offx,           /* offset of x direction */
                        double  mx,             /* magnif factor for x */
                        float   *WarpX,         /* deformed value in x direct */
                        double  ny,             /* height of Bout */
                        double  offy,           /* offset of y direction */
                        double  my,             /* magnif factor for y */
                        float   *WarpY,         /* deformed value in y direct */
                        long    SplineDegree    /* degree of spline model */
                );
