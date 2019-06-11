/*******************************************************************************

        Batch transpose interpolation for all locations of points
        with mirror and zero end condition

        Created by Se Young Chun, Jan 19, 2007, the University of Michigan

*******************************************************************************/

extern void BatchInterpolatedValue3DTranMirr (
          		float   *Bout,          /* output interpolated */
                        float   *Bcoeff,        /* input B-spline coeff array */
                        long    Width,          /* width of coeff array */
                        long    Height,         /* height of coeff array */
                        long    Slice,          /* slice of coeff array */
                        double  nx,             /* width of Bout */
                        double  offx,           /* offset of x direction */
                        double  mx,             /* magnif factor for x */
                        double  ny,             /* height of Bout */
                        double  offy,           /* offset of y direction */
                        double  my,             /* magnif factor for y */
                        double  nz,             /* slice of Bout */
                        double  offz,           /* offset of z direction */
                        double  mz,             /* magnif factor for z */
                        long    SplineDegree    /* degree of spline model */
		);

extern void     BatchInterpolatedValue3DTranMirrWarp (
			float   *Bout,          /* output interpolated */
                        float   *Bcoeff,        /* input B-spline coeff array */
                        long    Width,          /* width of coeff array */
                        long    Height,         /* height of coeff array */
                        long    Slice,          /* slice of coeff array */
                        double  nx,             /* width of Bout */
                        double  offx,           /* offset of x direction */
                        double  mx,             /* magnif facor for x */
                        float   *WarpX,         /* deformed value in x direct */
                        double  ny,             /* height of Bout */
                        double  offy,           /* offset of y direction */
                        double  my,             /* magnif facor for y */
                        float   *WarpY,         /* deformed value in y direct */
                        double  nz,             /* slice of Bout */
                        double  offz,           /* offset of z direction */
                        double  mz,             /* magnif facor for z */
                        float   *WarpZ,         /* deformed value in z direct */
                        long    SplineDegree    /* degree of spline model */
                );

extern void BatchInterpolatedGradX3DTranMirr (
          		float   *Bout,          /* output interpolated */
                        float   *Bcoeff,        /* input B-spline coeff array */
                        long    Width,          /* width of coeff array */
                        long    Height,         /* height of coeff array */
                        long    Slice,          /* slice of coeff array */
                        double  nx,             /* width of Bout */
                        double  offx,           /* offset of x direction */
                        double  mx,             /* magnif factor for x */
                        double  ny,             /* height of Bout */
                        double  offy,           /* offset of y direction */
                        double  my,             /* magnif factor for y */
                        double  nz,             /* slice of Bout */
                        double  offz,           /* offset of z direction */
                        double  mz,             /* magnif factor for z */
                        long    SplineDegree    /* degree of spline model */
		);

extern void     BatchInterpolatedGradX3DTranMirrWarp (
			float   *Bout,          /* output interpolated */
                        float   *Bcoeff,        /* input B-spline coeff array */
                        long    Width,          /* width of coeff array */
                        long    Height,         /* height of coeff array */
                        long    Slice,          /* slice of coeff array */
                        double  nx,             /* width of Bout */
                        double  offx,           /* offset of x direction */
                        double  mx,             /* magnif facor for x */
                        float   *WarpX,         /* deformed value in x direct */
                        double  ny,             /* height of Bout */
                        double  offy,           /* offset of y direction */
                        double  my,             /* magnif facor for y */
                        float   *WarpY,         /* deformed value in y direct */
                        double  nz,             /* slice of Bout */
                        double  offz,           /* offset of z direction */
                        double  mz,             /* magnif facor for z */
                        float   *WarpZ,         /* deformed value in z direct */
                        long    SplineDegree    /* degree of spline model */
                );

extern void BatchInterpolatedGradY3DTranMirr (
          		float   *Bout,          /* output interpolated */
                        float   *Bcoeff,        /* input B-spline coeff array */
                        long    Width,          /* width of coeff array */
                        long    Height,         /* height of coeff array */
                        long    Slice,          /* slice of coeff array */
                        double  nx,             /* width of Bout */
                        double  offx,           /* offset of x direction */
                        double  mx,             /* magnif factor for x */
                        double  ny,             /* height of Bout */
                        double  offy,           /* offset of y direction */
                        double  my,             /* magnif factor for y */
                        double  nz,             /* slice of Bout */
                        double  offz,           /* offset of z direction */
                        double  mz,             /* magnif factor for z */
                        long    SplineDegree    /* degree of spline model */
		);

extern void     BatchInterpolatedGradY3DTranMirrWarp (
			float   *Bout,          /* output interpolated */
                        float   *Bcoeff,        /* input B-spline coeff array */
                        long    Width,          /* width of coeff array */
                        long    Height,         /* height of coeff array */
                        long    Slice,          /* slice of coeff array */
                        double  nx,             /* width of Bout */
                        double  offx,           /* offset of x direction */
                        double  mx,             /* magnif facor for x */
                        float   *WarpX,         /* deformed value in x direct */
                        double  ny,             /* height of Bout */
                        double  offy,           /* offset of y direction */
                        double  my,             /* magnif facor for y */
                        float   *WarpY,         /* deformed value in y direct */
                        double  nz,             /* slice of Bout */
                        double  offz,           /* offset of z direction */
                        double  mz,             /* magnif facor for z */
                        float   *WarpZ,         /* deformed value in z direct */
                        long    SplineDegree    /* degree of spline model */
                );

extern void BatchInterpolatedGradZ3DTranMirr (
          		float   *Bout,          /* output interpolated */
                        float   *Bcoeff,        /* input B-spline coeff array */
                        long    Width,          /* width of coeff array */
                        long    Height,         /* height of coeff array */
                        long    Slice,          /* slice of coeff array */
                        double  nx,             /* width of Bout */
                        double  offx,           /* offset of x direction */
                        double  mx,             /* magnif factor for x */
                        double  ny,             /* height of Bout */
                        double  offy,           /* offset of y direction */
                        double  my,             /* magnif factor for y */
                        double  nz,             /* slice of Bout */
                        double  offz,           /* offset of z direction */
                        double  mz,             /* magnif factor for z */
                        long    SplineDegree    /* degree of spline model */
		);

extern void     BatchInterpolatedGradZ3DTranMirrWarp (
			float   *Bout,          /* output interpolated */
                        float   *Bcoeff,        /* input B-spline coeff array */
                        long    Width,          /* width of coeff array */
                        long    Height,         /* height of coeff array */
                        long    Slice,          /* slice of coeff array */
                        double  nx,             /* width of Bout */
                        double  offx,           /* offset of x direction */
                        double  mx,             /* magnif facor for x */
                        float   *WarpX,         /* deformed value in x direct */
                        double  ny,             /* height of Bout */
                        double  offy,           /* offset of y direction */
                        double  my,             /* magnif facor for y */
                        float   *WarpY,         /* deformed value in y direct */
                        double  nz,             /* slice of Bout */
                        double  offz,           /* offset of z direction */
                        double  mz,             /* magnif facor for z */
                        float   *WarpZ,         /* deformed value in z direct */
                        long    SplineDegree    /* degree of spline model */
                );

extern void BatchInterpolatedValue3DTranZero (
          		float   *Bout,          /* output interpolated */
                        float   *Bcoeff,        /* input B-spline coeff array */
                        long    Width,          /* width of coeff array */
                        long    Height,         /* height of coeff array */
                        long    Slice,          /* slice of coeff array */
                        double  nx,             /* width of Bout */
                        double  offx,           /* offset of x direction */
                        double  mx,             /* magnif factor for x */
                        double  ny,             /* height of Bout */
                        double  offy,           /* offset of y direction */
                        double  my,             /* magnif factor for y */
                        double  nz,             /* slice of Bout */
                        double  offz,           /* offset of z direction */
                        double  mz,             /* magnif factor for z */
                        long    SplineDegree    /* degree of spline model */
		);

extern void     BatchInterpolatedValue3DTranZeroWarp (
			float   *Bout,          /* output interpolated */
                        float   *Bcoeff,        /* input B-spline coeff array */
                        long    Width,          /* width of coeff array */
                        long    Height,         /* height of coeff array */
                        long    Slice,          /* slice of coeff array */
                        double  nx,             /* width of Bout */
                        double  offx,           /* offset of x direction */
                        double  mx,             /* magnif facor for x */
                        float   *WarpX,         /* deformed value in x direct */
                        double  ny,             /* height of Bout */
                        double  offy,           /* offset of y direction */
                        double  my,             /* magnif facor for y */
                        float   *WarpY,         /* deformed value in y direct */
                        double  nz,             /* slice of Bout */
                        double  offz,           /* offset of z direction */
                        double  mz,             /* magnif facor for z */
                        float   *WarpZ,         /* deformed value in z direct */
                        long    SplineDegree    /* degree of spline model */
                );

extern void BatchInterpolatedGradX3DTranZero (
          		float   *Bout,          /* output interpolated */
                        float   *Bcoeff,        /* input B-spline coeff array */
                        long    Width,          /* width of coeff array */
                        long    Height,         /* height of coeff array */
                        long    Slice,          /* slice of coeff array */
                        double  nx,             /* width of Bout */
                        double  offx,           /* offset of x direction */
                        double  mx,             /* magnif factor for x */
                        double  ny,             /* height of Bout */
                        double  offy,           /* offset of y direction */
                        double  my,             /* magnif factor for y */
                        double  nz,             /* slice of Bout */
                        double  offz,           /* offset of z direction */
                        double  mz,             /* magnif factor for z */
                        long    SplineDegree    /* degree of spline model */
		);

extern void     BatchInterpolatedGradX3DTranZeroWarp (
			float   *Bout,          /* output interpolated */
                        float   *Bcoeff,        /* input B-spline coeff array */
                        long    Width,          /* width of coeff array */
                        long    Height,         /* height of coeff array */
                        long    Slice,          /* slice of coeff array */
                        double  nx,             /* width of Bout */
                        double  offx,           /* offset of x direction */
                        double  mx,             /* magnif facor for x */
                        float   *WarpX,         /* deformed value in x direct */
                        double  ny,             /* height of Bout */
                        double  offy,           /* offset of y direction */
                        double  my,             /* magnif facor for y */
                        float   *WarpY,         /* deformed value in y direct */
                        double  nz,             /* slice of Bout */
                        double  offz,           /* offset of z direction */
                        double  mz,             /* magnif facor for z */
                        float   *WarpZ,         /* deformed value in z direct */
                        long    SplineDegree    /* degree of spline model */
                );

extern void BatchInterpolatedGradY3DTranZero (
          		float   *Bout,          /* output interpolated */
                        float   *Bcoeff,        /* input B-spline coeff array */
                        long    Width,          /* width of coeff array */
                        long    Height,         /* height of coeff array */
                        long    Slice,          /* slice of coeff array */
                        double  nx,             /* width of Bout */
                        double  offx,           /* offset of x direction */
                        double  mx,             /* magnif factor for x */
                        double  ny,             /* height of Bout */
                        double  offy,           /* offset of y direction */
                        double  my,             /* magnif factor for y */
                        double  nz,             /* slice of Bout */
                        double  offz,           /* offset of z direction */
                        double  mz,             /* magnif factor for z */
                        long    SplineDegree    /* degree of spline model */
		);

extern void     BatchInterpolatedGradY3DTranZeroWarp (
			float   *Bout,          /* output interpolated */
                        float   *Bcoeff,        /* input B-spline coeff array */
                        long    Width,          /* width of coeff array */
                        long    Height,         /* height of coeff array */
                        long    Slice,          /* slice of coeff array */
                        double  nx,             /* width of Bout */
                        double  offx,           /* offset of x direction */
                        double  mx,             /* magnif facor for x */
                        float   *WarpX,         /* deformed value in x direct */
                        double  ny,             /* height of Bout */
                        double  offy,           /* offset of y direction */
                        double  my,             /* magnif facor for y */
                        float   *WarpY,         /* deformed value in y direct */
                        double  nz,             /* slice of Bout */
                        double  offz,           /* offset of z direction */
                        double  mz,             /* magnif facor for z */
                        float   *WarpZ,         /* deformed value in z direct */
                        long    SplineDegree    /* degree of spline model */
                );

extern void BatchInterpolatedGradZ3DTranZero (
          		float   *Bout,          /* output interpolated */
                        float   *Bcoeff,        /* input B-spline coeff array */
                        long    Width,          /* width of coeff array */
                        long    Height,         /* height of coeff array */
                        long    Slice,          /* slice of coeff array */
                        double  nx,             /* width of Bout */
                        double  offx,           /* offset of x direction */
                        double  mx,             /* magnif factor for x */
                        double  ny,             /* height of Bout */
                        double  offy,           /* offset of y direction */
                        double  my,             /* magnif factor for y */
                        double  nz,             /* slice of Bout */
                        double  offz,           /* offset of z direction */
                        double  mz,             /* magnif factor for z */
                        long    SplineDegree    /* degree of spline model */
		);

extern void     BatchInterpolatedGradZ3DTranZeroWarp (
			float   *Bout,          /* output interpolated */
                        float   *Bcoeff,        /* input B-spline coeff array */
                        long    Width,          /* width of coeff array */
                        long    Height,         /* height of coeff array */
                        long    Slice,          /* slice of coeff array */
                        double  nx,             /* width of Bout */
                        double  offx,           /* offset of x direction */
                        double  mx,             /* magnif facor for x */
                        float   *WarpX,         /* deformed value in x direct */
                        double  ny,             /* height of Bout */
                        double  offy,           /* offset of y direction */
                        double  my,             /* magnif facor for y */
                        float   *WarpY,         /* deformed value in y direct */
                        double  nz,             /* slice of Bout */
                        double  offz,           /* offset of z direction */
                        double  mz,             /* magnif facor for z */
                        float   *WarpZ,         /* deformed value in z direct */
                        long    SplineDegree    /* degree of spline model */
                );
