/*******************************************************************************

	Batch interpolation for each location of points
        using point interpolation function, 3D case

        Created by Se Young Chun, Nov 28, 2006, the University of Michigan

*******************************************************************************/

extern void     BatchInterpolatedValue3DMirr (
                        float   *Bout,          /* output interpolated */
                        float   *Bcoeff,        /* input B-spline coeff array */
                        long    Width,          /* width of coeff array */
                        long    Height,         /* height of coeff array */
                        long    Slice,         	/* slice of coeff array */
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

extern void     BatchInterpolatedValue3DMirrWarp (
                        float   *Bout,          /* output interpolated */
                        float   *Bcoeff,        /* input B-spline coeff array */
                        long    Width,          /* width of coeff array */
                        long    Height,         /* height of coeff array */
                        long    Slice,         	/* slice of coeff array */
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

extern void     BatchInterpolatedGradX3DMirr (
                        float   *Bout,          /* output interpolated */
                        float   *Bcoeff,        /* input B-spline coeff array */
                        long    Width,          /* width of the coeff array */
                        long    Height,         /* height of the coeff array */
                        long    Slice,         	/* slice of the coeff array */
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

extern void     BatchInterpolatedGradX3DMirrWarp (
                        float   *Bout,          /* output interpolated */
                        float   *Bcoeff,        /* input B-spline coeff array */
                        long    Width,          /* width of coeff array */
                        long    Height,         /* height of coeff array */
                        long    Slice,         	/* slice of coeff array */
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

extern void     BatchInterpolatedGradY3DMirr (
                        float   *Bout,          /* output interpolated */
                        float   *Bcoeff,        /* input B-spline coeff array */
                        long    Width,          /* width of coeff array */
                        long    Height,         /* height of coeff array */
                        long    Slice,         	/* slice of coeff array */
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

extern void     BatchInterpolatedGradY3DMirrWarp (
                        float   *Bout,          /* output interpolated */
                        float   *Bcoeff,        /* input B-spline coeff array */
                        long    Width,          /* width of coeff array */
                        long    Height,         /* height of coeff array */
                        long    Slice,         	/* slice of coeff array */
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

extern void     BatchInterpolatedGradZ3DMirr (
                        float   *Bout,          /* output interpolated */
                        float   *Bcoeff,        /* input B-spline coeff array */
                        long    Width,          /* width of coeff array */
                        long    Height,         /* height of coeff array */
                        long    Slice,         	/* slice of coeff array */
                        double  nx,             /* width of Bout */
                        double  offx,           /* offset of x direction */
                        double  mx,             /* magnif facor for x */
                        double  ny,             /* height of Bout */
                        double  offy,           /* offset of y direction */
                        double  my,             /* magnif facor for y */
                        double  nz,             /* slice of Bout */
                        double  offz,           /* offset of z direction */
                        double  mz,             /* magnif facor for z */
                        long    SplineDegree    /* degree of spline model */
                );

extern void     BatchInterpolatedGradZ3DMirrWarp (
                        float   *Bout,          /* output interpolated */
                        float   *Bcoeff,        /* input B-spline coeff array */
                        long    Width,          /* width of coeff array */
                        long    Height,         /* height of coeff array */
                        long    Slice,         	/* slice of coeff array */
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

extern void     BatchInterpolatedAll3DMirr (
                        float   *Bo,          	/* output value interpolated */
                        float   *Bogx,          /* output grad x interpolated */
                        float   *Bogy,          /* output grad y interpolated */
                        float   *Bogz,          /* output grad z interpolated */
                        float   *Bcoeff,        /* input B-spline coeff array */
                        long    Width,          /* width of coeff array */
                        long    Height,         /* height of coeff array */
                        long    Slice,         	/* slice of coeff array */
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

extern void     BatchInterpolatedAll3DMirrWarp (
                        float   *Bo,          	/* output value interpolated */
                        float   *Bogx,          /* output grad x interpolated */
                        float   *Bogy,          /* output grad y interpolated */
                        float   *Bogz,          /* output grad z interpolated */
                        float   *Bcoeff,        /* input B-spline coeff array */
                        long    Width,          /* width of coeff array */
                        long    Height,         /* height of coeff array */
                        long    Slice,         	/* slice of coeff array */
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

extern void     BatchInterpolatedThree3DMirr (
                        float   *Bo1,          	/* output1 value interpolated */
                        float   *Bo1gx,         /* output1 gradx interpolated */
                        float   *Bo1gy,         /* output1 grady interpolated */
                        float   *Bo1gz,         /* output1 gradz interpolated */
                        float   *Bo2,         	/* output2 value interpolated */
                        float   *Bo2gx,         /* output2 gradx interpolated */
                        float   *Bo2gy,         /* output2 grady interpolated */
                        float   *Bo2gz,         /* output2 gradz interpolated */
                        float   *Bo3,         	/* output3 value interpolated */
                        float   *Bo3gx,         /* output3 gradx interpolated */
                        float   *Bo3gy,         /* output3 grady interpolated */
                        float   *Bo3gz,         /* output3 gradz interpolated */
                        float   *Bcoeff1,       /* input Bspline coeff1 array */
                        float   *Bcoeff2,       /* input Bspline coeff2 array */
                        float   *Bcoeff3,       /* input Bspline coeff3 array */
                        long    Width,          /* width of coeff array */
                        long    Height,         /* height of coeff array */
                        long    Slice,         	/* slice of coeff array */
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

extern void     BatchInterpolatedValue3DZero (
                        float   *Bout,          /* output interpolated */
                        float   *Bcoeff,        /* input B-spline coeff array */
                        long    Width,          /* width of coeff array */
                        long    Height,         /* height of coeff array */
                        long    Slice,         	/* slice of coeff array */
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

extern void     BatchInterpolatedValue3DZeroWarp (
                        float   *Bout,          /* output interpolated */
                        float   *Bcoeff,        /* input B-spline coeff array */
                        long    Width,          /* width of coeff array */
                        long    Height,         /* height of coeff array */
                        long    Slice,         	/* slice of coeff array */
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

extern void     BatchInterpolatedGradX3DZero (
                        float   *Bout,          /* output interpolated */
                        float   *Bcoeff,        /* input B-spline coeff array */
                        long    Width,          /* width of the coeff array */
                        long    Height,         /* height of the coeff array */
                        long    Slice,         	/* slice of the coeff array */
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

extern void     BatchInterpolatedGradX3DZeroWarp (
                        float   *Bout,          /* output interpolated */
                        float   *Bcoeff,        /* input B-spline coeff array */
                        long    Width,          /* width of coeff array */
                        long    Height,         /* height of coeff array */
                        long    Slice,         	/* slice of coeff array */
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

extern void     BatchInterpolatedGradY3DZero (
                        float   *Bout,          /* output interpolated */
                        float   *Bcoeff,        /* input B-spline coeff array */
                        long    Width,          /* width of coeff array */
                        long    Height,         /* height of coeff array */
                        long    Slice,         	/* slice of coeff array */
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

extern void     BatchInterpolatedGradY3DZeroWarp (
                        float   *Bout,          /* output interpolated */
                        float   *Bcoeff,        /* input B-spline coeff array */
                        long    Width,          /* width of coeff array */
                        long    Height,         /* height of coeff array */
                        long    Slice,         	/* slice of coeff array */
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

extern void     BatchInterpolatedGradZ3DZero (
                        float   *Bout,          /* output interpolated */
                        float   *Bcoeff,        /* input B-spline coeff array */
                        long    Width,          /* width of coeff array */
                        long    Height,         /* height of coeff array */
                        long    Slice,         	/* slice of coeff array */
                        double  nx,             /* width of Bout */
                        double  offx,           /* offset of x direction */
                        double  mx,             /* magnif facor for x */
                        double  ny,             /* height of Bout */
                        double  offy,           /* offset of y direction */
                        double  my,             /* magnif facor for y */
                        double  nz,             /* slice of Bout */
                        double  offz,           /* offset of z direction */
                        double  mz,             /* magnif facor for z */
                        long    SplineDegree    /* degree of spline model */
                );

extern void     BatchInterpolatedGradZ3DZeroWarp (
                        float   *Bout,          /* output interpolated */
                        float   *Bcoeff,        /* input B-spline coeff array */
                        long    Width,          /* width of coeff array */
                        long    Height,         /* height of coeff array */
                        long    Slice,         	/* slice of coeff array */
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

extern void     BatchInterpolatedAll3DZero (
                        float   *Bo,          	/* output value interpolated */
                        float   *Bogx,          /* output grad x interpolated */
                        float   *Bogy,          /* output grad y interpolated */
                        float   *Bogz,          /* output grad z interpolated */
                        float   *Bcoeff,        /* input B-spline coeff array */
                        long    Width,          /* width of coeff array */
                        long    Height,         /* height of coeff array */
                        long    Slice,         	/* slice of coeff array */
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

extern void     BatchInterpolatedAll3DZeroWarp (
                        float   *Bo,          	/* output value interpolated */
                        float   *Bogx,          /* output grad x interpolated */
                        float   *Bogy,          /* output grad y interpolated */
                        float   *Bogz,          /* output grad z interpolated */
                        float   *Bcoeff,        /* input B-spline coeff array */
                        long    Width,          /* width of coeff array */
                        long    Height,         /* height of coeff array */
                        long    Slice,         	/* slice of coeff array */
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
                );
